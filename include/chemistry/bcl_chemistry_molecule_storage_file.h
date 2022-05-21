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

#ifndef BCL_CHEMISTRY_MOLECULE_STORAGE_FILE_H_
#define BCL_CHEMISTRY_MOLECULE_STORAGE_FILE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "io/bcl_io_retrieve_interface.h"
#include "io/bcl_io_store_interface.h"
#include "sched/bcl_sched_mutex.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_wrapper_enum.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeStorageFile
    //! @brief is a handler for storing and retrieving smallmolecule files from a file system.
    //!
    //! @details SmallMolecules are stored in a file structure that follows that logic:
    //!        {path/{SOURCE}filename.sdf} where the filename can have a prefix component {SOURCE}
    //!
    //! SOURCE is a prefix which can be empty
    //!
    //! In an file the key corresponds to the index of molecules which is auto incremented and is only numeric
    //! when initialized as attached, the largest key is located
    //! when initialized as created or overwrite, the key starts with 0
    //!
    //! @see @link example_chemistry_molecule_storage_file.cpp @endlink
    //! @author mendenjl
    //! @date Mar 09, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeStorageFile :
      public io::StoreInterface< ConformationInterface>,
      public io::RetrieveInterface< MoleculeComplete, MoleculeEnsemble>
    {

    protected:

    //////////
    // data //
    //////////

      //! filename of file that contains storage
      std::string m_Filename;

      //! whether the enumerated instance is intended for storage (true) or retrieval (false)
      bool        m_Store;

      //! number of unique conformations in a given file
      mutable size_t      m_EnsembleSize;

      //! number of unique constitutions in a given file
      mutable size_t      m_ConstitutionsCount;

      //! map of molecule sha1 string to index (contains both constitution and conformation strings and respective indices as values)
      mutable util::SiPtr< storage::Vector< storage::Map< std::string, size_t> > > m_MoleculeHashIndexMap;

      //! file size
      mutable long m_FileSize;

      //! vector of molecule sizes; built the 1st time that GetSize is called
      mutable storage::Vector< size_t> m_MoleculeSizes;

      //! mutex that must be locked to make changes to MoleculeHashIndexMap
      mutable util::SiPtr< sched::Mutex> m_HashCacheMutex;

      //! Map of known hash caches, indexed by filename
      static storage::Map
      <
        std::string,
        storage::Pair
        <
          sched::Mutex,
          storage::Vector< storage::Map< std::string, size_t> >
        >
      > s_HashCache;

    public:

    //////////
    // data //
    //////////

      //! store and retrieve instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_StoreInstance;
      static const util::SiPtr< const util::ObjectInterface> s_RetrieveInstance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from bool; which indicates whether the enumerated instance is intended for storage or retrieval
      explicit MoleculeStorageFile( const bool &STORE = false);

      //! @brief Clone function
      //! @return pointer to new MoleculeStorageFile
      MoleculeStorageFile *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

      //! @brief get the filename
      //! @return the filename
      const std::string &GetFilename() const
      {
        return m_Filename;
      }

      //! @brief number of molecules in storage
      //! @return number of molecules in this storage
      size_t GetSize() const;

      //! @brief get all keys for this storage
      //! @return all keys for this storage
      storage::Vector< std::string> GetAllKeys() const;

      //! @brief get the # of atoms for the molecule of the specified key (units defined by t_Type)
      //! @param KEY the identifier for the specific object
      //! @return the number of t_Type-defined units possessed by the key
      size_t GetKeySize( const std::string &KEY) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief get stored MoleculeComplete from key
      //! @param KEY key identifier for specific object
      //! @return Molecule of interest
      MoleculeComplete Retrieve( const std::string &KEY) const;

      //! @brief get molecule ensemble for given keys
      //! @param KEYS vector of keys
      //! @return list of MoleculeComplete objects corresponding to KEYS
      MoleculeEnsemble RetrieveEnsemble( const storage::Vector< std::string> &KEYS) const;

      //! @brief get ensemble of stored molecules objects for a range of keys, specified by plain numbers
      //! @param RANGE range of keys
      //! @return list of MoleculeComplete objects corresponding to keys in RANGE
      MoleculeEnsemble RetrieveEnsemble( const math::Range< size_t> &RANGE) const;

      //! @brief store molecule
      //! @param MOLECULE molecule to store
      //! @return key associated with molecule
      std::string Store( const ConformationInterface &MOLECULE);

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      void CreateConstitutionConformationHashMaps( const ConformationInterface &MOLECULE) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief read from std::ostream
      //! @param OSTREAM input stream
      //! @param INDENT indentation
      //! @return ostream which was read from
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! Enum to indicate what find is searching for (constitution, configuration, conformation)
      enum FindType
      {
        e_Constitution,
        e_Conformation,
        s_NumberOfFindTypes
      };

      //! @brief Alignment as string
      //! @param ALIGNMENT the alignment
      //! @return the Alignment as string
      static const std::string &GetFindTypeString( const FindType &TYPE);

      //! SSEInfoTypeEnum simplifies the usage of the SSEInfoType enum of this class
      typedef util::WrapperEnum< FindType, &GetFindTypeString, s_NumberOfFindTypes> FindTypeEnum;

      //! @brief find a molecule in the file; must have exactly matching conformation, but misc properties are ignored
      //! @param MOLECULE the molecule to find
      //! @param TYPE whether the index is that of the first conformation or constitution
      //! @return pair; bool indicates whether the molecule is already in the file, size_t indicates position in the file,
      //! or, if the molecule is not in the file, the size of the file
      std::pair< size_t, bool> Find( const ConformationInterface &MOLECULE, const FindType &TYPE = e_Conformation) const;

      void CacheHashes( const FindType &TYPE) const;
    }; // class MoleculeStorageFile

  } // namespace chemistry

} // namespace bcl

#endif // BCL_CHEMISTRY_MOLECULE_STORAGE_FILE_H_
