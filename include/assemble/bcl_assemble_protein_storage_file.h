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

#ifndef BCL_ASSEMBLE_PROTEIN_STORAGE_FILE_H_
#define BCL_ASSEMBLE_PROTEIN_STORAGE_FILE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "io/bcl_io_directory.h"
#include "io/bcl_io_stream_buffer_classes.h"
#include "io/bcl_io_stream_buffer_interface.h"
#include "pdb/bcl_pdb_factory.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinStorageFile
    //! @brief class is a handler for storing to and retrieving files from the file system
    //! @details Proteins are stored in a directory structure that follows that logic:
    //! {INITIALIZER}/{SOURCE}{KEY}.pdb
    //!
    //! KEY is always 6 chars long "000001", "000002"
    //!
    //! the key is auto incremented an is only numeric
    //! when initialized as attached, the largest key is located
    //! when initialized as created or overwrite, the key starts with 0
    //!
    //! @see @link example_assemble_protein_storage_file.cpp @endlink
    //! @author woetzen
    //! @date Jul 30, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinStorageFile :
      public util::ObjectInterface
    {
    public:

      //! @enum InitializerType
      //! enumerator for initialization types
      enum InitializerType
      {
        e_Attach,    //!< Attach to storage (files not over-written, instead, appends new entries)
        e_Overwrite, //!< overwrite forces creation for price of deletion
        s_MaxInitializerType //!< max number of initializer types
      };

      //! @brief InitializerType as string
      //! @param INIT_TYPE the InitializerType
      //! @return the string for the INIT_TYPE
      static const std::string &GetInitializerTypeDescriptor( const InitializerType &INIT_TYPE);

      //! @brief TypeEnum is used for I/O of InitializerType
      typedef util::WrapperEnum< InitializerType, &GetInitializerTypeDescriptor, s_MaxInitializerType> TypeEnum;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief default storage as given in the storage flag used in the command line
      static util::ShPtr< ProteinStorageFile> GetDefaultStorage();

      //! @brief flag for default ProteinStorage
      //! @return ShPtr to flag that is used for default protein storage
      static const util::ShPtr< command::FlagInterface> &GetDefaultStorageFlag();

    private:

    //////////
    // data //
    //////////

      //! directory e.g. /home/user/run5/
      util::ShPtr< io::Directory> m_Directory;

      //! initializer type
      InitializerType m_InitializerType;

      //! PDB factory
      util::ShPtr< pdb::Factory> m_Factory;

      //! set of sources that have been used to store
      storage::Set< std::string> m_Sources;

      //! Default compression type
      io::StreamBufferClass m_CompressionType;

    public:

    //////////
    // data //
    //////////

      //! large key size (has stage)
      static const size_t s_LargeKeySize = 11;

      //! small key size (no stage)
      static const size_t s_SmallKeySize = 5;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief initialize with stream buffer class type
      ProteinStorageFile
      (
        const std::string &INITIALIZER = std::string(),
        const InitializerType &INIT_TYPE = e_Attach,
        const io::StreamBufferClass &COMPRESSION = io::StreamBufferClass()
      );

      //! @brief Clone function
      //! @return pointer to new ProteinStorageFile
      ProteinStorageFile *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief sets the pdb factory
      //! @param FACTORY PDB factory to use
      void SetFactory( const util::ShPtr< pdb::Factory> &FACTORY);

      //! @brief initialize the protein storage
      //! @param INITIALIZER path (including prefix) in which proteins are stored
      //! @param INIT_FLAG flag for type of initialization
      //! @return true if initialize was successful
      bool SetDirectory( const std::string &INITIALIZER);

      //! @brief get the initializer
      //! @return string used to initialize protein storage
      const std::string &GetInitializer() const;

      //! @brief get the initialization type
      //! @return true if existing files in the storage will be overwritten
      bool WillOverwrite() const;

      //! @brief number of proteins in source
      //! @param SOURCE corresponds to the prefix for the files e.g "run1"
      //! @return number of proteins with that source prefix
      size_t GetSize( const std::string &SOURCE) const;

      //! @brief gets all sources in this storage
      //! @return all sources in this storage
      storage::Set< std::string> GetAllSources() const;

      //! @brief get all keys for given source
      //! @param SOURCE prefix like "run1"
      //! @return all keys of given prefix
      storage::Vector< std::string> GetAllKeys( const std::string &SOURCE) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief get protein
      //! @param SOURCE prefix of protein eg. "run1"
      //! @param KEY key identifier for specific protein in given source like "00001"
      //! @return shptr to protein of interest, undefined if there is no such protein
      util::ShPtr< ProteinModel> Retrieve( const std::string &SOURCE, const std::string &KEY) const;

      //! @brief get ensemble of proteins
      //! @param SOURCE corresponds to the prefix for the files e.g "run1"
      //! @return shptr list of proteins from given source
      util::ShPtrList< ProteinModel> RetrieveEnsemble( const std::string &SOURCE) const;

      //! @brief get ensemble of proteins for given keys
      //! @param SOURCE corresponds to the prefix for the files e.g "run1"
      //! @param KEYS vector of identifiers for specific proteins in given source like "00001", "00002"
      //! @return shptr list of proteins from given source
      util::ShPtrList< ProteinModel> RetrieveEnsemble
      (
        const std::string &SOURCE,
        const storage::Vector< std::string> &KEYS
      ) const;

      //! @brief get ensemble of proteins for given keys which does not map to keys, following plain numbering
      //! @param SOURCE corresponds to the prefix for the files e.g "run1"
      //! @param RANGE range of proteins
      //! @return shptr list of proteins from given source
      util::ShPtrList< ProteinModel> RetrieveEnsemble
      (
        const std::string &SOURCE,
        const math::Range< size_t> &RANGE
      ) const;

      //! @brief get any PDB remark lines
      //! @param SOURCE prefix of protein eg. "run1"
      //! @param KEY key identifier for specific protein in given source like "00001"
      //! @return PDB remark lines
      util::ShPtrList< pdb::Line> RetrieveRemarkLines
      (
        const std::string &SOURCE,
        const std::string &KEY
      ) const;

      //! @brief store protein
      //! @param MOLECULE protein to store
      //! @param SOURCE corresponds to the prefix for the files e.g "run1"
      //! @return key associated with protein
      std::string Store( const ProteinModel &MOLECULE, const std::string &SOURCE);

      //! @brief store protein
      //! @param MOLECULE protein to store
      //! @param SOURCE corresponds to the prefix for the files e.g "run1"
      //! @param KEY key under which protein was stored
      //! @return true if store was successful
      bool Store( const ProteinModel &MOLECULE, const std::string &SOURCE, const std::string &KEY);

      //! @brief store ensemble of proteins
      //! @param ENSEMBLE shptr list of proteins
      //! @param SOURCE corresponds to the prefix for the files e.g "run1"
      //! @return vector of keys
      storage::Vector< std::string> Store
      (
        const util::ShPtrList< ProteinModel> &ENSEMBLE,
        const std::string &SOURCE
      );

      //! @brief deletes all entries in this storage
      void Delete();

    ///////////////
    // operators //
    ///////////////

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

    public:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief construct complete filename from source and key
      //! @brief SOURCE the source for the model
      //! @brief KEY the key for that protein
      //! @brief filename of form {initializer}/{SOURCE}{KEY}.pdb
      std::string Filename( const std::string &SOURCE, const std::string &KEY) const;

      //! @brief converts a key into a string
      //! @param KEY key to be converted
      //! @return converted string
      static std::string KeyToString( const size_t &KEY);

      //! @brief check if key is valid string
      //! @param KEY the key to be checked
      //! @return true if the key is valid
      static bool IsValidKey( const std::string &KEY);

    private:

      //! @brief gets all entries that have the given source
      //! @param SOURCE source of entries
      //! @return list of found entries
      storage::List< io::DirectoryEntry> GetAllEntries( const std::string &SOURCE) const;

      //! @brief removes all entries that have the given source
      //! @param SOURCE source of entries to be removed
      void RemoveEntries( const std::string &SOURCE) const;

    }; // class ProteinStorageFile

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_PROTEIN_STORAGE_FILE_H_
