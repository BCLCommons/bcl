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

#ifndef BCL_CONTACT_CORRELATION_STORAGE_FILE_H_
#define BCL_CONTACT_CORRELATION_STORAGE_FILE_H_

// include the namespace header
#include "bcl_contact.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "align/bcl_align_alignment_node.h"
#include "io/bcl_io_directory.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CorrelationStorageFile
    //! @brief is a handler for storing/retrieving target sequences and their corresponding MSA and correlation matrix
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_contact_correlation_storage_file.cpp @endlink
    //! @author teixeipl
    //! @date Aug 29, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! CorrelationStorageFile should inherit from a generalized version of InterfaceStorageInterface but it currently
    //! is specific for model storage
    class BCL_API CorrelationStorageFile :
      public util::ObjectInterface
    {
    //////////
    // data //
    //////////

    public:

      typedef storage::Pair< align::AlignmentNode< biol::AABase>, CorrelationMatrix> AlignmentMatrixPair;

      //! @enum InitializerType
      //! @brief enumerator for initialization types
      enum InitializerType
      {
        e_Create,    //!< create storage, if exists initialization should fail
        e_Attach,    //!< attach to storage, if does not exist initialization should fail
        e_Overwrite, //!< overwrite forces creation for price of deletion
        s_MaxInitializerType //!< max number of initializer types
      };

      static const util::Format s_KeyToString; //!< format to convert Key to string
      static const std::string s_FileExtension; //!< default file extension for correlation storage files
      static const std::string s_FilePrefix; //!< prefix of AlignmentMatrixPair files stored
      static const util::SiPtr< const util::ObjectInterface> s_Instance; //!< single instance of that class
      static const char s_Delim; //!< Default delimiter to be used for reading/writing alignments to file

    private:

      util::ShPtr< io::Directory> m_Directory; //!< directory e.g. /home/user/...

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:
      //! @brief Clone function
      //! @return pointer to new CorrelationStorageFile
      CorrelationStorageFile *Clone() const
      {
        return new CorrelationStorageFile( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get member variable directory
      //! @return shptr on current directory
      util::ShPtr< io::Directory> GetDirectory() const
      {
        return m_Directory;
      }

      //! @brief initialize the correlation storage
      //! @param INITIALIZER encodes where data is stored
      //! @param INIT_FLAG flag for type of initialization
      //! @return true if initialize was successful
      bool Initialize( const std::string &INITIALIZER, const InitializerType INIT_FLAG);

      //! @brief get the initializer
      //! @return string used to initialize object storage
      const std::string &GetInitializer() const;

      //! @brief number of data items in source
      //! @param SOURCE source of data file
      //! @return number of MSA correlation matrix pairs in storage at SOURCE
      size_t GetSize() const;

      //! @brief get all keys for given source
      //! @return all keys of given source
      storage::Vector< std::string> GetAllKeys() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief get AlignmentMatrixPair
      //! @param KEY key identifier for specific molecule in given source
      //! @return shptr to AlignmentMatrixPair of interest
      util::ShPtr< AlignmentMatrixPair> Retrieve( const std::string &KEY) const;

      //! @brief get entire ensemble of AlignmentMatrixPairs in directory
      //! @return shptr list of AlignmentMatrixPair from given source
      util::ShPtrList< AlignmentMatrixPair> RetrieveEnsemble() const;

      //! @brief get ensemble of AlignmentMatrixPair for given keys
      //! @param KEYS vector of keys
      //! @return shptr list of molecules from given source
      util::ShPtrList< AlignmentMatrixPair> RetrieveEnsemble( const storage::Vector< std::string> &KEYS) const;

//      //! @brief get ensemble of AlignmentMatrixPair for given keys which does not map to keys, following plain numbering
//      //! @param RANGE range of AlignmentMatrixPair
//      //! @return shptr list of AlignmentMatrixPair from given source
//      util::ShPtrList< AlignmentMatrixPair> RetrieveEnsemble( const math::Range< size_t> &RANGE) const;

      //! @brief store AlignmentMatrixPairs
      //! @param ALIGNMENT_MATRIX_PAIR to store
      //! @return key associated with molecule
      std::string Store( const AlignmentMatrixPair &ALIGNMENT_MATRIX_PAIR);

      //! @brief store AlignmentMatrixPairs
      //! @param ALIGNMENT_MATRIX_PAIR to store
      //! @param KEY key under which AlignmentMatrixPair was stored
      //! @return true if store was successful
      bool Store( const AlignmentMatrixPair &ALIGNMENT_MATRIX_PAIR, const std::string &KEY);

      //! @brief store ensemble of AlignmentMatrixPair
      //! @param ENSEMBLE shptr list of AlignmentMatrixPair
      //! @return vector of keys
      storage::Vector< std::string> Store
      (
        const util::ShPtrList< AlignmentMatrixPair> &ENSEMBLE
      );

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

//      //! @brief check if key is valid string
//      //! @param KEY the key to be checked
//      //! @return true if the key is valid
//      bool IsValidKey( const std::string &KEY) const
//      {
//        // check that the key is of size 6 and that it is an unsigned integer
//        return KEY.length() == 6 && util::LengthOfUnsignedIntegerType( KEY) == 6;
//      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief construct complete filename from source and key
      //! @brief KEY the key for that protein
      //! @brief filename of form {initializer}/s_FilePrefix{KEY}.s_FileExtension
      std::string StorageFilename( const std::string &KEY) const
      {
        return m_Directory->AppendFilename( KEY + s_FilePrefix + io::File::GetExtensionDelimiter() + s_FileExtension);
      }

    }; // class CorrelationStorageFile

  } // namespace contact
} // namespace bcl

#endif // BCL_CONTACT_CORRELATION_STORAGE_FILE_H_
