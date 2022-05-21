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

#ifndef BCL_ASSEMBLE_PROTEIN_WITH_CACHE_STORAGE_FILE_H_
#define BCL_ASSEMBLE_PROTEIN_WITH_CACHE_STORAGE_FILE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "pdb/bcl_pdb.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_assemble_protein_model_with_cache.h"
#include "io/bcl_io_retrieve_interface.h"
#include "io/bcl_io_store_interface.h"
#include "pdb/bcl_pdb_residue_simple.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinWithCacheStorageFile
    //! @brief is a handler for retrieving pdb or fasta files from a file system.
    //!
    //! @details Proteins are stored in a file structure like {path/{pdbid}{suffix}
    //!        pdb ids may also be explicitly given in a list file in which case
    //!        pdbs and/or fastas matching the given pdb id (and optional suffix)
    //!        will be searched for starting from the base directory, and going up to 3 directories below it
    //!
    //! @see @link example_assemble_protein_with_cache_storage_file.cpp @endlink
    //! @author mendenjl
    //! @date Jan 09, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinWithCacheStorageFile :
      public RetrieveProteinModelWithCache
    {

    protected:

    //////////
    // data //
    //////////

      //! directory where proteins are stored
      std::string m_Directory;

      //! optional suffix to look for
      std::string m_Suffix;

      //! optional pdb id list file
      std::string m_PdbIdListFile;

      //! optional pdb path list file
      std::string m_PdbPathFile;

      //! All protein files in the directory
      storage::Vector< std::string> m_ProteinFiles;

      //! # of AAs for each of the protein files, if known
      mutable storage::Vector< size_t> m_ProteinFileSizes;

      //! Map from key (e.g. 1UBIA) to path index in m_ProteinFiles
      storage::Map< std::string, size_t> m_KeyToPathId;

      //! Class of AAs to produce
      biol::AAClass m_Class;

      //! whether to only retrieve AAs that have defined coordinates (implies that only pdb files will be considered)
      bool m_RequireCoordinates;

      //! Whether to auto-load blast and any available SSPred files
      bool m_Autoload;

      //! Whether to automatically perform DSSP on every protein loaded in (applies if m_RequireCoordinates is true)
      bool m_DoDSSP;

      //! alias, based on whether to exclude loops or not and the AAClass
      std::string m_Alias;

      //! thickness of transition region
      double m_TransitionRegionThickness;

      //! thickness of gap region
      double m_GapThickness;

      //! blast extension; if given, an assert will be thrown if no blast profile with the given extension is found
      std::string m_BlastExtension;

      //! Bool; if set, autodetect all proteins in any directory or subdirectory (default is to look only in the given directory)
      bool m_RecursiveFind;

      //! last protein model loaded via Retrieve
      mutable util::ShPtr< ProteinModelWithCache> m_LastModel;

      //! key for last model loaded
      mutable std::string                         m_LastModelKey;

      //! true if the index file has been loaded or attempted to be loaded
      mutable bool                                m_HaveLoadedIndexFile;

      //! true if the index file has been loaded or attempted to be loaded
      mutable bool                                m_CanWriteIndexFile;

    public:

    //////////
    // data //
    //////////

      //! retrieve instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_SequenceInstance;
      static const util::SiPtr< const util::ObjectInterface> s_ProteinInstance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor, takes whether or not to require coordinates
      //! @param REQUIRE_COORDINATES if true, only load PDB files, and only consider AAs with coordinates
      //! @param AUTOLOAD if true, automatically load any blast and secondary structure prediction files
      ProteinWithCacheStorageFile
      (
        const bool &REQUIRE_COORDINATES,
        const bool &AUTOLOAD = true
      );

      //! @brief copy constructor; ignores all cached values
      ProteinWithCacheStorageFile( const ProteinWithCacheStorageFile &PARENT);

      //! @brief Clone function
      //! @return pointer to new ProteinWithCacheStorageFile
      ProteinWithCacheStorageFile *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

      //! @brief get the directory name
      //! @return the directory name
      const std::string &GetDirectory() const
      {
        return m_Directory;
      }

      //! @brief number of proteins in storage
      //! @return number of proteins in this storage
      size_t GetSize() const;

      //! @brief get all keys for this storage
      //! @return all keys for this storage
      storage::Vector< std::string> GetAllKeys() const;

      //! @brief get the full path for all keys
      //! @return get the full path for all keys
      const storage::Vector< std::string> &GetAllLocations() const
      {
        return m_ProteinFiles;
      }

      //! @brief get the type of AA class that will be produced
      //! @return the type of AA class that will be produced
      const biol::AAClass &GetAAClass() const;

      //! @brief count the number of AAs in a given key, without actually creating a complete model
      //! @param KEY the identifier for the protein of interest
      //! @return the number of AAs (or AAs with coordinates if m_RequireCoordinates is set) for the given key
      size_t GetKeySize( const std::string &KEY) const;

      //! @brief return true if this class works with proteins (and thus requires coordinates), false if we are just working
      //!        with sequences
      bool GetRequireCoordinates() const
      {
        return m_RequireCoordinates;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief get the filename for a given key
      //! @param KEY key identifier for specific object
      //! @return the filename for the given key (empty if the key is invalid for this storage)
      const std::string &GetFilename( const std::string &KEY) const;

      //! @brief get stored ProteinModelWithCache from key
      //! @param KEY key identifier for specific object
      //! @return Protein of interest
      util::ShPtr< ProteinModelWithCache> Retrieve( const std::string &KEY) const;

      //! @brief get molecule ensemble for given keys
      //! @param KEYS vector of keys
      //! @return list of ProteinModelWithCache objects corresponding to KEYS
      util::ShPtrVector< ProteinModelWithCache> RetrieveEnsemble( const storage::Vector< std::string> &KEYS) const;

      //! @brief get ensemble of stored proteins objects for a range of keys, specified by plain numbers
      //! @param RANGE range of keys
      //! @return list of ProteinModelWithCache objects corresponding to keys in RANGE
      util::ShPtrVector< ProteinModelWithCache> RetrieveEnsemble( const math::Range< size_t> &RANGE) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief given a fasta or pdb filename, generate a protein model
      //! @param FILENAME filename of interest
      //! @param REQUIRE_COORDINATES true if coordinates are required
      static ProteinModel GenerateProteinModelFromFile
      (
        const std::string &FILENAME,
        const bool &REQUIRE_COORDINATES,
        const biol::AAClass &CLASS
      );

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

    protected:

      //! @brief add protein keys and paths from the given directory
      //! @param DIRECTORY the directory of interest
      //! @param ERR_STREAM output stream for errors
      //! @return true if no duplicate keys were found
      bool AddKeys( const io::Directory &DIRECTORY, std::ostream &ERR_STREAM);

      //! @brief add protein keys and paths from the given directory if they are in the given set, remove them from the set as each is added
      //! @param DIRECTORY the directory of interest
      //! @param DESIRED_KEYS set of keys of interest
      //! @param ERR_STREAM output stream for errors
      //! @return false if any duplicates are found
      bool AddKeys( const io::Directory &DIRECTORY, storage::Set< std::string> &DESIRED_KEYS, std::ostream &ERR_STREAM);

      //! @brief load all blast and secondary structure prediction information for the model
      //! @param MODEL protein model of interest
      //! @param PREFIX filename of the model, without the pdb/fasta extension
      void LoadBlastSspred( ProteinModel &MODEL, const std::string &PREFIX) const;

      //! @brief helper function to get the key from a filename
      //! @param FILENAME filename of interest
      //! @param SUFFIX suffix that can be removed from the filename
      //! @return key
      static std::string GetKey( const std::string &FILENAME, const std::string &SUFFIX);

      //! @brief helper function to determine the total number of residues in a map returned from pdb::Head
      //! @param MAP map, usually obtained from pdb::Head::GetSEQRESProteinChains or pdb::Head::GetMissingResidues
      //! @return the total number of residues in the map
      static size_t GetMapSize( const storage::Map< char, storage::List< pdb::ResidueSimple> > &MAP);

      //! @brief helper function to count the # of residues in a given PDB file from the atom lines
      //! @param ISTREAM input file
      static size_t CountResiduesInAtomLines( std::istream &INPUT);

    private:

      //! undefined assignment operator
      ProteinWithCacheStorageFile &operator=( const ProteinWithCacheStorageFile &);

      //! @brief get the name of the index file for this protein cache storage file
      std::string GetIndexFileName() const;

    }; // class ProteinWithCacheStorageFile

  } // namespace assemble

} // namespace bcl

#endif // BCL_ASSEMBLE_PROTEIN_WITH_CACHE_STORAGE_FILE_H_
