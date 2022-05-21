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
#include "assemble/bcl_assemble_protein_with_cache_storage_file.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_blast_profile_handler.h"
#include "biol/bcl_biol_dssp.h"
#include "biol/bcl_biol_membrane.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_head.h"
#include "sspred/bcl_sspred_mahssmi.h"
#include "sspred/bcl_sspred_method_handler.h"
#include "sspred/bcl_sspred_pdb.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    const util::SiPtr< const util::ObjectInterface> ProteinWithCacheStorageFile::s_SequenceInstance
    (
      util::Enumerated< RetrieveProteinModelWithCache>::AddInstance( new ProteinWithCacheStorageFile( false))
    );
    const util::SiPtr< const util::ObjectInterface> ProteinWithCacheStorageFile::s_ProteinInstance
    (
      util::Enumerated< RetrieveProteinModelWithCache>::AddInstance( new ProteinWithCacheStorageFile( true))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor, takes whether or not to require coordinates
    //! @param REQUIRE_COORDINATES if true, only load PDB files, and only consider AAs with coordinates
    //! @param AUTOLOAD if true, automatically load any blast and secondary structure prediction files
    ProteinWithCacheStorageFile::ProteinWithCacheStorageFile
    (
      const bool &REQUIRE_COORDINATES,
      const bool &AUTOLOAD
    ) :
      m_Class( REQUIRE_COORDINATES ? biol::GetAAClasses().e_AABackBone : biol::GetAAClasses().e_AA),
      m_RequireCoordinates( REQUIRE_COORDINATES),
      m_Autoload( AUTOLOAD),
      m_DoDSSP( false),
      m_Alias( std::string( m_RequireCoordinates ? "Protein" : "Sequence") + "Directory"),
      m_TransitionRegionThickness( 5.0),
      m_GapThickness( 0.0),
      m_BlastExtension(),
      m_RecursiveFind( false),
      m_HaveLoadedIndexFile( false),
      m_CanWriteIndexFile( true)
    {
    }

    //! @brief copy constructor; ignores all cached values
    ProteinWithCacheStorageFile::ProteinWithCacheStorageFile( const ProteinWithCacheStorageFile &PARENT) :
      m_Directory( PARENT.m_Directory),
      m_Suffix( PARENT.m_Suffix),
      m_PdbIdListFile( PARENT.m_PdbIdListFile),
      m_PdbPathFile( PARENT.m_PdbPathFile),
      m_ProteinFiles( PARENT.m_ProteinFiles),
      m_ProteinFileSizes( PARENT.m_ProteinFileSizes),
      m_KeyToPathId( PARENT.m_KeyToPathId),
      m_Class( PARENT.m_Class),
      m_RequireCoordinates( PARENT.m_RequireCoordinates),
      m_Autoload( PARENT.m_Autoload),
      m_DoDSSP( PARENT.m_DoDSSP),
      m_Alias( PARENT.m_Alias),
      m_TransitionRegionThickness( PARENT.m_TransitionRegionThickness),
      m_GapThickness( PARENT.m_GapThickness),
      m_BlastExtension( PARENT.m_BlastExtension),
      m_RecursiveFind( PARENT.m_RecursiveFind),
      m_HaveLoadedIndexFile( PARENT.m_HaveLoadedIndexFile),
      m_CanWriteIndexFile( PARENT.m_CanWriteIndexFile)
    {
    }

    //! @brief Clone function
    //! @return pointer to new ProteinWithCacheStorageFile
    ProteinWithCacheStorageFile *ProteinWithCacheStorageFile::Clone() const
    {
      return new ProteinWithCacheStorageFile( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinWithCacheStorageFile::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &ProteinWithCacheStorageFile::GetAlias() const
    {
      return m_Alias;
    }

    //! @brief number of proteins in storage
    //! @return number of proteins in this storage
    size_t ProteinWithCacheStorageFile::GetSize() const
    {
      return m_ProteinFiles.GetSize();
    }

    //! @brief get all keys for given source
    //! @return all keys of given prefix
    storage::Vector< std::string> ProteinWithCacheStorageFile::GetAllKeys() const
    {
      // make a vector of strings big enough to hold the keys
      storage::Vector< std::string> keys( GetSize());

      for( size_t key( 0), max_key( keys.GetSize()); key < max_key; ++key)
      {
        keys( key) = GetKey( m_ProteinFiles( key), m_Suffix);
      }

      // return all created keys
      return keys;
    }

    //! @brief count the number of AAs in a given key, without actually creating a complete model
    //! @param KEY the identifier for the protein of interest
    //! @return the number of AAs (or AAs with coordinates if m_RequireCoordinates is set) for the given key
    size_t ProteinWithCacheStorageFile::GetKeySize( const std::string &KEY) const
    {
      // find the filename id for the given key
      const std::string cleaned_key( GetKey( KEY, m_Suffix));
      storage::Map< std::string, size_t>::const_iterator itr( m_KeyToPathId.Find( cleaned_key));

      // if no id was found
      if( itr == m_KeyToPathId.End())
      {
        // return an empty model
        return util::GetUndefined< size_t>();
      }

      // find the size in the vector
      size_t &number_aas( m_ProteinFileSizes( itr->second));

      // if the value is known, it will be defined
      if( util::IsDefined( number_aas))
      {
        // so just return it
        return number_aas;
      }

      // try to read the index file, if it has not been read already
      if( !m_HaveLoadedIndexFile)
      {
        m_HaveLoadedIndexFile = true;
        std::string index_filename( GetIndexFileName());
        if( io::DirectoryEntry( index_filename).DoesExist())
        {
          io::IFStream index_read;
          if( io::File::TryOpenIFStream( index_read, index_filename))
          {
            util::ChopHeader( index_read);
            storage::Vector< storage::Vector< std::string> > path_to_counts
            (
              util::SplittedStringLineListFromIStream( index_read)
            );
            io::File::CloseClearFStream( index_read);
            for
            (
              storage::Vector< storage::Vector< std::string> >::const_iterator
                itr_index( path_to_counts.Begin()), itr_index_end( path_to_counts.End());
              itr_index != itr_index_end;
              ++itr_index
            )
            {
              if( itr_index->GetSize() != size_t( 2))
              {
                continue;
              }
              size_t key_index( m_ProteinFiles.Find( itr_index->FirstElement()));
              if( key_index < m_ProteinFiles.GetSize())
              {
                m_ProteinFileSizes( key_index) = util::ConvertStringToNumericalValue< size_t>( ( *itr_index)( 1));
              }
            }
          }
        }
      }

      // if the value is known, it will be defined
      if( util::IsDefined( number_aas))
      {
        // so just return it
        return number_aas;
      }

      // otherwise, the key will have to be read
      number_aas = Retrieve( KEY)->GetSize();

      if( m_CanWriteIndexFile)
      {
        std::string index_filename( GetIndexFileName());
        if( !index_filename.empty())
        {
          io::OFStream index_out;
          if( !io::DirectoryEntry( index_filename).DoesExist())
          {
            if( io::File::TryOpenOFStream( index_out, index_filename))
            {
              index_out << "# BCL Index file\n# Filename\tNumber Residues\n";
              index_out << "#This file is used to cache the number of residues in input files for BCL descriptor\n";
              index_out << "#generation and PDB file access.  It can be deleted at any time to force recalculation of sizes\n";
            }
            else
            {
              m_CanWriteIndexFile = false;
            }
          }
          else if( !io::File::TryOpenOFStream( index_out, index_filename, std::ios::app))
          {
            m_CanWriteIndexFile = false;
          }
          if( m_CanWriteIndexFile)
          {
            index_out << m_ProteinFiles( itr->second) << '\t' << number_aas << '\n';
            io::File::CloseClearFStream( index_out);
          }
        }
        else
        {
          m_CanWriteIndexFile = false;
        }
      }
      return number_aas;
    }

    //! @brief get the type of AA class that will be produced
    //! @return the type of AA class that will be produced
    const biol::AAClass &ProteinWithCacheStorageFile::GetAAClass() const
    {
      return m_Class;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get the filename for a given key
    //! @param KEY key identifier for specific object
    //! @return the filename for the given key (empty if the key is invalid for this storage)
    const std::string &ProteinWithCacheStorageFile::GetFilename( const std::string &KEY) const
    {
      static std::string s_empty;
      // look for the key
      storage::Map< std::string, size_t>::const_iterator itr( m_KeyToPathId.Find( GetKey( KEY, m_Suffix)));

      // get the complete filename, or empty if no filename was found
      return itr == m_KeyToPathId.End() ? s_empty : m_ProteinFiles( itr->second);
    }

    //! @brief get stored MoleculeComplete from key
    //! @param KEY key identifier for specific object
    //! @return Molecule of interest
    util::ShPtr< ProteinModelWithCache> ProteinWithCacheStorageFile::Retrieve( const std::string &KEY) const
    {
      // get the complete filename
      const std::string cleaned_key( GetKey( KEY, m_Suffix));

      // look for the model in the cache
      if( cleaned_key == m_LastModelKey)
      {
        return m_LastModel;
      }

      // look for the key
      storage::Map< std::string, size_t>::const_iterator itr( m_KeyToPathId.Find( cleaned_key));

      // handle unavailable keys
      if( itr == m_KeyToPathId.End())
      {
        BCL_MessageCrt( "Warning: Could not retrieve requested protein key: " + KEY);
        return util::ShPtr< ProteinModelWithCache>();
      }

      const std::string &filename( m_ProteinFiles( itr->second));

      BCL_MessageVrb( "Reading key: " + cleaned_key + " from " + filename);

      // create the model
      ProteinModel model( GenerateProteinModelFromFile( filename, m_RequireCoordinates, m_Class));

      // handle autoload and dssp parameters
      if( m_Autoload)
      {
        if( m_DoDSSP)
        {
          biol::DSSP().SetPDBSSPred( model);
        }
        if( m_Suffix.empty() || !util::EndsWith( m_Suffix, filename))
        {
          const std::string fname( io::File::RemoveLastExtension( filename));
          LoadBlastSspred( model, fname);
          if( cleaned_key.size() >= size_t( 5))
          {
            LoadBlastSspred( model, fname.substr( 0, fname.size() - 1));
          }
        }
        else
        {
          if( cleaned_key.size() >= size_t( 5))
          {
            LoadBlastSspred( model, filename.substr( 0, filename.size() - m_Suffix.size() - 1));
          }
          LoadBlastSspred( model, filename.substr( 0, filename.size() - m_Suffix.size()));
        }
      }
      m_LastModel = util::ShPtr< ProteinModelWithCache>( new ProteinModelWithCache( model, m_RequireCoordinates));
      if( m_Autoload && !m_BlastExtension.empty())
      {
        // check that the blast profile was read
        BCL_Assert
        (
          m_LastModel->IsEmpty() || m_LastModel->GetIterator()->GetBlastProfilePtr().IsDefined(),
          "Could not read blast profile for " + model.GetIdentification() + " @ " + filename
          + " with extension " + m_BlastExtension
        );
      }
      m_LastModelKey = cleaned_key;
      return m_LastModel;
    }

    //! @brief get molecule ensemble for given keys
    //! @param KEYS vector of keys
    //! @return list of MoleculeComplete objects corresponding to KEYS
    util::ShPtrVector< ProteinModelWithCache> ProteinWithCacheStorageFile::RetrieveEnsemble( const storage::Vector< std::string> &KEYS) const
    {
      // create ensemble with final proteins
      util::ShPtrVector< ProteinModelWithCache> list_proteins;
      list_proteins.AllocateMemory( KEYS.GetSize());
      for
      (
        storage::Vector< std::string>::const_iterator itr( KEYS.Begin()), itr_end( KEYS.End());
        itr != itr_end;
        ++itr
      )
      {
        list_proteins.PushBack( Retrieve( *itr));
      }

      // end
      return list_proteins;
    }

    //! @brief get ensemble of stored proteins objects for a range of keys, specified by plain numbers
    //! @param RANGE range of keys
    //! @return list of MoleculeComplete objects corresponding to keys in RANGE
    util::ShPtrVector< ProteinModelWithCache> ProteinWithCacheStorageFile::RetrieveEnsemble( const math::Range< size_t> &RANGE) const
    {
      // create a normal range
      if( !util::IsDefined( RANGE.GetMax()))
      {
        return this->RetrieveEnsemble( m_ProteinFiles);
      }

      math::Range< size_t> std_range( RANGE.StandardizeRange());

      // create ensemble with final proteins
      util::ShPtrVector< ProteinModelWithCache> list_proteins;
      list_proteins.AllocateMemory( std_range.GetWidth());
      for( size_t i( std_range.GetMin()), mx( std::min( m_ProteinFiles.GetSize(), std_range.GetMax())); i < mx; ++i)
      {
        list_proteins.PushBack( Retrieve( m_ProteinFiles( i)));
      }

      // end
      return list_proteins;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool ProteinWithCacheStorageFile::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      if( !m_PdbPathFile.empty())
      {
        m_Directory = io::File::MakeAbsolutePath( m_Directory);
        if( m_Directory.size() && m_Directory[ m_Directory.size() - 1] != '/')
        {
          m_Directory += '/';
        }
        if( !m_PdbIdListFile.empty())
        {
          BCL_MessageCrt( "key file parameter ignored because pdb paths were given");
        }
        io::IFStream input;
        io::File::MustOpenIFStream( input, m_PdbPathFile);
        storage::Vector< std::string> paths_of_interest( util::StringListFromIStream( input));
        io::File::CloseClearFStream( input);
        // remove any empty lines
        m_ProteinFiles.Reset();
        m_ProteinFiles.AllocateMemory( paths_of_interest.GetSize());
        for
        (
          storage::Vector< std::string>::const_iterator itr( paths_of_interest.Begin()), itr_end( paths_of_interest.End());
          itr != itr_end;
          ++itr
        )
        {
          const std::string trimmed_path( m_Directory + util::TrimString( *itr));
          const std::string trimmed_key( GetKey( trimmed_path, m_Suffix));
          if( !trimmed_key.empty())
          {
            if( !io::DirectoryEntry( trimmed_path).DoesExist())
            {
              ERR_STREAM << "Path from file does not exist: " << trimmed_path;
              return false;
            }
            if( m_KeyToPathId.Has( trimmed_key))
            {
              ERR_STREAM << "Unique key violation: " << trimmed_key << " encountered at "
                         << trimmed_path << " and " << m_ProteinFiles( m_KeyToPathId[ trimmed_key]);
              return false;
            }
            m_KeyToPathId[ trimmed_key] = m_ProteinFiles.GetSize();
            m_ProteinFiles.PushBack( trimmed_path);
          }
        }
      }
      // if the list file was set, filter the keys with the list
      else if( !m_PdbIdListFile.empty())
      {
        io::IFStream input;
        io::File::MustOpenIFStream( input, m_PdbIdListFile);
        storage::Vector< std::string> keys_of_interest( util::StringListFromIStream( input));
        io::File::CloseClearFStream( input);
        // remove any empty lines
        storage::Vector< std::string> valid_keys;
        valid_keys.AllocateMemory( keys_of_interest.GetSize());
        for
        (
          storage::Vector< std::string>::const_iterator itr( keys_of_interest.Begin()), itr_end( keys_of_interest.End());
          itr != itr_end;
          ++itr
        )
        {
          const std::string trimmed( GetKey( util::TrimString( *itr), m_Suffix));
          if( !trimmed.empty())
          {
            m_KeyToPathId[ trimmed] = valid_keys.GetSize();
            valid_keys.PushBack( trimmed);
          }
        }
        m_ProteinFiles.Resize( valid_keys.GetSize());

        storage::Set< std::string> desired_key_set( valid_keys.Begin(), valid_keys.End());

        // try finding all the desired keys
        if( !AddKeys( io::Directory( m_Directory), desired_key_set, ERR_STREAM))
        {
          return false;
        }

        // check that all desired keys were found
        if( !desired_key_set.IsEmpty())
        {
          ERR_STREAM << "Could not find the following keys: " << desired_key_set << '\n';
          return false;
        }
      }
      else if( !m_Directory.empty())
      {
        io::DirectoryEntry entry( m_Directory);
        if( entry.DoesExist() && entry.IsType( io::Directory::e_Dir))
        {
          // get all the keys
          if( !AddKeys( io::Directory( m_Directory), ERR_STREAM))
          {
            return false;
          }
        }
        else if( entry.IsType( io::Directory::e_File))
        {
          const std::string trimmed_key( GetKey( entry.GetName(), m_Suffix));
          m_KeyToPathId[ trimmed_key] = 0;
          m_ProteinFiles.PushBack( entry.GetFullName());
        }
      }

      // create the sizes vector, which will only be updated whenever AA count is requested for the given protein
      m_ProteinFileSizes.Resize( m_ProteinFiles.GetSize());
      m_ProteinFileSizes.SetAllElements( util::GetUndefined< size_t>());

      if( m_RequireCoordinates)
      {
        biol::Membrane::GetParameterTransitionThickness()->SetParameter( util::Format()( m_TransitionRegionThickness), ERR_STREAM);
        biol::Membrane::GetParameterGapThickness()->SetParameter( util::Format()( m_GapThickness), ERR_STREAM);
      }
      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProteinWithCacheStorageFile::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        m_RequireCoordinates
        ? "Retrieves PDB amino acids that have coordinates from a directory or subdirectory"
        : "Retrieves amino acid sequences from PDBs or FASTAs from a directory or subdirectory"
      );

      // set the parameter check accordingly
      parameters.AddOptionalInitializer
      (
        "",
        "Protein file or directory which contains proteins or folders containing proteins",
        io::Serialization::GetAgentInputFilename( &m_Directory)
      );

      parameters.AddOptionalInitializer
      (
        "suffix",
        "Protein file suffix; if not given, use .pdb if it is available, .fasta otherwise",
        io::Serialization::GetAgent( &m_Suffix)
      );

      parameters.AddOptionalInitializer
      (
        "key file",
        "File containing pdb ids, one per line, to look for in dir (by default, all files with the correct suffix are used)",
        io::Serialization::GetAgentInputFilename( &m_PdbIdListFile)
      );

      parameters.AddOptionalInitializer
      (
        "path file",
        "File containing absolute or relative paths to files that are to be loaded. "
        "If not given, the directory will be searched recursively for the keys in the key file, if it was given, or "
        "all files with the given suffix will be taken from the directory. "
        "This option is useful on files systems (e.g. NFS) where scanning directories is slow, also when there may be  "
        "multiple files in the directories that would yield the same key",
        io::Serialization::GetAgentInputFilename( &m_PdbPathFile)
      );

      if( m_RequireCoordinates)
      {
        parameters.AddInitializer
        (
          "aa class",
          "Level of detail at which to load amino acids",
          io::Serialization::GetAgent( &m_Class),
          m_RequireCoordinates ? biol::GetAAClasses().e_AABackBone.GetName() : biol::GetAAClasses().e_AA.GetName()
        );
        parameters.AddInitializer
        (
          "dssp",
          "Flag to enable re-determining Helix Strand and Coil regions via the DSSP algorithm",
          io::Serialization::GetAgent( &m_DoDSSP),
          "No"
        );
        parameters.AddInitializer
        (
          "membrane gap thickness",
          "Thickness (in angstrom) of gap between membrane and transition region",
          io::Serialization::GetAgentWithMin( &m_GapThickness, 0.0),
          "0"
        );
        parameters.AddInitializer
        (
          "membrane transition thickness",
          "Thickness of transition region; 5 A is thicness of typical phosholipid head group",
          io::Serialization::GetAgentWithMin( &m_TransitionRegionThickness, 0.0),
          "5"
        );
      }
      parameters.AddInitializer
      (
        "blast extension",
        "Extension for blast pssm files; if empty, choose .ascii6 if available, .ascii otherwise",
        io::Serialization::GetAgent( &m_BlastExtension),
        ""
      );
      parameters.AddInitializer
      (
        "recursive",
        "Recursively search all directories and subdirectory for the given keys. "
        "If this flag is not set, and no paths are given, only the given directory (or explicit file) will be searched for keys",
        io::Serialization::GetAgent( &m_RecursiveFind),
        "False"
      );
      parameters.AddDataMember( "key map", io::Serialization::GetAgent( &m_KeyToPathId));
      parameters.AddDataMember( "file list", io::Serialization::GetAgent( &m_ProteinFiles));
      return parameters;
    }

    //! @brief add protein keys and paths from the given directory
    //! @param DIRECTORY the directory of interest
    //! @param ERR_STREAM output stream for errors
    //! @return true if no duplicate keys were found
    bool ProteinWithCacheStorageFile::AddKeys( const io::Directory &DIRECTORY, std::ostream &ERR_STREAM)
    {
      static const std::string pdb_extension( "." + pdb::GetDefaultFileExtension());

      // determine the default extension
      const std::string def_extension( m_Suffix.empty() ? m_RequireCoordinates ? pdb_extension : ".fasta" : m_Suffix);
      {
        // get all files with the desired extension in the directory
        storage::List< io::DirectoryEntry> entries( DIRECTORY.ListFiles( "", def_extension, m_RecursiveFind));

        // add all the entries
        for
        (
          storage::List< io::DirectoryEntry>::const_iterator itr( entries.Begin()), itr_end( entries.End());
          itr != itr_end;
          ++itr
        )
        {
          // add the new file
          m_ProteinFiles.PushBack( itr->GetFullName());

          // get the key for this protein
          const std::string key( GetKey( itr->GetName(), m_Suffix));

          // ensure that the key is not already present could be added
          if( !m_KeyToPathId.Insert( std::make_pair( key, m_ProteinFiles.GetSize() - 1)).second)
          {
            ERR_STREAM << "Found " << key << " at " << m_ProteinFiles.LastElement()
                       << " and " << m_ProteinFiles( m_KeyToPathId[ key]);
            return false;
          }
        }
      }

      // if looking for sequences, also consider any pdb files
      if( m_Suffix.empty() && !m_RequireCoordinates)
      {
        // get all pdb files in the directory
        storage::List< io::DirectoryEntry> entries( DIRECTORY.ListFiles( "", pdb_extension, m_RecursiveFind));
        // add all the pdb entries, if a fasta has not already been found
        for
        (
          storage::List< io::DirectoryEntry>::const_iterator itr( entries.Begin()), itr_end( entries.End());
          itr != itr_end;
          ++itr
        )
        {
          // get the key for this pdb
          const std::string key( GetKey( itr->GetName(), m_Suffix));

          // ensure that the key is not already present and could be added
          if( m_KeyToPathId.Insert( std::make_pair( key, m_ProteinFiles.GetSize())).second)
          {
            // insert succeeded, so add the pdb key to the protein file
            m_ProteinFiles.PushBack( itr->GetFullName());
          }
        }
      }
      return true;
    }

    //! @brief add protein keys and paths from the given directory if they are in the given set, remove them from the set as each is added
    //! @param DIRECTORY the directory of interest
    //! @param DESIRED_KEYS set of keys of interest
    //! @param ERR_STREAM output stream for errors
    //! @return false if any duplicates are found
    bool ProteinWithCacheStorageFile::AddKeys
    (
      const io::Directory &DIRECTORY,
      storage::Set< std::string> &DESIRED_KEYS,
      std::ostream &ERR_STREAM
    )
    {
      static const std::string pdb_extension( "." + pdb::GetDefaultFileExtension());

      // determine the default extension
      const std::string def_extension( m_Suffix.empty() ? m_RequireCoordinates ? pdb_extension : ".fasta" : m_Suffix);
      {
        // get all files with the desired extension in the directory
        storage::List< io::DirectoryEntry> entries( DIRECTORY.ListFiles( "", def_extension, m_RecursiveFind));

        // add all the entries that are also in the desired keys
        for
        (
          storage::List< io::DirectoryEntry>::const_iterator itr( entries.Begin()), itr_end( entries.End());
          itr != itr_end;
          ++itr
        )
        {
          // get the key for this pdb
          const std::string key( GetKey( itr->GetName(), m_Suffix));

          // look for the key in the set
          storage::Set< std::string>::iterator itr_set( DESIRED_KEYS.Find( key));

          if( itr_set != DESIRED_KEYS.End())
          {
            m_ProteinFiles( m_KeyToPathId[ key]) = itr->GetFullName();
            DESIRED_KEYS.RemoveElement( itr_set);
          }
          else if( m_KeyToPathId.Has( key))
          {
            ERR_STREAM << "Found " << key << " at " << itr->GetFullName()
                       << " and " << m_ProteinFiles( m_KeyToPathId[ key]);
            return false;
          }
        }
      }

      // if looking for sequences, also consider any pdb files
      if( m_Suffix.empty() && !m_RequireCoordinates && !DESIRED_KEYS.IsEmpty())
      {
        // get all pdb files in the directory
        storage::List< io::DirectoryEntry> entries( DIRECTORY.ListFiles( "", pdb_extension, m_RecursiveFind));
        // add all the fasta entries, if a pdb has not already been found
        for
        (
          storage::List< io::DirectoryEntry>::const_iterator itr( entries.Begin()), itr_end( entries.End());
          itr != itr_end;
          ++itr
        )
        {
          // get the key for this fasta
          const std::string key( GetKey( itr->GetName(), m_Suffix));

          // look for the key in the set
          storage::Set< std::string>::iterator itr_set( DESIRED_KEYS.Find( key));

          if( itr_set != DESIRED_KEYS.End())
          {
            m_ProteinFiles( m_KeyToPathId[ key]) = itr->GetFullName();
            DESIRED_KEYS.RemoveElement( itr_set);
          }
        }
      }
      return true;
    }

    //! @brief given a fasta or pdb filename, generate a protein model
    //! @param FILENAME filename of interest
    //! @param REQUIRE_COORDINATES true if coordinates are required
    //! @param CLASS AA class to use
    ProteinModel ProteinWithCacheStorageFile::GenerateProteinModelFromFile
    (
      const std::string &FILENAME,
      const bool &REQUIRE_COORDINATES,
      const biol::AAClass &CLASS
    )
    {
      // get the complete filename
      const std::string cleaned_key( GetKey( FILENAME, "." + pdb::GetDefaultFileExtension()));

      BCL_MessageVrb( "Reading key: " + cleaned_key + " from " + FILENAME);

      // create a pdb factory
      pdb::Factory factory( CLASS);

      // Minimum SSE sizes (no filter)
      static storage::Map< biol::SSType, size_t> sse_sizes
      (
        storage::Map< biol::SSType, size_t>::Create
        (
          std::make_pair( biol::GetSSTypes().HELIX, 0),
          std::make_pair( biol::GetSSTypes().STRAND, 0),
          std::make_pair( biol::GetSSTypes().COIL, 0)
        )
      );

      // handle loading of pdbs
      if( REQUIRE_COORDINATES || io::File::GetLastExtension( FILENAME) == pdb::GetDefaultFileExtension())
      {
        // normal pdb file; read it in as usual, and add the identification from the key
        ProteinModel model( factory.ProteinModelFromPDBFilename( FILENAME, sse_sizes, true));
        util::ShPtr< ProteinModelData> sp_data( model.GetProteinModelData().HardCopy());
        sp_data->Insert
        (
          ProteinModelData::e_Identification,
          util::ShPtr< util::Wrapper< std::string> >( new util::Wrapper< std::string>( cleaned_key.substr( 0, 5)))
        );
        model.SetProteinModelData( sp_data);
        return model;
      }

      // fasta file, need to create a sequence and insert it into the protein model

      io::IFStream read;
      io::File::MustOpenIFStream( read, FILENAME);

      // protein model
      ProteinModel model;

      // chain id for sequence from last char of fasta filename without extension
      const char chain_id( cleaned_key.size() <= 4 || cleaned_key[ 4] == '_' || cleaned_key.size() > 5 ? ' ' : cleaned_key[ 4]);

      model.Insert( factory.ChainFromFastaStream( chain_id, read));
      io::File::CloseClearFStream( read);

      util::ShPtr< ProteinModelData> sp_data( model.GetProteinModelData().HardCopy());
      sp_data->Insert
      (
        ProteinModelData::e_PDBFile,
        util::ShPtr< util::Wrapper< std::string> >( new util::Wrapper< std::string>( FILENAME))
      );
      sp_data->Insert
      (
        ProteinModelData::e_Identification,
        util::ShPtr< util::Wrapper< std::string> >
        (
          new util::Wrapper< std::string>( cleaned_key.substr( 0, std::min( cleaned_key.size(), size_t( 6))))
        )
      );
      model.SetProteinModelData( sp_data);
      return model;
    }

    //! @brief load all extra data available for the model
    //! @param MODEL protein model of interest
    //! @param PREFIX filename of the model, without the pdb/fasta extension
    void ProteinWithCacheStorageFile::LoadBlastSspred( ProteinModel &MODEL, const std::string &PREFIX) const
    {
      biol::BlastProfileHandler::TryReadProfileForProteinModel( MODEL, PREFIX, m_BlastExtension);
      sspred::MethodHandler::ReadAllPredictionsForProteinModel( MODEL, PREFIX, "");
      sspred::PDB::SetEnvironmentTypes( MODEL, true);
      if( MODEL.GetChains().IsEmpty() || MODEL.GetChains().FirstElement()->GetSequence()->GetSize() == size_t( 0))
      {
        return;
      }
      const biol::AABase &first_aa( **( *MODEL.GetChains().FirstElement()->GetSequence()).Begin());
      if
      (
        first_aa.GetSSPrediction( sspred::GetMethods().e_StrideDSSP).IsDefined()
        && first_aa.GetSSPrediction( sspred::GetMethods().e_PALSSE).IsDefined()
        && !first_aa.GetSSPrediction( sspred::GetMethods().e_MAHSSMI).IsDefined()
        &&
        (
          first_aa.GetAAClass() == biol::GetAAClasses().e_AACaCb
          || first_aa.GetAAClass() == biol::GetAAClasses().e_AAComplete
          || first_aa.GetAAClass() == biol::GetAAClasses().e_AABackBone
        )
      )
      {
        sspred::Mahssmi::Calculate( MODEL, true);
      }
    }

    //! @brief helper function to get the key from a filename
    //! @param FILENAME filename of interest
    //! @param SUFFIX suffix that can be removed from the filename
    //! @return key
    std::string ProteinWithCacheStorageFile::GetKey( const std::string &FILENAME, const std::string &SUFFIX)
    {
      // create the basic key by removing the final extension, the full path
      std::string key( FILENAME);
      if( key.find( PATH_SEPARATOR) < std::string::npos)
      {
        key = io::File::RemovePath( key);
      }
      if( util::EndsWith( key, ".gz") || util::EndsWith( key, ".bz2"))
      {
        key = io::File::RemoveLastExtension( key);
      }
      if( !SUFFIX.empty() && util::EndsWith( key, SUFFIX))
      {
        key.erase( key.size() - SUFFIX.size());
      }
      if( util::EndsWith( key, ".fasta") || util::EndsWith( key, ".pdb"))
      {
        key = io::File::RemoveLastExtension( key);
      }
      // remove any trailing '_'
      size_t last_char( key.size() - 1);
      while( last_char > size_t( 3) && key[ last_char] == '_')
      {
        --last_char;
      }
      if( ++last_char != key.size())
      {
        key.erase( last_char);
      }
      if( key.size() < size_t( 5)) // no chain id
      {
        return util::ToLower( key);
      }

      // if a chain id is present, preserve its case for size 5 keys (presumably pdb id)
      return key.size() == size_t( 5) ? util::ToLower( key.substr( 0, 4)) + key.substr( 4) : key;
    }

    //! @brief helper function to determine the total number of residues in a map returned from pdb::Head
    //! @param MAP map, usually obtained from pdb::Head::GetSEQRESProteinChains or pdb::Head::GetMissingResidues
    //! @return the total number of residues in the map
    size_t ProteinWithCacheStorageFile::GetMapSize( const storage::Map< char, storage::List< pdb::ResidueSimple> > &MAP)
    {
      size_t total_size( 0);
      for
      (
        storage::Map< char, storage::List< pdb::ResidueSimple> >::const_iterator itr( MAP.Begin()), itr_end( MAP.End());
        itr != itr_end;
        ++itr
      )
      {
        total_size += itr->second.GetSize();
      }
      return total_size;
    }

    //! @brief helper function to determine the # of residues in a given PDB file from the atom lines
    //! @param ISTREAM input file
    size_t ProteinWithCacheStorageFile::CountResiduesInAtomLines( std::istream &INPUT)
    {
      // keep track of the total # of residues found
      size_t number_aas( 0);

      // create a vector (indexed by chain id) of valid seq ids for that chain
      storage::Vector< storage::Set< int> > valid_pdb_ids( 128);

      // skip everything up until the atom lines
      std::string line;
      while( INPUT.good() && !util::StartsWith( line, pdb::GetLineTypes().ATOM.GetName()))
      {
        std::getline( INPUT, line);
      }

      // get the end column of the occupancy column for the pdb
      const size_t chain_pos( pdb::GetEntryTypes().ATOMChainID->GetStart());
      const size_t seq_id_start( pdb::GetEntryTypes().ATOMResidueSequenceID->GetStart());
      const size_t seq_id_size( pdb::GetEntryTypes().ATOMResidueSequenceID->GetLength());
      const size_t het_three_letter_code_start( pdb::GetEntryTypes().HETATMResidueName->GetStart());
      const size_t seq_id_end( seq_id_start + seq_id_size);

      // temporary string for sequence id
      std::string seq_id_str( seq_id_size, ' ');

      // temporary error stream
      std::ostringstream err_stream;

      // scan each atom line
      while
      (
        util::StartsWith( line, pdb::GetLineTypes().ATOM.GetName())
        ||
        (
          // also allow HETATM entries that correspond to valid amino acid types, including UND
          // sometimes unnatural AAs like MSE are labeled with HETATM label rather than ATOM
          util::StartsWith( line, pdb::GetLineTypes().HETATM.GetName())
          &&
          (
            biol::GetAATypes().AATypeFromThreeLetterCode( line.substr( het_three_letter_code_start, 3)).IsDefined()
            || line.substr( het_three_letter_code_start, 3) == "UND"
          )
        )
      )
      {
        // if the line has all the necessary, members, extract the chain id and sequence id
        if( line.size() > seq_id_end)
        {
          // get the chain id
          const char chain_id( line[ chain_pos]);
          std::copy( line.begin() + seq_id_start, line.begin() + seq_id_end, seq_id_str.begin());
          int seq_id;
          if( util::TryConvertFromString( seq_id, seq_id_str, err_stream))
          {
            // sequence id was converted successfully
            // increment the # of AAs seen so far
            if( valid_pdb_ids( int( chain_id)).Insert( seq_id).second)
            {
              ++number_aas;
            }
          }
        }

        std::getline( INPUT, line);

        // skip empty lines and ANISOU lines
        while( INPUT.good() && ( line.empty() || util::StartsWith( line, pdb::GetLineTypes().ANISOU.GetName())))
        {
          std::getline( INPUT, line);
        }

        // catch the end of the file; in case this is a bad PDB file that does not have a proper END or TER or MASTER line
        if( !INPUT.good())
        {
          break;
        }
      }
      return number_aas;
    }

    //! @brief get the name of the index file for this protein cache storage file
    std::string ProteinWithCacheStorageFile::GetIndexFileName() const
    {
      if( io::DirectoryEntry( m_Directory).IsType( io::Directory::e_File))
      {
        return std::string();
      }
      std::string index_file_name( m_Directory + "bcl");

      // handle PDBs
      if( io::File::GetLastExtension( m_Suffix) == pdb::GetDefaultFileExtension())
      {
        index_file_name += ".seqres.counts";
      }
      else
      {
        // fastas are fast to load, so just read the file into an AASequence
        index_file_name += ".fasta.aa.counts";
      }
      return index_file_name;
    }

  } // namespace assemble
} // namespace bcl
