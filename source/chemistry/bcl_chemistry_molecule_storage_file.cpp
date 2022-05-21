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
#include "chemistry/bcl_chemistry_molecule_storage_file.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "chemistry/bcl_chemistry_molecule_ensemble.h"
#include "command/bcl_command_parameter_check_extension.h"
#include "command/bcl_command_parameter_check_or.h"
#include "crypt/bcl_crypt_sha1.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "sdf/bcl_sdf_mdl_entry_types.h"

// external includes - sorted alphabetically
#include <fstream>

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    // storage instance of this class
    const util::SiPtr< const util::ObjectInterface> MoleculeStorageFile::s_StoreInstance
    (
      util::Enumerated< io::StoreInterface< ConformationInterface> >::AddInstance
      (
        new MoleculeStorageFile( true)
      )
    );

    // retrieval instance of this class
    const util::SiPtr< const util::ObjectInterface> MoleculeStorageFile::s_RetrieveInstance
    (
      util::Enumerated< io::RetrieveInterface< MoleculeComplete, MoleculeEnsemble> >::AddInstance
      (
        new MoleculeStorageFile( false)
      )
    );

    //! Map of known hash caches, indexed by filename
    storage::Map
    <
      std::string,
      storage::Pair
      <
        sched::Mutex,
        storage::Vector< storage::Map< std::string, size_t> >
      >
    > MoleculeStorageFile::s_HashCache
     =
      storage::Map
      <
        std::string,
        storage::Pair
        <
          sched::Mutex,
          storage::Vector< storage::Map< std::string, size_t> >
        >
      >();

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from bool; which indicates whether the enumerated instance is intended for storage or retrieval
    MoleculeStorageFile::MoleculeStorageFile( const bool &STORE) :
      m_Store( STORE),
      m_EnsembleSize( 0),
      m_ConstitutionsCount( 0),
      m_FileSize( 0)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MoleculeStorageFile
    MoleculeStorageFile *MoleculeStorageFile::Clone() const
    {
      return new MoleculeStorageFile( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MoleculeStorageFile::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &MoleculeStorageFile::GetAlias() const
    {
      static const std::string s_name( "File");
      return s_name;
    }

    //! @brief number of molecules in storage
    //! @return number of molecules in this storage
    size_t MoleculeStorageFile::GetSize() const
    {
      long file_size( io::File::Size( m_Filename).Second());
      // if the file has been truncated, reset the file size
      if( file_size < m_FileSize)
      {
        m_FileSize = 0;
        m_MoleculeSizes.Reset();
      }
      if( file_size > m_FileSize)
      {
        // read the data set from the file
        io::IFStream input;
        io::File::MustOpenIFStream( input, m_Filename);

        // if the file did not have a terminal new line, call getline once to remove the terminal new line
        std::string temp_string;

        // fast forward to the old end of the file
        input.seekg( std::ifstream::off_type( m_FileSize), std::ios::beg);

        if( !input)
        {
          // seekg does not work on this type; reread the whole file
          io::File::CloseClearFStream( input);
          io::File::MustOpenIFStream( input, m_Filename);
          m_MoleculeSizes.Reset();
        }

        m_FileSize = file_size;
        // count the # of molecules in the file
        long lines_till_header_line( 4);
        while( input.good() && !input.eof())
        {
          std::getline( input, temp_string);
          if( sdf::MdlHandler::IsTerminalLine( temp_string))
          {
            lines_till_header_line = 4;
          }
          else if( !--lines_till_header_line)
          {
            m_MoleculeSizes.PushBack( sdf::GetMdlEntryTypes().Header_NumberAtomLines->GetUnsignedInt( temp_string));
          }
        }
      }

      return m_MoleculeSizes.GetSize();
    }

    //! @brief get all keys for given source
    //! @return all keys of given prefix
    storage::Vector< std::string> MoleculeStorageFile::GetAllKeys() const
    {
      // make a vector of strings big enough to hold the keys
      storage::Vector< std::string> keys( GetSize());

      for( size_t key( 0), max_key( keys.GetSize()); key < max_key; ++key)
      {
        keys( key) = util::Format()( key);
      }

      // return all created keys
      return keys;
    }

    //! @brief get the # of atoms for the molecule of the specified key (units defined by t_Type)
    //! @param KEY the identifier for the specific object
    //! @return the number of t_Type-defined units possessed by the key
    size_t MoleculeStorageFile::GetKeySize( const std::string &KEY) const
    {
      // call get size to force refreshing of the key sizes, if necessary
      GetSize();

      // index of molecule in sdf file
      const size_t molecule_index( util::ConvertStringToNumericalValue< size_t>( KEY));

      return
        molecule_index < m_MoleculeSizes.GetSize() ? m_MoleculeSizes( molecule_index) : util::GetUndefined< size_t>();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get stored MoleculeComplete from key
    //! @param KEY key identifier for specific object
    //! @return Molecule of interest
    MoleculeComplete MoleculeStorageFile::Retrieve( const std::string &KEY) const
    {
      // IFstream for reading in source sdf file
      io::IFStream read;
      // open sdf file for reading
      io::File::MustOpenIFStream( read, m_Filename);

      // index of molecule in sdf file
      const size_t molecule_index( util::ConvertStringToNumericalValue< size_t>( KEY));

      // create an ensemble containing only that molecule, add hydrogens as necessary
      MoleculeEnsemble chosen_molecule( read, math::Range< size_t>( molecule_index, molecule_index));

      // close sdf file stream
      io::File::CloseClearFStream( read);

      // index is out of range
      BCL_Assert
      (
        chosen_molecule.GetSize() == size_t( 1),
        "Could not retrieve molecule with key " + KEY + " from " + m_Filename
      );

      // return molecule
      return *chosen_molecule.Begin();
    }

    //! @brief get molecule ensemble for given keys
    //! @param KEYS vector of keys
    //! @return list of MoleculeComplete objects corresponding to KEYS
    MoleculeEnsemble MoleculeStorageFile::RetrieveEnsemble( const storage::Vector< std::string> &KEYS) const
    {
      // create ensemble with final molecules
      MoleculeEnsemble list_molecules;

      // if no keys were desired, just return the empty list
      if( KEYS.GetSize() == 0)
      {
        return list_molecules;
      }

      // IFstream for reading in source sdf file
      io::IFStream read;
      // open sdf file for reading
      io::File::MustOpenIFStream( read, m_Filename);

      // convert the strings to a vector of size_t's, which are the indices of the molecules in the sdf
      linal::Vector< size_t> molecule_keys( KEYS.GetSize());
      size_t *itr_sizet_keys( molecule_keys.Begin());

      // keep track of the min and max index seen so far
      // store the first element in the key
      *itr_sizet_keys = util::ConvertStringToNumericalValue< size_t>( KEYS.FirstElement());
      size_t min_key( *itr_sizet_keys);
      size_t max_key( *itr_sizet_keys);
      ++itr_sizet_keys;
      for
      (
        storage::Vector< std::string>::const_iterator itr_keys( KEYS.Begin() + 1), itr_keys_end( KEYS.End());
        itr_keys != itr_keys_end;
        ++itr_keys, ++itr_sizet_keys
      )
      {
        // store the index of molecule in sdf file
        *itr_sizet_keys = util::ConvertStringToNumericalValue< size_t>( *itr_keys);

        // update the min/max
        if( *itr_sizet_keys > max_key)
        {
          max_key = *itr_sizet_keys;
        }
        else if( *itr_sizet_keys < min_key)
        {
          min_key = *itr_sizet_keys;
        }
      }

      // create ensemble containing all keys in the range of keys found
      MoleculeEnsemble ensemble( read, math::Range< size_t>( min_key, max_key));

      // close sdf file stream
      io::File::CloseClearFStream( read);

      // iterator on end of molecule ensemble
      util::SiPtrVector< MoleculeComplete> molecules;
      molecules.AllocateMemory( ensemble.GetSize());
      for
      (
        storage::List< MoleculeComplete>::iterator itr( ensemble.Begin()), itr_end( ensemble.End());
        itr != itr_end;
        ++itr
      )
      {
        molecules.PushBack( *itr);
      }

      // determine molecule according to KEYs
      for
      (
        const size_t *itr_keys( molecule_keys.Begin()), *itr_keys_end( molecule_keys.End());
        itr_keys != itr_keys_end;
        ++itr_keys
      )
      {
        // fill list with relevant mdl handler
        list_molecules.PushBack( *molecules( *itr_keys - min_key));
      }

      // end
      return list_molecules;
    }

    //! @brief get ensemble of stored molecules objects for a range of keys, specified by plain numbers
    //! @param RANGE range of keys
    //! @return list of MoleculeComplete objects corresponding to keys in RANGE
    MoleculeEnsemble MoleculeStorageFile::RetrieveEnsemble( const math::Range< size_t> &RANGE) const
    {
      // IFstream for reading in source sdf file
      io::IFStream read;
      // open sdf file for reading
      io::File::MustOpenIFStream( read, m_Filename);
      // create ensemble containing all keys in the range of keys found
      MoleculeEnsemble ensemble( read, RANGE);
      // close sdf file stream
      io::File::CloseClearFStream( read);
      // return the molecules that were loaded in
      return ensemble;
    }

    //! @brief store molecule
    //! @param MOLECULE molecule to store
    //! @return key associated with molecule
    std::string MoleculeStorageFile::Store( const ConformationInterface &MOLECULE)
    {
      // call get size to populate the sizes of the size array if it was empty
      if( m_MoleculeSizes.IsEmpty())
      {
        GetSize();
      }

      // find whether the molecule is already in the sdf file, and if so, where it is at
      const std::pair< size_t, bool> pos_and_existence( Find( MOLECULE));
      const std::string key( util::Format()( pos_and_existence.first));
      // if the molecule is already in the datafile, just return its key
      if( pos_and_existence.second)
      {
        return key;
      }

      if
      (
        ( *m_MoleculeHashIndexMap)( size_t( e_Constitution)).Insert
        (
            std::make_pair
          (
            crypt::Sha1().Hash( sdf::MdlHandler::CreateConstitutionalHashString( MOLECULE.GetAtomInfo(), MOLECULE.GetBondInfo())),
            m_ConstitutionsCount
          )
        ).second
      )
      {
        ++m_ConstitutionsCount;
      }

      ( *m_MoleculeHashIndexMap)( size_t( e_Conformation)).Insert
      (
        std::make_pair
        (
          crypt::Sha1().Hash( sdf::MdlHandler::CreateConformationalHashString( MOLECULE.GetAtomInfo(), MOLECULE.GetBondInfo())),
          m_EnsembleSize
        )
      );
      const std::string key_new( util::Format()( m_EnsembleSize));
      ++m_EnsembleSize;

      // OFstream for writing into sdf file
      io::OFStream write;
      // open sdf file for reading
      io::File::MustOpenOFStream( write, m_Filename, std::ios::app);

      // if there are any valences, construct a new molecule
      if( MOLECULE.GetNumberValences())
      {
        MoleculeComplete new_mol( MOLECULE);
        m_MoleculeSizes.PushBack( new_mol.GetNumberAtoms());
        new_mol.WriteMDL( write);
      }
      else
      {
        // write molecule directly
        m_MoleculeSizes.PushBack( MOLECULE.GetNumberAtoms());
        MOLECULE.WriteMDL( write);
      }
      // close sdf file stream
      io::File::CloseClearFStream( write);

      // increment key
      return key_new;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool MoleculeStorageFile::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      io::DirectoryEntry file( m_Filename);
      m_Filename = file.GetFullName();
      // if storing molecules, create the file if it does not already exist
      if( m_Store)
      {
        if( !file.DoesExist())
        {
          io::OFStream output;
          io::File::MustOpenOFStream( output, m_Filename);
          io::File::CloseClearFStream( output);
        }
      }
      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeStorageFile::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        m_Store ? "Stores molecules in a file" : "Retrieves molecules from a file"
      );

      // set the parameter check accordingly
      if( m_Store)
      {
        // for the storage, the file may or may not already exist, but the extension should be valid
        parameters.AddInitializer
        (
          "",
          "Filename",
          io::Serialization::GetAgentWithCheck
          (
            &m_Filename,
            command::ParameterCheckOr
            (
              command::ParameterCheckExtension( ".sdf"),
              command::ParameterCheckExtension( ".sdf.bz2"),
              command::ParameterCheckExtension( ".sdf.gz")
            )
          )
        );
      }
      else
      {
        parameters.AddInitializer
        (
          "",
          "Filename",
          io::Serialization::GetAgentInputFilename( &m_Filename)
        );
      }

      return parameters;
    }
  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MoleculeStorageFile::Read( std::istream &ISTREAM)
    {
      io::StoreInterface< ConformationInterface>::Read( ISTREAM);
      return ISTREAM;
    }

    //! @brief read from std::ostream
    //! @param OSTREAM input stream
    //! @param INDENT indentation
    //! @return ostream which was read from
    std::ostream &MoleculeStorageFile::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::StoreInterface< ConformationInterface>::Write( OSTREAM, INDENT);
      return OSTREAM;
    }

    //! @brief Alignment as string
    //! @param ALIGNMENT the alignment
    //! @return the Alignment as string
    const std::string &MoleculeStorageFile::GetFindTypeString( const FindType &TYPE)
    {
      static std::string s_names[ s_NumberOfFindTypes + 1] =
      {
        "Constitution",
        "Conformation",
        ""
      };
      return s_names[ TYPE];
    }

    //! @brief find a molecule in the file; must have exactly matching conformation, but misc properties are ignored
    //! @return pair; bool indicates whether the molecule is already in the file, size_t indicates position in the file,
    //! or, if the molecule is not in the file, the size of the file
    std::pair< size_t, bool> MoleculeStorageFile::Find( const ConformationInterface &MOLECULE, const FindType &TYPE) const
    {
      if( !m_MoleculeHashIndexMap.IsDefined())
      {
        CacheHashes( TYPE);
      }
      storage::Map< std::string, size_t>::const_iterator itr_map
      (
        ( *m_MoleculeHashIndexMap)( size_t( TYPE)).Find
        (
          crypt::Sha1().Hash
          (
            TYPE == e_Conformation
            ? sdf::MdlHandler::CreateConformationalHashString( MOLECULE.GetAtomInfo(), MOLECULE.GetBondInfo())
            : sdf::MdlHandler::CreateConstitutionalHashString( MOLECULE.GetAtomInfo(), MOLECULE.GetBondInfo())
          )
        )
      );
      if( itr_map != ( *m_MoleculeHashIndexMap)( size_t( TYPE)).End())
      {
        // return the index, which is the number of molecules in the file, and false, to indicate the molecule
        // is not already in the file
        return std::make_pair( itr_map->second, true);
      }
      return std::make_pair( TYPE == e_Conformation ? m_EnsembleSize : m_ConstitutionsCount, false);
    }

    void MoleculeStorageFile::CacheHashes( const FindType &TYPE) const
    {
      if( !m_HashCacheMutex.IsDefined())
      {
        static sched::Mutex s_MapMutex;
        s_MapMutex.Lock();
        auto &res( s_HashCache[ m_Filename]);
        m_HashCacheMutex = util::ToSiPtrNonConst( res.First());
        m_MoleculeHashIndexMap = util::ToSiPtrNonConst( res.Second());
        s_MapMutex.Unlock();
      }
      m_HashCacheMutex->Lock();
      if( m_MoleculeHashIndexMap->IsEmpty())
      {
        m_MoleculeHashIndexMap->Resize( 2);
      }
      else if( !( *m_MoleculeHashIndexMap)( size_t( TYPE)).IsEmpty())
      {
        // already constructed hashes
        m_EnsembleSize = ( *m_MoleculeHashIndexMap)( size_t( e_Conformation)).GetSize();
        m_ConstitutionsCount = ( *m_MoleculeHashIndexMap)( size_t( e_Constitution)).GetSize();
        m_HashCacheMutex->Unlock();
        return;
      }
      // open sdf file for reading
      FragmentFeed feed( storage::Vector< std::string>::Create( m_Filename), sdf::e_Remove);
      if( TYPE == e_Constitution)
      {
        m_ConstitutionsCount = size_t( 0);
      }
      ( *m_MoleculeHashIndexMap)( size_t( TYPE)).Reset();

      // read in from the stream until we reach the end of the file or the last index
      for( ; feed.NotAtEnd(); ++feed)
      {
        storage::Vector< sdf::AtomInfo> atom_info( feed->GetAtomInfo());
        storage::Vector< sdf::BondInfo> bond_info( feed->GetBondInfo());

        if( TYPE == e_Conformation)
        {
          const std::string conformation_sha1_hash
          (
            crypt::Sha1().Hash( sdf::MdlHandler::CreateConformationalHashString( atom_info, bond_info))
          );
          ( *m_MoleculeHashIndexMap)( size_t( TYPE)).Insert( std::make_pair( conformation_sha1_hash, feed.GetPosition()));
        }

        if( TYPE == e_Constitution)
        {
          const std::string constitutional_sha1_hash
          (
            crypt::Sha1().Hash( sdf::MdlHandler::CreateConstitutionalHashString( atom_info, bond_info))
          );

          if( ( *m_MoleculeHashIndexMap)( size_t( TYPE)).Insert( std::make_pair( constitutional_sha1_hash, m_ConstitutionsCount)).second)
          {
            ++m_ConstitutionsCount;
          }
        }
      }
      m_EnsembleSize = feed.GetPosition();
      m_HashCacheMutex->Unlock();
    }
  } // namespace chemistry
} // namespace bcl
