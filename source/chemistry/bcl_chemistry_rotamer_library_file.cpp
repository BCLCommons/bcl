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
#include "chemistry/bcl_chemistry_rotamer_library_file.h"

// includes from bcl - sorted alphabetically
#include "bcl_version.h"
#include "chemistry/bcl_chemistry_fragment_constitution_shared.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "command/bcl_command_command_state.h"
#include "command/bcl_command_parameter_check_or.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_loggers.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically

//path for models
#if defined (__MINGW32__)
  #define BCL_ROTAMER_LIBRARY_PATH "../../../rotamer_library/"
#elif defined (__GNUC__)
  #define BCL_ROTAMER_LIBRARY_PATH "/dors/meilerlab/apps/bcl/rotamer_library/rev_5374/"
#elif defined (_MSC_VER)
  #define BCL_ROTAMER_LIBRARY_PATH "../../../rotamer_library/"
#endif

#define BCL_PROFILE_RotamerLibraryFile

namespace bcl
{
  namespace chemistry
  {

    //! @brief create the command line file search path object
    command::ParameterCheckFileInSearchPath RotamerLibraryFile::GetRotamerFinder()
    {
      return command::ParameterCheckFileInSearchPath
             (
               "rotamer_library",
               GetVersion().IsLicense()
               ? GetVersion().GetInstallPrefix() + "/rotamer_library/"
               : BCL_ROTAMER_LIBRARY_PATH,
               io::Directory::e_Dir
             );
    }

  //////////
  // data //
  //////////

    // storage instance of this class
    const util::SiPtr< const util::ObjectInterface> RotamerLibraryFile::s_StoreInstance
    (
      util::Enumerated< RotamerLibraryInterface>::AddInstance( new RotamerLibraryFile())
    );

    // storage instance of this class
    const util::SiPtr< const util::ObjectInterface> RotamerLibraryFile::s_CODInstance
    (
      util::Enumerated< RotamerLibraryInterface>::AddInstance( new RotamerLibraryFile( "cod"))
    );

    util::Stopwatch RotamerLibraryFile::s_ReadingTimer
    (
      "RotamerLibraryFile reading files",
      util::Message::e_Standard,
      true,
      false
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor, try to find the rotamer library
    RotamerLibraryFile::RotamerLibraryFile()
      : m_FileNumber( size_t( 100)),
        m_Alias( "")
    {
    }

    //! @brief constructor, try to find the rotamer library
    RotamerLibraryFile::RotamerLibraryFile( const std::string &ALIAS)
    : m_FileNumber( size_t( 100)),
      m_Alias( ALIAS)
    {
    }

    //! @brief Clone function
    //! @return pointer to new RotamerLibraryFile
    RotamerLibraryFile *RotamerLibraryFile::Clone() const
    {
      return new RotamerLibraryFile( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RotamerLibraryFile::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &RotamerLibraryFile::GetAlias() const
    {
      static const std::string s_name( "File");

      return m_Alias.empty() ? s_name : m_Alias;
    }

    //! @brief get the filename
    //! @return the filename
    std::string RotamerLibraryFile::GetFileExtension() const
    {
      static const std::string s_gz( "GZ"), s_uncompressed( "Uncompressed"), s_bz2( "BZ2");
      static const std::string s_gz_lower( ".gz"), s_empty( ""), s_bz2_lower( ".bz2");
      static const std::string s_preferred_compression
      (
        // prefer compressed files in the order: GZ > BZ2 > None. GZ is faster at decompression and the size difference
        // doesn't warrant BZ2s much slower decompression rate
        io::GetStreamBufferClasses().HaveEnumWithName( s_gz)
        ? s_gz_lower
        : ( io::GetStreamBufferClasses().HaveEnumWithName( s_bz2) ? s_bz2_lower : s_empty)
      );
      const std::string ret
      (
        m_Compression.IsDefined() && m_Compression->IsDefined()
          ? ( *m_Compression)->AddExtension( "")
          : s_preferred_compression
      );
      return ret;
    }

    //! @brief detect whether the rotamer library is defined
    bool RotamerLibraryFile::IsDefined() const
    {
      if( !m_Lines.IsEmpty())
      {
        return true;
      }
      if( m_Filename.empty())
      {
        m_Filename = GetRotamerFinder().FindFile( "") + GetPrefix();
      }
      return io::DirectoryEntry( m_Filename + ".constitutions.txt" + GetFileExtension()).DoesExist();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get stored MoleculeComplete from key
    //! @param KEY key identifier for specific object
    //! @return Molecule of interest
    graph::ConstGraph< size_t, size_t> RotamerLibraryFile::RetrieveRotamerSubstructureTree() const
    {
      s_ReadingTimer.Start();
      if( m_Lines.IsEmpty())
      {
        LoadFiles();
      }
      // IFstream for reading in source sdf file
      io::IFStream read;
      // open sdf file for reading
      io::File::MustOpenIFStream( read, m_Filename + ".substructure.txt" + GetFileExtension());

      graph::ConstGraph< size_t, size_t> constitution_connectivity;
      constitution_connectivity.ReadBasicConnectivity( read);

      // close sdf file stream
      io::File::CloseClearFStream( read);
      s_ReadingTimer.Stop();
      // return molecule
      return constitution_connectivity;
    }

    //! @brief get the directed-graph whose vertices are constitution ids and child nodes are substructures of parent node
    graph::ConstGraph< size_t, size_t> RotamerLibraryFile::RetrieveRotamerRequirementsGraph() const
    {
      s_ReadingTimer.Start();
      if( m_Lines.IsEmpty())
      {
        LoadFiles();
      }
      // IFstream for reading in source sdf file
      io::IFStream read;
      // open sdf file for reading
      graph::ConstGraph< size_t, size_t> constitution_connectivity;
      if( io::File::TryOpenIFStream( read, m_Filename + ".requirements.txt" + GetFileExtension()))
      {
        constitution_connectivity.ReadBasicConnectivity( read);

        // close sdf file stream
        io::File::CloseClearFStream( read);
      }
      s_ReadingTimer.Stop();
      // return molecule
      return constitution_connectivity;
    }

    //! @brief get molecule ensemble for given keys
    //! @param KEYS vector of keys
    //! @return list of MoleculeComplete objects corresponding to KEYS
    storage::Vector< graph::ConstGraph< size_t, size_t> > RotamerLibraryFile::RetrieveConstitution( const storage::Vector< size_t> &IDS) const
    {
      s_ReadingTimer.Start();
      storage::Vector< graph::ConstGraph< size_t, size_t> > constitution_graphs;
      constitution_graphs.AllocateMemory( IDS.GetSize());

      for
      (
        storage::Vector< size_t>::const_iterator itr( IDS.Begin()), itr_end( IDS.End());
          itr != itr_end;
        ++itr
      )
      {
        size_t desired_line( *itr * 2 + size_t( 1));
        constitution_graphs.PushBack( CreateConstitutionGraphFromString( m_Lines( desired_line), m_Lines( desired_line + 1)));
      }
      s_ReadingTimer.Stop();
      return constitution_graphs;
    }

    //! @brief get constitutions that are root in the substructure tree i.e. they are not contained in any other fragments
    storage::Vector< graph::ConstGraph< size_t, size_t> > RotamerLibraryFile::GetRootConstitutions() const
    {
      s_ReadingTimer.Start();
      if( m_Lines.IsEmpty())
      {
        LoadFiles();
      }
      size_t number_root_nodes( util::ConvertStringToNumericalValue< size_t>( m_Lines( 0)( 0)));

      storage::Vector< graph::ConstGraph< size_t, size_t> > constitution_graphs;
      constitution_graphs.AllocateMemory( number_root_nodes);

      for( size_t root_node_id( 0); root_node_id < number_root_nodes; ++root_node_id)
      {
        const size_t desired_line( root_node_id * 2 + size_t( 1));
        constitution_graphs.PushBack( CreateConstitutionGraphFromString( m_Lines( desired_line), m_Lines( desired_line + 1)));
      }
      s_ReadingTimer.Stop();
      return constitution_graphs;
    }

    //! @brief get ensemble of stored molecules objects for a range of keys, specified by plain numbers
    //! @param RANGE range of keys
    //! @return list of MoleculeComplete objects corresponding to keys in RANGE
    storage::Vector< storage::Set< size_t> > RotamerLibraryFile::RetrieveConfigurationMapping() const
    {
      s_ReadingTimer.Start();
      if( m_Lines.IsEmpty())
      {
        LoadFiles();
      }
      // IFstream for reading in source sdf file
      io::IFStream read;
      // open sdf file for reading
      io::File::MustOpenIFStream( read, m_Filename + ".configuration_mapping.txt" + GetFileExtension());

      storage::Vector< storage::Vector< std::string> > lines( util::SplittedStringLineListFromIStream( read, " <->,"));
      const size_t number_constiutions( lines.GetSize());

      storage::Vector< storage::Set< size_t> > mapping;
      mapping.AllocateMemory( number_constiutions);

      // iterate through the vertices
      for( size_t index_a( 0); index_a < number_constiutions; index_a++)
      {
        storage::Set< size_t> associated_configurations;
        for
        (
          storage::Vector< std::string>::const_iterator
            itr( lines( index_a).Begin() + 1), itr_end( lines( index_a).End());
           itr != itr_end;
           ++itr
        )
        {
          associated_configurations.Insert( util::ConvertStringToNumericalValue< size_t>( *itr));
        }
        mapping.PushBack( associated_configurations);
      }

      // close sdf file stream
      io::File::CloseClearFStream( read);
      s_ReadingTimer.Stop();
      return mapping;
    }

    //! @brief store molecule
    //! @param MOLECULE molecule to store
    //! @return key associated with molecule
    void RotamerLibraryFile::CreateImpl
    (
      const size_t &NUMBER_OF_ROOT_NODES,
      const storage::Vector< storage::Set< size_t> > &ROTAMER_SUBSTRUCTURE_TREE,
      const util::ShPtrVector< FragmentConstitutionShared> &CONSTITUTIONS,
      const storage::Vector< storage::Set< size_t> > &CONSTITUTION_TO_CONFORMATION_MAPPING,
      const FragmentEnsemble &CONFORMATIONS,
      const storage::Vector< size_t> &PARENTS
    ) const
    {
      size_t number_unique_constitutions( CONSTITUTIONS.GetSize());

      storage::Vector< storage::Vector< size_t> > child_nodes_from_parents( number_unique_constitutions);
      if( !PARENTS.IsEmpty())
      {
        for( size_t i( 0), sz( PARENTS.GetSize()); i < sz; ++i)
        {
          if( util::IsDefined( PARENTS( i)))
          {
            child_nodes_from_parents( PARENTS( i)).PushBack( i);
          }
        }
      }

      // OFstream for writing into sdf file
      io::OFStream write;

      // open sdf file for reading
      io::File::MustOpenOFStream( write, m_Filename + ".substructure.txt" + GetFileExtension());
      for( size_t constitution_index( 0); constitution_index < number_unique_constitutions; ++constitution_index)
      {
        write << constitution_index << " <-> ";
        if( PARENTS.IsEmpty())
        {
          const storage::Set< size_t> &child_nodes( ROTAMER_SUBSTRUCTURE_TREE( constitution_index));
          for
          (
            storage::Set< size_t>::const_iterator itr_set( child_nodes.Begin()), itr_set_end( child_nodes.End());
              itr_set != itr_set_end;
            ++itr_set
          )
          {
            write << *itr_set << ", ";
          }
        }
        else
        {
          const storage::Vector< size_t> &child_nodes( child_nodes_from_parents( constitution_index));
          for
          (
            storage::Vector< size_t>::const_iterator itr_set( child_nodes.Begin()), itr_set_end( child_nodes.End());
              itr_set != itr_set_end;
            ++itr_set
          )
          {
            write << *itr_set << ", ";
          }
        }
        write << '\n';
      }

      // close sdf file stream
      io::File::CloseClearFStream( write);
      if( !PARENTS.IsEmpty())
      {
        io::File::MustOpenOFStream( write, m_Filename + ".requirements.txt" + GetFileExtension());

        for( size_t constitution_index( 0); constitution_index < number_unique_constitutions; ++constitution_index)
        {
          write << constitution_index << " <-> ";
          const storage::Set< size_t> &child_nodes( ROTAMER_SUBSTRUCTURE_TREE( constitution_index));
          for
          (
            storage::Set< size_t>::const_iterator itr_set( child_nodes.Begin()), itr_set_end( child_nodes.End());
              itr_set != itr_set_end;
            ++itr_set
          )
          {
            write << *itr_set << ", ";
          }
          write << '\n';
        }

        io::File::CloseClearFStream( write);
      }

      const ConfigurationalBondTypeData::Data bond_data
      (
        ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness
      );

      io::File::MustOpenOFStream( write, m_Filename + ".constitutions.txt" + GetFileExtension());

      write << util::Format()( NUMBER_OF_ROOT_NODES) << "\n";

      // Iterate through all molecules in the input ensemble
      for
      (
        util::ShPtrVector< FragmentConstitutionShared>::const_iterator
          itr_fragment( CONSTITUTIONS.Begin()), itr_fragment_end( CONSTITUTIONS.End());
        itr_fragment != itr_fragment_end;
        ++itr_fragment
      )
      {
        std::stringstream molecules_string;
        molecules_string << util::Format()( ( *itr_fragment)->GetNumberAtoms()) << " ";
        molecules_string << util::Format()( ( *itr_fragment)->GetNumberBonds()) << " ";
        molecules_string << "\n";

        const AtomVector< AtomComplete> atom_vector( ( *itr_fragment)->GetAtomInfo(), ( *itr_fragment)->GetBondInfo());
        storage::Vector< graph::UndirectedEdge< size_t> > adjacency_list( atom_vector.GetAdjacencyList( bond_data));
        for
        (
          storage::Vector< graph::UndirectedEdge< size_t> >::const_iterator
            itr_adj( adjacency_list.Begin()), itr_adj_end( adjacency_list.End());
          itr_adj != itr_adj_end;
          ++itr_adj
        )
        {
          size_t index_low( itr_adj->GetIndexLow());
          size_t index_high( itr_adj->GetIndexHigh());
          size_t edge_data( itr_adj->GetEdgeData());
          molecules_string << util::Format()( index_low) << " ";
          molecules_string << util::Format()( atom_vector( index_low).GetAtomType()->GetIndex()) << " ";
          molecules_string << util::Format()( index_high) << " ";
          molecules_string << util::Format()( atom_vector( index_high).GetAtomType()->GetIndex()) << " ";
          molecules_string << util::Format()( edge_data) << " ";
        }

        molecules_string << "\n";
        write << molecules_string.str();
      }
      io::File::CloseClearFStream( write);

      // Iterate through all molecules in the input ensemble
      io::File::MustOpenOFStream( write, m_Filename + ".configuration_mapping.txt" + GetFileExtension());
      size_t constitution( 0);
      for
      (
        storage::Vector< storage::Set< size_t> >::const_iterator
          itr( CONSTITUTION_TO_CONFORMATION_MAPPING.Begin()), itr_end( CONSTITUTION_TO_CONFORMATION_MAPPING.End());
        itr != itr_end;
        ++itr, ++constitution
      )
      {
        const storage::Set< size_t> &child_nodes( *itr);
        write << constitution << " <-> ";
        for
        (
          storage::Set< size_t>::const_iterator itr_set( child_nodes.Begin()), itr_set_end( child_nodes.End());
            itr_set != itr_set_end;
          ++itr_set
        )
        {
          write << *itr_set << ", ";
        }
        write << '\n';
      }
      io::File::CloseClearFStream( write);

      std::string dir_prefix( m_Filename + "_conformations");
      io::Directory new_directory( dir_prefix);
      new_directory.Make();

      // OFstream for writing into sdf file
      io::OFStream write_conformations;

      size_t molecule_index( 0);
      storage::List< FragmentComplete> molecule_chunk;
      size_t total_conformations( CONFORMATIONS.GetSize() - 1);
      std::string second_dir_prefix;
      // Iterate through all conformations, write them out in chunks of 100
      for
      (
        storage::List< FragmentComplete>::const_iterator
          itr_fragment( CONFORMATIONS.Begin()), itr_fragment_end( CONFORMATIONS.End());
        itr_fragment != itr_fragment_end;
        ++itr_fragment, ++molecule_index
      )
      {
        if( molecule_index % size_t( m_FileNumber * 100) == 0)
        {
          size_t dir_chunk( molecule_index / size_t( m_FileNumber * 100));
          second_dir_prefix = dir_prefix + "/" + util::Format()( dir_chunk);
          io::Directory new_directory( second_dir_prefix);
          new_directory.Make();
        }
        FragmentComplete tmp_mol( *itr_fragment);
        molecule_chunk.PushBack( tmp_mol);
        if( molecule_index % size_t( 100) == size_t( 99) || molecule_index == total_conformations)
        {
          FragmentEnsemble molecules( molecule_chunk);
          size_t prefix( ( molecule_index - size_t( 99)) / size_t( 100));

          if( molecule_index == total_conformations)
          {
            prefix += 1;
          }

          if( molecule_index < 100)
          {
            prefix = size_t( 0);
          }

          std::string file_prefix( second_dir_prefix + "/" + util::Format()( prefix) + ".sdf" + GetFileExtension());

          // open sdf file for reading
          io::File::MustOpenOFStream( write_conformations, file_prefix);
          molecules.WriteMDL( write_conformations);
          io::File::CloseClearFStream( write_conformations);
          molecule_chunk.Reset();
        }
      }
      // IFstream for reading in source sdf file containing constitutions
      io::IFStream read;
      // open sdf file for reading
      io::File::MustOpenIFStream( read, m_Filename + ".constitutions.txt" + GetFileExtension());

      // split the file into a format so that it can be iterated
      m_Lines = util::SplittedStringLineListFromIStream( read, " ");

      io::File::CloseClearFStream( read);
    }

    FragmentEnsemble RotamerLibraryFile::RetrieveAssociatedConfigurations( const storage::Set< size_t> &IDS) const
    {
      s_ReadingTimer.Start();
      FragmentEnsemble configurations;

      if( !io::DirectoryEntry( m_Filename + "_conformations/").DoesExist())
      {
        m_Filename = GetRotamerFinder().FindFile( "") + GetPrefix();
      }

      // OFstream for writing into sdf file
      FragmentEnsemble temp_ensemble;
      storage::Vector< FragmentComplete> conformation_vector;
      for
      (
        storage::Set< size_t>::const_iterator itr( IDS.Begin()), itr_end( IDS.End());
        itr != itr_end;
        //++itr
      )
      {
        size_t configuration_id( *itr);
        size_t prefix( ( *itr) / size_t( 100));
        storage::Set< size_t>::const_iterator itr_start( itr);
        for( ++itr; itr != itr_end && prefix == ( *itr) / size_t( 100); ++itr)
        {
        }
        const size_t max_index( *--itr);
        ++itr;
        io::IFStream conformation_file;
        size_t dir_prefix( prefix / size_t( m_FileNumber));

        std::string file_prefix( m_Filename + "_conformations/" + util::Format()( dir_prefix) + "/" + util::Format()( prefix) + ".sdf" + GetFileExtension());
        const size_t molecule_start( configuration_id % size_t( 100));
        const size_t molecule_end( max_index % size_t( 100));
        math::Range< size_t> range( molecule_start, molecule_end);
        io::File::MustOpenIFStream( conformation_file, file_prefix);
        FragmentEnsemble temp_ensemble( conformation_file, sdf::e_Maintain, range, sdf::e_None);
        conformation_vector = storage::Vector< FragmentComplete>( temp_ensemble.Begin(), temp_ensemble.End());
        io::File::CloseClearFStream( conformation_file);
        for( ; itr_start != itr; ++itr_start)
        {
          const size_t new_mol_index( ( *itr_start % size_t( 100)) - molecule_start);
          // Iterate through all molecules in the input ensemble
          configurations.PushBack( conformation_vector( new_mol_index));
        }
        // close sdf file stream
      }
      s_ReadingTimer.Stop();
      // return the molecules that were loaded in
      return configurations;
    }

    //! @brief get the bond angle map
    RotamerLibraryFile::t_BondAngleMap RotamerLibraryFile::GetBondAngleMap() const
    {
      s_ReadingTimer.Start();
      if( m_Lines.IsEmpty())
      {
        LoadFiles();
      }
      std::string fname( m_Filename + ".bond_angles.txt" + GetFileExtension());
      if( !io::DirectoryEntry( fname).DoesExist())
      {
        BCL_MessageCrt( "Bond angle data is unavailable! Searched: " + fname);
        s_ReadingTimer.Stop();
        return t_BondAngleMap();
      }
      t_BondAngleMap map;
      io::IFStream input;
      io::File::MustOpenIFStream( input, fname);
      while( input.good())
      {
        std::string atombondline;
        std::getline( input, atombondline);
        if( !input.good())
        {
          break;
        }
        atombondline = util::TrimString( atombondline);
        // check for empty line or end of file
        while( atombondline.empty())
        {
          std::getline( input, atombondline);
          atombondline = util::TrimString( atombondline);
          if( !input.good())
          {
            break;
          }
        }
        if( !input.good())
        {
          break;
        }
        storage::Vector< std::string> atombondlist( util::SplitString( atombondline));
        storage::Triplet
        <
          ConformationGraphConverter::AtomComparisonTypeEnum,
          size_t,
          storage::Vector< storage::Pair< size_t, size_t> >
        > key;
        key.First() = ConformationGraphConverter::AtomComparisonTypeEnum( atombondlist( 1));
        key.Second() = AtomType( atombondlist( 0)).GetIndex();
        const size_t nbonds( ( atombondlist.GetSize() - 2) / 2);
        key.Third().Resize( nbonds);
        for( size_t bondn( 0), pos( 2); bondn < nbonds; ++bondn, pos += 2)
        {
          key.Third()( bondn).First() = util::ConvertStringToNumericalValue< size_t>( atombondlist( pos));
          key.Third()( bondn).Second() = key.First() == ConformationGraphConverter::e_AtomType
                                         ? AtomType( atombondlist( pos + 1)).GetIndex()
                                         : key.First() == ConformationGraphConverter::e_ElementType
                                           ? ElementType( atombondlist( pos + 1)).GetIndex()
                                           : size_t( 1); // identity
        }
        io::Serialize::Read( map[ key], input);
        auto &mapkey( map[ key]);
        double max_cnt( 0.0);
        for( auto itr_mk( mapkey.Begin()), itr_mk_end( mapkey.End()); itr_mk != itr_mk_end; ++itr_mk)
        {
          max_cnt = std::max( itr_mk->Second().GetWeight(), max_cnt);
        }
        double threshold( max_cnt / 100.0);
        for( auto itr_mk( mapkey.Begin()), itr_mk_end( mapkey.End()); itr_mk != itr_mk_end;)
        {
          if( itr_mk->Second().GetWeight() > threshold)
          {
            ++itr_mk;
          }
          else
          {
            auto itr_mk_cpy( itr_mk);
            ++itr_mk;
            mapkey.RemoveElement( itr_mk_cpy);
          }
        }
      }
      io::File::CloseClearFStream( input);
      s_ReadingTimer.Stop();
      return map;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RotamerLibraryFile::GetSerializer() const
    {
      io::Serializer parameters;
      if( m_Alias.empty())
      {
        // for the storage, the file may or may not already exist, but the extension should be valid
        parameters.AddInitializer
        (
          "prefix",
          "File path and prefix to the rotamer library. If not provided, uses the default rotamer library, located in "
          "the bcl directory that is home to this executable",
          io::Serialization::GetAgent( &m_Filename),
          ""
        );
        parameters.AddInitializer
        (
          "number_of_files",
          "Split the rotamer library into directories containing files. Each file contains 100 rotamers. Set the number of files ",
          io::Serialization::GetAgent( &m_FileNumber),
          "100"
        );
        parameters.AddInitializer
        (
          "compression",
          "compression to use for files",
          io::Serialization::GetAgent( &m_Compression),
          io::GetStreamBufferClasses().HaveEnumWithName( "GZ") ? "GZ" : "Uncompressed"
        );
      }
      else
      {
        parameters.SetClassDescription( "Load fragments from the " + m_Alias + " library");
      }

      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool RotamerLibraryFile::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      if( m_Filename.empty())
      {
        m_Filename = GetRotamerFinder().FindFile( "") + GetPrefix();
        GetRotlibLicenseText();
      }
      return true;
    }

    //! @brief Load in necessary information, if needed
    void RotamerLibraryFile::LoadFiles() const
    {
      if( m_Filename.empty())
      {
        m_Filename = GetRotamerFinder().FindFile( "") + GetPrefix();
      }
      if( io::DirectoryEntry( m_Filename + ".constitutions.txt" + GetFileExtension()).DoesExist())
      {
        // IFstream for reading in source sdf file containing constitutions
        io::IFStream read;
        // open sdf file for reading
        io::File::MustOpenIFStream( read, m_Filename + ".constitutions.txt" + GetFileExtension());

        // split the file into a format so that it can be iterated
        m_Lines = util::SplittedStringLineListFromIStream( read, " ");

        io::File::CloseClearFStream( read);
      }
    }

    //! @brief Load in necessary information, if needed
    void RotamerLibraryFile::GetRotlibLicenseText() const
    {
      static bool s_DisplayMessage( false);

      if( !command::CommandState::GetInMainCommandLineParsing() && !s_DisplayMessage)
      {
        // IFstream for reading in source sdf file containing constitutions
        io::IFStream read;

        std::string license_file( GetRotamerFinder().FindFile( "") + GetPrefix() + ".license.txt");

        if( io::DirectoryEntry( license_file).DoesExist())
        {
          // open sdf file for reading
          io::File::MustOpenIFStream( read, license_file);
          std::string line;
          std::string file;
          while( std::getline( read, line))
          {
            file += '\n' + line;
          }
          util::GetLogger() << file << '\n' << std::endl;
          s_DisplayMessage = true;
        }
      }
      return;
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace chemistry
} // namespace bcl
