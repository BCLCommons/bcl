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
#include "descriptor/bcl_descriptor_aa_dssp_info.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_statistics.h"
#include "sspred/bcl_sspred_method_handler.h"
#include "util/bcl_util_string_functions.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    //! @brief DSSPInfoType as string
    //! @param DSSP_INFO_TYPE the message level
    //! @return the DSSPInfoType as string
    const std::string &AADSSPInfo::GetDSSPInfoTypeString( const DSSPInfoType &DSSP_INFO_TYPE)
    {
      static const std::string s_types[] =
      {
        "DSSP_Accessibility",
        "DSSP_TotalEnergy",
        "DSSP_MaxHBondPotential",
        "DSSP_MaxHBondNeighborOffset",
        "DSSP_IsInSheetCore",
        "DSSP_IsOnSheetEdge",
        "DSSP_NonBondedStrandResidue",
        "DSSP_BondedStrandResidue",
        "DSSP_BondsParallel",
        "DSSP_BondsAntiparallel",
        "DSSP_BondsBridge",
        GetStaticClassName< DSSPInfoType>()
      };
      return s_types[ DSSP_INFO_TYPE];
    }

    //! @brief DSSPInfoType description
    //! @param DSSP_INFO_TYPE the description
    //! @return the DSSPInfoType's description
    const std::string &AADSSPInfo::GetDSSPInfoTypeDescription( const DSSPInfoType &SSE_INFO_TYPE)
    {
      static const std::string s_types[] =
      {
        "Residue water exposed surface area in A^2",
        "Total h-bond energy determined for this AA",
        "Best (lowest) energy for any neighbor",
        "Offset of the neighbor with max hbond energy",
        "1 if residue is h-bonded to another strand and at least one neighboring residue is too",
        "1 if the residue is labeled a strand but neither neighbor is bonded to a strand, though this residue is",
        "1 if the residue is labeled a strand but this residue is not hbonded to a strand",
        "1 if the residue is labeled a strand and this residue is hbonded to a strand",
        "True if the residue bonds to another strand oriented parallel to this strand",
        "True if the residue bonds to another strand oriented antiparallel to this strand",
        "True if the residue bonds to another strand only at this residue",
        GetStaticClassName< AADSSPInfo::DSSPInfoType>()
      };

      return s_types[ SSE_INFO_TYPE];
    }

    namespace
    {
      util::SiPtr< const util::ObjectInterface> AddInstances()
      {
        util::SiPtr< const util::ObjectInterface> last_instance;
        for( AADSSPInfo::DSSPInfoTypeEnum i; i < AADSSPInfo::s_NumberDSSPInfoTypes; ++i)
        {
          last_instance = util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AADSSPInfo( i));
        }
        return last_instance;
      }
    }

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> AADSSPInfo::s_Instance
    (
      AddInstances()
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AADSSPInfo::AADSSPInfo( const DSSPInfoTypeEnum &INFO_TYPE_ENUM) :
      m_InfoTypeEnum( INFO_TYPE_ENUM),
      m_DsspExtension( ( *sspred::GetMethods().e_DSSP)->GetFileExtension())
    {
    }

    //! @brief Clone function
    //! @return pointer to new BaseElement
    AADSSPInfo *AADSSPInfo::Clone() const
    {
      return new AADSSPInfo( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AADSSPInfo::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &AADSSPInfo::GetAlias() const
    {
      return m_InfoTypeEnum.GetString();
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t AADSSPInfo::GetNormalSizeOfFeatures() const
    {
      return 1;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT_A the element of interest
    //! @param STORAGE storage for the descriptor
    void AADSSPInfo::Calculate
    (
      const iterate::Generic< const biol::AABase> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      // if the sequence maps are empty, this indicates that this object is now operating over a new sequence,
      // so it is necessary to reload the files
      if( m_Info.IsEmpty())
      {
        LoadFiles();
      }

      const int seq_id( ELEMENT->GetData()->GetSeqID());
      storage::Map< int, size_t>::const_iterator map_itr( m_SeqIDToIndexMap.Find( seq_id));

      if( map_itr == m_SeqIDToIndexMap.End())
      {
        STORAGE( 0) = util::GetUndefined< float>();
        return;
      }

      const size_t vec_index( map_itr->second);
      STORAGE( 0) = m_Info( vec_index);

      return;
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void AADSSPInfo::SetObjectHook()
    {
      m_SeqIDToIndexMap.Reset();
      m_Info = linal::Vector< float>();
    }

    //! @brief function to load files; should only be called the first time Calculate is called with a new sequence
    //! since the results are often in the cache
    void AADSSPInfo::LoadFiles()
    {
      // Get protein model
      util::SiPtr< const assemble::ProteinModelWithCache> sp_protein_model( this->GetCurrentObject());

      // Check to make sure there is only a single chain in the given PDB file
      BCL_Assert
      (
        sp_protein_model.IsDefined() && sp_protein_model->GetChains().GetSize() == 1,
        "More than one chain in the given protein model: "
         + util::Format()( sp_protein_model->GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile))
      );

      // Get the filename
      util::ShPtr< util::Wrapper< std::string> > sp_filename_wrapper
      (
        sp_protein_model->GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile)
      );
      const std::string pdb_filenamebase( io::File::RemoveLastExtension( io::File::RemoveCompressionExtension( sp_filename_wrapper->GetData())));
      const std::string pdb_filenamedirectory( io::File::SplitToPathAndFileName( pdb_filenamebase).First());
      std::string pdb_filename_nodirectory( io::File::SplitToPathAndFileName( pdb_filenamebase).Second());

      // now, guess the dssp filename.  If the base of the filename without extension is 4 or fewer letters, then just
      // look for the filename base
      const char chainid( sp_protein_model->GetChainIDs()[ 0]);
      storage::Vector< io::DirectoryEntry> dssp_files;
      if( pdb_filename_nodirectory.size() > size_t( 4))
      {
        if( pdb_filename_nodirectory[ pdb_filename_nodirectory.size() - 1] == chainid)
        {
          pdb_filename_nodirectory.erase( pdb_filename_nodirectory.size() - 1);
          dssp_files =
            sspred::MethodHandler::PossibleFileEntries
            (
              m_DsspExtension,
              chainid,
              pdb_filenamedirectory + pdb_filename_nodirectory,
              ""
            );
          pdb_filename_nodirectory += chainid;
        }
        if
        (
          pdb_filename_nodirectory.size() != size_t( 5)
          && pdb_filename_nodirectory[ 4] == chainid && dssp_files.IsEmpty()
        )
        {
          std::string pdb_code( pdb_filename_nodirectory.substr( 0, 4));
          dssp_files.Append
          (
            sspred::MethodHandler::PossibleFileEntries
            (
              m_DsspExtension,
              chainid,
              pdb_filenamedirectory + pdb_code,
              ""
            )
          );
        }
      }
      if( dssp_files.IsEmpty())
      {
        dssp_files.Append
        (
          sspred::MethodHandler::PossibleFileEntries
          (
            m_DsspExtension,
            chainid,
            pdb_filenamedirectory + pdb_filename_nodirectory,
            ""
          )
        );
      }

      // read through the dssp file until the proper chain is found
      if( dssp_files.IsEmpty())
      {
        BCL_MessageCrt( "No dssp file found for " + pdb_filenamebase + "; assuming it has no secondary structure");
        m_Info = linal::Vector< float>( size_t( 1), float( 0.0));
        return;
      }

      // initialize necessary variables
      std::string line;
      const char sequence_chain_id( chainid);

      // track all chain ids seen in dssp that are not == the sequence chain id
      std::string dssp_chain_ids;

      // track the number of ss types found for the chain
      size_t number_aas_found( 0);

      io::IFStream ifstream;
      io::File::MustOpenIFStream( ifstream, dssp_files.FirstElement().GetFullName());

      // read through the file until the residue lines are reached
      while( ifstream.good())
      {
        std::getline( ifstream, line);

        // does the line contain SS summary information?
        if( util::StartsWith( line, "  #  RESIDUE"))
        {
          // Reached residue lines
          break;
        }
      }

      // keep track of all current data
      storage::Vector< float> current_data;
      storage::Vector< int>   max_energy_neighbor_offset;
      storage::Vector< int>   max_energy_neighbor_offset_alt;
      std::string             ss_type_str;
      std::string             is_h_bonded;

      // read through residue lines
      while( ifstream.good())
      {
        std::getline( ifstream, line);

        // skip lines that are too short
        if( line.size() < 17)
        {
          continue;
        }

        // the split string contains items in the following order
        // characters 0-4: Residue ID
        // characters 5-9: PDB_ID
        // character   11: Chain ID
        // character   13: 1-letter AA code
        // character   16: 1-letter SS type
        // followed by a bunch of irrelevant (for our purposes) information

        const char ss_type_chr( line[ 16]);
        const std::string pdb_id_str( line.substr( 5, 5));

        if( pdb_id_str[ 4] == ' ' && pdb_id_str[ 3] == ' ')
        {
          // empty pdb id, happens with gaps, continue
          continue;
        }

        // get the chain id.  Dssp automatically changes blank chain ids to - for easier parsing, so map them back
        const char chain_id( line[ 11] == '-' ? ' ' : line[ 11]);

        // skip undesired chain ids
        if( chain_id != sequence_chain_id)
        {
          if( dssp_chain_ids.find( chain_id) == std::string::npos)
          {
            // new chain id
            dssp_chain_ids += chain_id;
          }
          continue;
        }

        // check that the pdb ids are numeric
        int pdb_id( util::GetUndefined< int>());
        if( !util::TryConvertFromString( pdb_id, pdb_id_str, util::GetLogger()))
        {
          BCL_MessageCrt( "Non-numeric PDB id in " + line);
          continue;
        }
        m_SeqIDToIndexMap[ pdb_id] = number_aas_found++;
        if( m_InfoTypeEnum == e_Accessibility)
        {
          // handle the trivial case where only the accessibility is required
          current_data.PushBack( util::ConvertStringToNumericalValue< float>( line.substr( 35, 3)));
          continue;
        }

        // determine SS-type
        const bool ss_type_is_strand( ss_type_chr == 'E');

        // find the start of the neighbor indice / energy list; normally this is at 45, but it can be a few later if
        // there are insertion codes
        size_t n_start( line.find( ',', 44) - 6);

        // get all neighbors and energies
        storage::Vector< int> neighbor_offset;
        storage::Vector< float> neighbor_energy;
        for
        (
          size_t neighbor_number( 0), max_neighbors( 4);
          neighbor_number < max_neighbors;
          ++neighbor_number, n_start += 11
        )
        {
          const int neighbor( util::ConvertStringToNumericalValue< int>( line.substr( n_start, 6)));
          const float energy( util::ConvertStringToNumericalValue< float>( line.substr( n_start + 7, 4)));
          size_t neighbor_offset_id( neighbor_offset.Find( neighbor));
          if( neighbor_offset_id < neighbor_offset.GetSize())
          {
            neighbor_energy( neighbor_offset_id) += energy;
          }
          else
          {
            neighbor_offset_id = neighbor_offset.GetSize();
            neighbor_offset.PushBack( neighbor);
            neighbor_energy.PushBack( energy);
          }
        }
        const float total_e( math::Statistics::Sum( neighbor_energy.Begin(), neighbor_energy.End()));
        if( m_InfoTypeEnum == e_TotalEnergy)
        {
          // just get the sum of energies
          current_data.PushBack( total_e);
          continue;
        }
        std::vector< std::pair< float, int> > sorted_neighbor_energies;
        for
        (
          size_t neighbor_number( 0), number_neighbors( neighbor_offset.GetSize());
          neighbor_number < number_neighbors;
          ++neighbor_number
        )
        {
          sorted_neighbor_energies.push_back
          (
            std::make_pair( neighbor_energy( neighbor_number), neighbor_offset( neighbor_number))
          );
        }
        std::sort( sorted_neighbor_energies.begin(), sorted_neighbor_energies.end());
//        for( size_t i( 0); i < sorted_neighbor_energies.size(); ++i)
//        {
//          BCL_MessageStd
//          (
//            "pdb id: " + util::Format()( pdb_id)
//            + " " + util::Format()( i) + "th highest energy: "
//            + util::Format()( sorted_neighbor_energies[ i].first)
//            + " @ offset: " + util::Format()( sorted_neighbor_energies[ i].second)
//          );
//        }
        if( m_InfoTypeEnum == e_MaxHBondPotential)
        {
          // just get the min energy
          current_data.PushBack( sorted_neighbor_energies.begin()->first);
          continue;
        }
        else if( m_InfoTypeEnum == e_MaxHBondNeighborOffset)
        {
          // just get the min energy neighbor
          current_data.PushBack( sorted_neighbor_energies.begin()->second);
          continue;
        }

        const bool is_making_hbond( total_e < -1.7);
        if( m_InfoTypeEnum == e_BondedStrandResidue)
        {
          current_data.PushBack( is_making_hbond && ss_type_is_strand ? 1.0 : 0.0);
          continue;
        }
        else if( m_InfoTypeEnum == e_NonBondedStrandResidue)
        {
          current_data.PushBack( !is_making_hbond && ss_type_is_strand ? 1.0 : 0.0);
          continue;
        }
        else if( m_InfoTypeEnum == e_IsInSheetCore)
        {
          current_data.PushBack( 0.0);
          if
          (
            ss_type_is_strand
            && is_making_hbond
            && is_h_bonded.size()
            && is_h_bonded[ is_h_bonded.size() - 1] == '1'
          )
          {
            current_data( current_data.GetSize() - 1) = current_data( current_data.GetSize() - 2) = 1.0;
          }
        }
        else if( m_InfoTypeEnum == e_IsOnSheetEdge)
        {
          current_data.PushBack( ss_type_is_strand && is_making_hbond ? 1.0 : 0.0);
          if
          (
            current_data.LastElement() > 0.5
            && is_h_bonded.size()
            && is_h_bonded[ is_h_bonded.size() - 1] == '1'
          )
          {
            current_data( current_data.GetSize() - 1) = current_data( current_data.GetSize() - 2) = 0.0;
          }
        }
        is_h_bonded += is_making_hbond && ss_type_is_strand ? '1' : '0';
        ss_type_str += ss_type_chr;
        if( m_InfoTypeEnum == e_IsInSheetCore || m_InfoTypeEnum == e_IsOnSheetEdge)
        {
          continue;
        }
        // only types left are e_BondsParallel, e_BondsAntiparallel, e_BondsBridge, all of which need the bonding pattern
        if( is_making_hbond)
        {
          max_energy_neighbor_offset.PushBack( sorted_neighbor_energies[ 0].second);
        }
        else
        {
          max_energy_neighbor_offset.PushBack( 0);
        }
        if( sorted_neighbor_energies.size() > size_t( 1))
        {
          const std::pair< float, int> second_best_energy( sorted_neighbor_energies[ 1]);
          if( second_best_energy.first < -1.0)
          {
            max_energy_neighbor_offset_alt.PushBack( second_best_energy.second);
          }
          else
          {
            max_energy_neighbor_offset_alt.PushBack( 0);
          }
        }
        else
        {
          max_energy_neighbor_offset_alt.PushBack( 0);
        }
      }

      io::File::CloseClearFStream( ifstream);

      // discern the bonding pattern.  for parallel, allow no more than 1 residue difference
      // for anti-parallel, allow 3-5 residues difference when moving two residues away
      // all others -- classify as bridge
      if( m_InfoTypeEnum == e_BondsParallel || m_InfoTypeEnum == e_BondsAntiparallel || m_InfoTypeEnum == e_BondsBridge)
      {
        for( size_t res_id( 0), n_res( max_energy_neighbor_offset.GetSize()); res_id < n_res; ++res_id)
        {
          if( ss_type_str[ res_id] != 'E' || is_h_bonded[ res_id] != '1')
          {
            current_data.PushBack( 0.0);
            continue;
          }
          const int neighbor_a( max_energy_neighbor_offset( res_id));
          const int neighbor_b( max_energy_neighbor_offset_alt( res_id));
          size_t parallel_neighbor_counts( 0), antiparallel_neighbor_counts( 0);
          {
            storage::Vector< int> residues_to_check;
            // determine distance till the end of the strand
            const size_t strand_right_length
            (
              std::min( ss_type_str.size(), ss_type_str.find_first_not_of( 'E', res_id) - size_t( res_id + 1))
            );
            // determine distance till beginning of the strand; or 3, whichever is lower
            size_t strand_left_length( 0);
            for
            (
              int res_id_prev( res_id - 1);
              res_id_prev >= 0 && strand_left_length < 3 && ss_type_str[ res_id_prev] == 'E';
              ++strand_left_length, --res_id_prev
            )
            {
            }
            // order the residues to check by distance from the non-twisted strand alignment, in which residues that
            // connect to the same strand have a single residue in between them
            // When strands kink or do not align perfectly, there may be 2 residues between them instead, so check
            // these afterwards
            if( strand_left_length >= 2)
            {
              residues_to_check.PushBack( res_id - 2);
            }
            if( strand_right_length >= 2)
            {
              residues_to_check.PushBack( res_id + 2);
            }
            if( strand_left_length >= 3)
            {
              residues_to_check.PushBack( res_id - 3);
            }
            if( strand_right_length >= 3)
            {
              residues_to_check.PushBack( res_id + 3);
            }
            // test each of the residues in the order given
            for
            (
              storage::Vector< int>::const_iterator
                itr_neighbor( residues_to_check.Begin()), itr_neighbor_end( residues_to_check.End());
              itr_neighbor != itr_neighbor_end;
              ++itr_neighbor
            )
            {
              // skip residues that are not h-bonded
              if( is_h_bonded[ *itr_neighbor] != '1')
              {
                continue;
              }
              // determine strand partners neighbor offsets
              const int neighbor_a_prev( max_energy_neighbor_offset( *itr_neighbor));
              const int neighbor_b_prev( max_energy_neighbor_offset_alt( *itr_neighbor));

              // collect all valid strand neighbor offset differences between the central residue and this residue
              // to a vector
              storage::Vector< int> nearby_offsets;
              nearby_offsets.PushBack( math::Absolute( neighbor_a - neighbor_a_prev));
              if( neighbor_b)
              {
                nearby_offsets.PushBack( math::Absolute( neighbor_b - neighbor_a_prev));
                if( neighbor_b_prev)
                {
                  nearby_offsets.PushBack( math::Absolute( neighbor_a - neighbor_b_prev));
                  nearby_offsets.PushBack( math::Absolute( neighbor_b - neighbor_b_prev));
                }
              }
              else if( neighbor_b_prev)
              {
                nearby_offsets.PushBack( math::Absolute( neighbor_a - neighbor_b_prev));
              }
              for
              (
                storage::Vector< int>::const_iterator
                  itr_offsets( nearby_offsets.Begin()), itr_offsets_end( nearby_offsets.End());
                itr_offsets != itr_offsets_end;
                ++itr_offsets
              )
              {
                if( *itr_offsets < 2) // normally, parallel residues will have an offset of 0, but 1 can happen if there is minor twist
                {
                  // increment parallel neighbor
                  ++parallel_neighbor_counts;
                }
                else if( *itr_offsets > 2 && *itr_offsets < 6)
                {
                  // normally, antiparallel residues will have an offset of 4, but 3-5 can happen if there is minor twist or bulge
                  // increment anti-parallel neighbor
                  ++antiparallel_neighbor_counts;
                }
              }
              // stop if this residue is in either a parallel or antiparallel configuration
              if( parallel_neighbor_counts != antiparallel_neighbor_counts)
              {
                break;
              }
            }
          }
          // test whether the condition was met
          if
          (
            ( parallel_neighbor_counts > antiparallel_neighbor_counts && m_InfoTypeEnum == e_BondsParallel)
            || ( antiparallel_neighbor_counts > parallel_neighbor_counts && m_InfoTypeEnum == e_BondsAntiparallel)
            || ( parallel_neighbor_counts == antiparallel_neighbor_counts && m_InfoTypeEnum == e_BondsBridge)
          )
          {
            current_data.PushBack( 1.0);
          }
          else
          {
            current_data.PushBack( 0.0);
          }
        }
      }

      if
      (
        ( m_InfoTypeEnum == e_IsInSheetCore || m_InfoTypeEnum == e_IsOnSheetEdge)
        && is_h_bonded.size() >= size_t( 2)
      )
      {
        // the simple windowing approach used in the previous for-loop fails on bridges of size 2, so manually fix these
        size_t bridge_pos( 0);
        if( is_h_bonded[ 0] == '1' && is_h_bonded[ 1] == '1' && is_h_bonded[ 2] == '0')
        {
          bridge_pos = 0;
        }
        else
        {
          bridge_pos = is_h_bonded.find( "0110");
        }
        const float new_data( m_InfoTypeEnum == e_IsInSheetCore ? 0.0 : 1.0);
        while( bridge_pos < is_h_bonded.size())
        {
          current_data( bridge_pos + 1) = new_data;
          current_data( bridge_pos + 2) = new_data;
          bridge_pos = is_h_bonded.find( "0110", bridge_pos + 1);
        }
        if
        (
          is_h_bonded[ is_h_bonded.size() - 1] == '1'
          && is_h_bonded[ is_h_bonded.size() - 2] == '1'
          && is_h_bonded[ is_h_bonded.size() - 3] == '0'
        )
        {
          current_data( is_h_bonded.size() - 1) = new_data;
          current_data( is_h_bonded.size() - 2) = new_data;
        }
      }

      m_Info = current_data;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AADSSPInfo::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        GetDSSPInfoTypeDescription( m_InfoTypeEnum) + " based on DSSP output. "
      );

      parameters.AddOptionalInitializer
      (
        "extension",
        "extension of the dssp file, defaults to .dssp",
        io::Serialization::GetAgent( &m_DsspExtension)
      );

      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
