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
#include "internal/biology/bcl_app_generate_aa_pair_statistics.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_with_cache_storage_file.h"
#include "biol/bcl_biol_aa_classes.h"
#include "biol/bcl_biol_atom.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math.h"
#include "storage/bcl_storage_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {

    const ApplicationType GenerateAAPairStatistics::GenerateAAPairStatistics_Instance
    (
      GetAppGroups().AddAppToGroup( new GenerateAAPairStatistics(), GetAppGroups().e_InternalBiol)
    );

  ///////////////////////////////////
  // construction and destruction //
  ///////////////////////////////////

    //! @brief standard constructor
    GenerateAAPairStatistics::GenerateAAPairStatistics() :
      m_OutputFilenameBase
      (
        new command::FlagStatic
        (
          "output",
          "base name for output of histograms",
          command::Parameter( "output", "base name for output of histograms")
        )
      ),
      m_MaxDistanceFlag
      (
        new command::FlagStatic
        (
          "max_distance",
          "Maximum distance between residues to consider",
          command::Parameter( "distance", "", command::ParameterCheckRanged< size_t>( 1, 100))
        )
      ),
      m_InputFlag
      (
        new command::FlagDynamic
        (
          "input",
          "input flag",
          command::Parameter( "input", "", command::ParameterCheckSerializable( assemble::ProteinWithCacheStorageFile( true, true))),
          1
        )
      ),
      m_SSPredFlag
      (
        new command::FlagStatic
        (
          "native_method",
          "Method to use to determine the native SSEs",
          command::Parameter( "native_method", "", command::ParameterCheckSerializable( sspred::Method()))
        )
      ),
      m_ConsiderEnvironmentFlag
      (
        new command::FlagStatic
        (
          "environment",
          "Use this flag to calculate pair statistics for the environment as well as the ss types"
        )
      ),
      m_RestrictEnvironmentFlag
      (
        new command::FlagStatic
        (
          "restrict_environment",
          "Use this flag to calculate statistics for residues only in a given environment",
          command::Parameter
          (
            "environment",
            "environment to take statistics for",
            command::ParameterCheckSerializable( biol::EnvironmentType()),
            ""
          )
        )
      )
    {
    }

    //! @brief Clone function
    //! @return pointer to new GenerateAAPairStatistics
    GenerateAAPairStatistics *GenerateAAPairStatistics::Clone() const
    {
      return new GenerateAAPairStatistics( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &GenerateAAPairStatistics::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> GenerateAAPairStatistics::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

    ////////////////////
    // common options //
    ////////////////////

      // add flags for input / output
      sp_cmd->AddFlag( m_InputFlag);

      // Output filename base
      sp_cmd->AddFlag( m_OutputFilenameBase);

      // SSPred type that should be used to determine the real statistics
      sp_cmd->AddFlag( m_SSPredFlag);

      // max distance
      sp_cmd->AddFlag( m_MaxDistanceFlag);

      // environment
      sp_cmd->AddFlag( m_ConsiderEnvironmentFlag);

      // restricted environment
      sp_cmd->AddFlag( m_RestrictEnvironmentFlag);

    ///////////////////
    // default flags //
    ///////////////////

      // default flags are unnecessary for this application, but message level and read-me are useful
      command::GetAppDefaultFlags().AddRequiredCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief compute the statistics requested
    //! @param MODEL the model of interest
    void GenerateAAPairStatistics::Process( const assemble::ProteinModelWithCache &MODEL) const
    {
      // generate vectors for the SS type and aa types for this protein
      const size_t protein_size( MODEL.GetSize());
      storage::Vector< biol::SSType> ss_types( protein_size, biol::GetSSTypes().COIL);
      storage::Vector< biol::EnvironmentType> env_types( protein_size, biol::GetEnvironmentTypes().e_Solution);
      linal::Vector< size_t> aa_type_ids_first( protein_size, util::GetUndefined< size_t>());
      linal::Vector< size_t> aa_type_ids_second( protein_size, util::GetUndefined< size_t>());
      linal::Vector< int> pdb_ids( protein_size, util::GetUndefined< int>());
      storage::Vector< linal::Vector3D> aa_positions( protein_size);

      // iterate over all residues
      size_t position( 0);
      for
      (
        iterate::Generic< const biol::AABase> aa_iterator( MODEL.GetIterator());
        aa_iterator.NotAtEnd();
        ++aa_iterator, ++position
      )
      {
        // get the AA type
        biol::AAType natural_aa_type( aa_iterator->GetType());
        if( natural_aa_type.IsDefined() && !natural_aa_type->IsNaturalAminoAcid())
        {
          natural_aa_type = natural_aa_type->GetParentType();
        }
        // get index
        const size_t natural_index
        (
          natural_aa_type.IsDefined() && natural_aa_type->IsNaturalAminoAcid()
          ? natural_aa_type.GetIndex()
          : util::GetUndefined< size_t>()
        );
        if( util::IsDefined( natural_index))
        {
          aa_type_ids_first( position) = natural_index * ( biol::AATypes::s_NumberStandardAATypes + 1);
          aa_type_ids_second( position) = natural_index;
        }
        pdb_ids( position) = aa_iterator->GetPdbID();

        // get the ss type
        util::SiPtr< const sspred::MethodInterface> method_prediction( aa_iterator->GetSSPrediction( m_Method));
        if( method_prediction.IsDefined())
        {
          ss_types( position) = method_prediction->GetOneStateSSPrediction();
        }
        else
        {
          BCL_MessageCrt( "No ss type found for residue: " + aa_iterator->GetIdentification());
        }
        if( m_ConsiderEnvironment || m_EnvironmentType.IsDefined())
        {
          util::SiPtr< const sspred::MethodInterface> actual_environment( aa_iterator->GetSSPrediction( sspred::GetMethods().e_PDB));
          if( actual_environment.IsDefined())
          {
            env_types( position) = actual_environment->GetOneStateTMPrediction();
          }
          else
          {
            BCL_MessageCrt( "No environment type found for residue: " + aa_iterator->GetIdentification());
          }
        }
        aa_positions( position) = aa_iterator->GetAtom( biol::GetAtomTypes().N).GetCoordinates();
      }

      // iterate through the sspredictions, accumulate all relevant data
      for( size_t distance_size( 0); distance_size < m_MaxDistance; ++distance_size)
      {
        // get references to the statistics vectors
        linal::Vector< size_t> &strand_strand_counts( m_StrandStrandCounts( distance_size));
        linal::Vector< size_t> &helix_helix_counts( m_HelixHelixCounts( distance_size));
        linal::Vector< size_t> &coil_coil_counts( m_CoilCoilCounts( distance_size));
        linal::Vector< size_t> &start_helix_counts( m_StartsHelixCounts( distance_size));
        linal::Vector< size_t> &start_strand_counts( m_StartsStrandCounts( distance_size));
        linal::Vector< size_t> &start_coil_counts( m_StartsCoilCounts( distance_size));
        linal::Vector< size_t> &end_helix_counts( m_EndsHelixCounts( distance_size));
        linal::Vector< size_t> &end_strand_counts( m_EndsStrandCounts( distance_size));
        linal::Vector< size_t> &end_coil_counts( m_EndsCoilCounts( distance_size));
        linal::Vector< size_t> &aa_pair_counts( m_AAPairTypeCounts( distance_size));

        linal::Vector< size_t> &trans_trans_counts( m_TransitionTransitionCounts( distance_size));
        linal::Vector< size_t> &membrane_membrane_counts( m_MembraneMembraneCounts( distance_size));
        linal::Vector< size_t> &sol_sol_counts( m_SolubleSolubleCounts( distance_size));
        linal::Vector< size_t> &start_membrane_counts( m_StartsMembraneCounts( distance_size));
        linal::Vector< size_t> &start_trans_counts( m_StartsTransitionCounts( distance_size));
        linal::Vector< size_t> &start_sol_counts( m_StartsSolubleCounts( distance_size));
        linal::Vector< size_t> &end_membrane_counts( m_EndsMembraneCounts( distance_size));
        linal::Vector< size_t> &end_trans_counts( m_EndsTransitionCounts( distance_size));
        linal::Vector< size_t> &end_sol_counts( m_EndsSolubleCounts( distance_size));

        storage::Vector< math::RunningAverage< double> > &separations( m_Separation( distance_size));

        for
        (
          size_t position_first( 0), position_second( distance_size + 1);
          position_second < protein_size;
          ++position_first, ++position_second
        )
        {
          const size_t id_first( aa_type_ids_first( position_first)), id_second( aa_type_ids_second( position_second));
          if( !util::IsDefined( id_first) || !util::IsDefined( id_second))
          {
            continue;
          }
          int target_pdb_id( pdb_ids( position_first) + distance_size + 1);
          size_t position_sec( position_second);
          while( position_sec > position_first && pdb_ids( position_sec) > target_pdb_id)
          {
            --position_sec;
          }
          if( pdb_ids( position_sec) != target_pdb_id)
          {
            continue;
          }

          bool first_in_correct_env( true);
          bool second_in_correct_env( true);

          // handle environment restriction
          if( m_EnvironmentType.IsDefined())
          {
            if( env_types( position_first) != m_EnvironmentType)
            {
              first_in_correct_env = false;
            }
            if( env_types( position_second) != m_EnvironmentType)
            {
              second_in_correct_env = false;
            }
            if( !first_in_correct_env && !second_in_correct_env)
            {
              continue;
            }
          }

          // get the basic aa-type index
          const size_t aa_type_index( id_first + id_second);

          ++aa_pair_counts( aa_type_index);
          const double distance( linal::Distance( aa_positions( position_first), aa_positions( position_sec)));
          if( !util::IsDefined( distance))
          {
            BCL_MessageStd
            (
              "Undefined distance between aa positions "
              + util::Format()( position_first) + " and " + util::Format()( position_sec)
            );
          }
          else
          {
            separations( aa_type_index) += distance;
          }
          if( ss_types( position_first) == ss_types( position_sec))
          {
            bool is_start_end( false);
            biol::SSType nominal_type( ss_types( position_sec));
            for( size_t p( position_first + 1); p < position_sec; ++p)
            {
              if( ss_types( p) != nominal_type)
              {
                is_start_end = true;
                break;
              }
            }
            if( nominal_type == biol::GetSSTypes().HELIX)
            {
              if( is_start_end)
              {
                if( second_in_correct_env)
                {
                  ++start_helix_counts( aa_type_index);
                }
                if( first_in_correct_env)
                {
                  ++end_helix_counts( aa_type_index);
                }
              }
              else
              {
                // two helices
                ++helix_helix_counts( aa_type_index);
              }
            }
            else if( nominal_type == biol::GetSSTypes().STRAND)
            {
              if( is_start_end)
              {
                if( second_in_correct_env)
                {
                  ++start_strand_counts( aa_type_index);
                }
                if( first_in_correct_env)
                {
                  ++end_strand_counts( aa_type_index);
                }
              }
              else
              {
                // two strands
                ++strand_strand_counts( aa_type_index);
              }
            }
            else if( nominal_type == biol::GetSSTypes().COIL)
            {
              if( is_start_end)
              {
                if( second_in_correct_env)
                {
                  ++start_coil_counts( aa_type_index);
                }
                if( first_in_correct_env)
                {
                  ++end_coil_counts( aa_type_index);
                }
              }
              else
              {
                // two coils
                ++coil_coil_counts( aa_type_index);
              }
            }
          }
          else if( ss_types( position_first) == biol::GetSSTypes().STRAND)
          {
            if( first_in_correct_env)
            {
              ++end_strand_counts( aa_type_index);
            }
            if( second_in_correct_env)
            {
              if( ss_types( position_sec) == biol::GetSSTypes().HELIX)
              {
                ++start_helix_counts( aa_type_index);
              }
              else
              {
                ++start_coil_counts( aa_type_index);
              }
            }
          }
          else if( ss_types( position_first) == biol::GetSSTypes().HELIX)
          {
            if( first_in_correct_env)
            {
              ++end_helix_counts( aa_type_index);
            }
            if( second_in_correct_env)
            {
              if( ss_types( position_sec) == biol::GetSSTypes().STRAND)
              {
                ++start_strand_counts( aa_type_index);
              }
              else
              {
                ++start_coil_counts( aa_type_index);
              }
            }
          }
          else
          {
            // first is coil
            if( first_in_correct_env)
            {
              ++end_coil_counts( aa_type_index);
            }
            if( second_in_correct_env)
            {
              if( ss_types( position_sec) == biol::GetSSTypes().STRAND)
              {
                ++start_strand_counts( aa_type_index);
              }
              else
              {
                ++start_helix_counts( aa_type_index);
              }
            }
          }

          if( !m_ConsiderEnvironment)
          {
            continue;
          }

          // handle environment
          if( env_types( position_first) == env_types( position_sec))
          {
            bool is_start_end( false);
            biol::EnvironmentType nominal_type( env_types( position_sec));
            for( size_t p( position_first + 1); p < position_sec; ++p)
            {
              if( env_types( p) != nominal_type)
              {
                is_start_end = true;
                break;
              }
            }
            if( nominal_type == biol::GetEnvironmentTypes().e_MembraneCore)
            {
              if( is_start_end)
              {
                ++start_membrane_counts( aa_type_index);
                ++end_membrane_counts( aa_type_index);
              }
              else
              {
                // two helices
                ++membrane_membrane_counts( aa_type_index);
              }
            }
            else if( nominal_type == biol::GetEnvironmentTypes().e_Transition)
            {
              if( is_start_end)
              {
                ++start_trans_counts( aa_type_index);
                ++end_trans_counts( aa_type_index);
              }
              else
              {
                // two transs
                ++trans_trans_counts( aa_type_index);
              }
            }
            else if( nominal_type == biol::GetEnvironmentTypes().e_Solution)
            {
              if( is_start_end)
              {
                ++start_sol_counts( aa_type_index);
                ++end_sol_counts( aa_type_index);
              }
              else
              {
                // two sols
                ++sol_sol_counts( aa_type_index);
              }
            }
          }
          else if( env_types( position_first) == biol::GetEnvironmentTypes().e_Transition)
          {
            ++end_trans_counts( aa_type_index);
            if( env_types( position_sec) == biol::GetEnvironmentTypes().e_MembraneCore)
            {
              ++start_membrane_counts( aa_type_index);
            }
            else
            {
              ++start_sol_counts( aa_type_index);
            }
          }
          else if( env_types( position_first) == biol::GetEnvironmentTypes().e_MembraneCore)
          {
            ++end_membrane_counts( aa_type_index);
            if( env_types( position_sec) == biol::GetEnvironmentTypes().e_Transition)
            {
              ++start_trans_counts( aa_type_index);
            }
            else
            {
              ++start_sol_counts( aa_type_index);
            }
          }
          else
          {
            // first is sol
            ++end_sol_counts( aa_type_index);
            if( env_types( position_sec) == biol::GetEnvironmentTypes().e_Transition)
            {
              ++start_trans_counts( aa_type_index);
            }
            else
            {
              ++start_membrane_counts( aa_type_index);
            }
          }
        }
      }
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int GenerateAAPairStatistics::Main() const
    {
      io::OFStream output;

      m_StrandStrandCounts.Reset();
      m_HelixHelixCounts.Reset();
      m_CoilCoilCounts.Reset();
      m_StartsHelixCounts.Reset();
      m_StartsStrandCounts.Reset();
      m_StartsCoilCounts.Reset();
      m_EndsHelixCounts.Reset();
      m_EndsStrandCounts.Reset();
      m_EndsCoilCounts.Reset();
      m_AAPairTypeCounts.Reset();
      m_Separation.Reset();
      m_TransitionTransitionCounts.Reset();
      m_MembraneMembraneCounts.Reset();
      m_SolubleSolubleCounts.Reset();
      m_StartsMembraneCounts.Reset();
      m_StartsTransitionCounts.Reset();
      m_StartsSolubleCounts.Reset();
      m_EndsMembraneCounts.Reset();
      m_EndsTransitionCounts.Reset();
      m_EndsSolubleCounts.Reset();

      m_ConsiderEnvironment = m_ConsiderEnvironmentFlag->GetFlag();

      if( m_RestrictEnvironmentFlag->GetFlag())
      {
        m_EnvironmentType = biol::EnvironmentType( m_RestrictEnvironmentFlag->GetFirstParameter()->GetValue());
      }

      // get the maximum distance
      m_MaxDistance = m_MaxDistanceFlag->GetFirstParameter()->GetNumericalValue< size_t>();

      // set the method
      m_Method = sspred::Method( m_SSPredFlag->GetFirstParameter()->GetValue());

      // number of potential AA type pairs
      const size_t number_pairs( math::Sqr( size_t( biol::AATypes::s_NumberStandardAATypes) + 1));
      // initialize all counts and statistics
      m_StrandStrandCounts.Resize( m_MaxDistance, linal::Vector< size_t>( number_pairs, size_t( 0)));
      m_HelixHelixCounts.Resize( m_MaxDistance, linal::Vector< size_t>( number_pairs, size_t( 0)));
      m_CoilCoilCounts.Resize( m_MaxDistance, linal::Vector< size_t>( number_pairs, size_t( 0)));
      m_StartsHelixCounts.Resize( m_MaxDistance, linal::Vector< size_t>( number_pairs, size_t( 0)));
      m_StartsStrandCounts.Resize( m_MaxDistance, linal::Vector< size_t>( number_pairs, size_t( 0)));
      m_StartsCoilCounts.Resize( m_MaxDistance, linal::Vector< size_t>( number_pairs, size_t( 0)));
      m_EndsHelixCounts.Resize( m_MaxDistance, linal::Vector< size_t>( number_pairs, size_t( 0)));
      m_EndsStrandCounts.Resize( m_MaxDistance, linal::Vector< size_t>( number_pairs, size_t( 0)));
      m_EndsCoilCounts.Resize( m_MaxDistance, linal::Vector< size_t>( number_pairs, size_t( 0)));
      m_AAPairTypeCounts.Resize( m_MaxDistance, linal::Vector< size_t>( number_pairs, size_t( 0)));
      m_Separation.Resize( m_MaxDistance, storage::Vector< math::RunningAverage< double> >( number_pairs));
      m_TransitionTransitionCounts.Resize( m_MaxDistance, linal::Vector< size_t>( number_pairs, size_t( 0)));
      m_MembraneMembraneCounts.Resize( m_MaxDistance, linal::Vector< size_t>( number_pairs, size_t( 0)));
      m_SolubleSolubleCounts.Resize( m_MaxDistance, linal::Vector< size_t>( number_pairs, size_t( 0)));
      m_StartsMembraneCounts.Resize( m_MaxDistance, linal::Vector< size_t>( number_pairs, size_t( 0)));
      m_StartsTransitionCounts.Resize( m_MaxDistance, linal::Vector< size_t>( number_pairs, size_t( 0)));
      m_StartsSolubleCounts.Resize( m_MaxDistance, linal::Vector< size_t>( number_pairs, size_t( 0)));
      m_EndsMembraneCounts.Resize( m_MaxDistance, linal::Vector< size_t>( number_pairs, size_t( 0)));
      m_EndsTransitionCounts.Resize( m_MaxDistance, linal::Vector< size_t>( number_pairs, size_t( 0)));
      m_EndsSolubleCounts.Resize( m_MaxDistance, linal::Vector< size_t>( number_pairs, size_t( 0)));

      for
      (
        util::ShPtrVector< command::ParameterInterface>::const_iterator
          itr_source( m_InputFlag->GetParameterList().Begin()), itr_source_end( m_InputFlag->GetParameterList().End());
        itr_source != itr_source_end;
        ++itr_source
      )
      {
        assemble::ProteinWithCacheStorageFile protein_reader( true, true);
        BCL_Assert
        (
          protein_reader.TryRead( ( *itr_source)->GetValue(), util::GetLogger()),
          "Could not initialize storage"
        );

        // get all keys from the reader
        storage::Vector< std::string> protein_keys( protein_reader.GetAllKeys());

        // iterate over the fragments from whatever feeds were given
        size_t model_number( 0);
        float model_fraction( 1.0 / float( protein_keys.GetSize()));
        for
        (
          storage::Vector< std::string>::const_iterator itr_key( protein_keys.Begin()), itr_key_end( protein_keys.End());
          itr_key != itr_key_end;
          ++itr_key, ++model_number
        )
        {
          // log status
          util::GetLogger().LogStatus
          (
            "Protein model "
            + util::Format()( model_number) + " / " + util::Format()( protein_keys.GetSize())
            + " " + util::Format().FFP( 2)( model_number * model_fraction) + " " + *itr_key
          );

          // read the key
          util::ShPtr< assemble::ProteinModelWithCache> protein_model_ptr( protein_reader.Retrieve( *itr_key));

          // process the protein
          Process( *protein_model_ptr);
        }
      }

      const std::string file_prefix( m_OutputFilenameBase->GetFirstParameter()->GetValue());

      std::string one_letter_codes;
      const size_t number_standard_aa_types( biol::AATypes::s_NumberStandardAATypes);
      for( size_t aa_type( 0); aa_type < number_standard_aa_types; ++aa_type)
      {
        one_letter_codes += biol::AAType( aa_type)->GetOneLetterCode();
      }
      one_letter_codes += 'X';
      const size_t number_aa_types( biol::AATypes::s_NumberStandardAATypes + 1);

      util::SiPtrVector< storage::Vector< linal::Vector< size_t> > > counts_vectors;
      counts_vectors.PushBack( m_AAPairTypeCounts);
      counts_vectors.PushBack( m_StrandStrandCounts);
      counts_vectors.PushBack( m_HelixHelixCounts);
      counts_vectors.PushBack( m_CoilCoilCounts);
      counts_vectors.PushBack( m_StartsHelixCounts);
      counts_vectors.PushBack( m_StartsStrandCounts);
      counts_vectors.PushBack( m_StartsCoilCounts);
      counts_vectors.PushBack( m_EndsHelixCounts);
      counts_vectors.PushBack( m_EndsStrandCounts);
      counts_vectors.PushBack( m_EndsCoilCounts);
      counts_vectors.PushBack( m_TransitionTransitionCounts);
      counts_vectors.PushBack( m_MembraneMembraneCounts);
      counts_vectors.PushBack( m_SolubleSolubleCounts);
      counts_vectors.PushBack( m_StartsMembraneCounts);
      counts_vectors.PushBack( m_StartsTransitionCounts);
      counts_vectors.PushBack( m_StartsSolubleCounts);
      counts_vectors.PushBack( m_EndsMembraneCounts);
      counts_vectors.PushBack( m_EndsTransitionCounts);
      counts_vectors.PushBack( m_EndsSolubleCounts);

      // compute averages for each type
      for( size_t aa_pair_type( 0); aa_pair_type < number_pairs; ++aa_pair_type)
      {
        const size_t first_aa_type_id( aa_pair_type / number_aa_types);
        const size_t second_aa_type_id( aa_pair_type % number_aa_types);
        if( first_aa_type_id == number_standard_aa_types || second_aa_type_id == number_standard_aa_types)
        {
          // skip unknown types
          continue;
        }
        const size_t first_x_id( first_aa_type_id * number_aa_types + number_standard_aa_types);
        const size_t second_x_id( number_standard_aa_types * number_aa_types + second_aa_type_id);
        for( size_t distance( 0); distance < m_MaxDistance; ++distance)
        {
          for
          (
            util::SiPtrVector< storage::Vector< linal::Vector< size_t> > >::iterator
              itr_counts( counts_vectors.Begin()), itr_counts_end( counts_vectors.End());
            itr_counts != itr_counts_end;
            ++itr_counts
          )
          {
            const size_t count( ( **itr_counts)( distance)( aa_pair_type));
            ( **itr_counts)( distance)( first_x_id) += count;
            ( **itr_counts)( distance)( second_x_id) += count;
          }
          m_Separation( distance)( first_x_id) += m_Separation( distance)( aa_pair_type);
          m_Separation( distance)( second_x_id) += m_Separation( distance)( aa_pair_type);
        }
      }
      const size_t last_x_row( number_standard_aa_types * number_aa_types);
      const size_t last_x_id( last_x_row + number_standard_aa_types);
      for( size_t distance( 0); distance < m_MaxDistance; ++distance)
      {
        for
        (
          util::SiPtrVector< storage::Vector< linal::Vector< size_t> > >::iterator
            itr_counts( counts_vectors.Begin()), itr_counts_end( counts_vectors.End());
          itr_counts != itr_counts_end;
          ++itr_counts
        )
        {
          ( **itr_counts)( distance)( last_x_id)
            = linal::VectorConstReference< size_t>
              (
                number_standard_aa_types,
                ( **itr_counts)( distance).Begin() + last_x_row
              ).Sum();
        }
        for( size_t i( 0); i < number_standard_aa_types; ++i)
        {
          m_Separation( distance)( last_x_row + number_standard_aa_types) += m_Separation( distance)( last_x_row + i);
        }
      }

      io::File::MustOpenOFStream( output, file_prefix + ".txt");
      for( size_t aa_pair_type( 0); aa_pair_type < number_pairs; ++aa_pair_type)
      {
        for( size_t distance( 0); distance < m_MaxDistance; ++distance)
        {
          const float count( std::max( m_AAPairTypeCounts( distance)( aa_pair_type), size_t( 1)));
          const float dist_count( count * ( distance + 1));
          output
            << one_letter_codes[ aa_pair_type / number_aa_types]
            << ' '
            << one_letter_codes[ aa_pair_type % number_aa_types]
            << ' '
            << distance + 1
            << ' '
            << m_AAPairTypeCounts( distance)( aa_pair_type)
            << ' '
            << m_HelixHelixCounts( distance)( aa_pair_type) / count
            << ' '
            << m_StartsHelixCounts( distance)( aa_pair_type) / dist_count
            << ' '
            << m_EndsHelixCounts( distance)( aa_pair_type) / dist_count
            << ' '
            << m_StrandStrandCounts( distance)( aa_pair_type) / count
            << ' '
            << m_StartsStrandCounts( distance)( aa_pair_type) / dist_count
            << ' '
            << m_EndsStrandCounts( distance)( aa_pair_type) / dist_count
            << ' '
            << m_CoilCoilCounts( distance)( aa_pair_type) / count
            << ' '
            << m_StartsCoilCounts( distance)( aa_pair_type) / dist_count
            << ' '
            << m_EndsCoilCounts( distance)( aa_pair_type) / dist_count
            << ' '
            << m_Separation( distance)( aa_pair_type).GetAverage()
            << '\n';
        }
      }
      io::File::CloseClearFStream( output);

      if( !m_ConsiderEnvironment)
      {
        return 0;
      }

      io::File::MustOpenOFStream( output, file_prefix + ".env.txt");
      for( size_t aa_pair_type( 0); aa_pair_type < number_pairs; ++aa_pair_type)
      {
        for( size_t distance( 0); distance < m_MaxDistance; ++distance)
        {
          const float count( std::max( m_AAPairTypeCounts( distance)( aa_pair_type), size_t( 1)));
          const float dist_count( count * ( distance + 1));
          output
            << one_letter_codes[ aa_pair_type / number_aa_types]
            << ' '
            << one_letter_codes[ aa_pair_type % number_aa_types]
            << ' '
            << distance + 1
            << ' '
            << m_AAPairTypeCounts( distance)( aa_pair_type)
            << ' '
            << m_MembraneMembraneCounts( distance)( aa_pair_type) / count
            << ' '
            << m_StartsMembraneCounts( distance)( aa_pair_type) / dist_count
            << ' '
            << m_EndsMembraneCounts( distance)( aa_pair_type) / dist_count
            << ' '
            << m_TransitionTransitionCounts( distance)( aa_pair_type) / count
            << ' '
            << m_StartsTransitionCounts( distance)( aa_pair_type) / dist_count
            << ' '
            << m_EndsTransitionCounts( distance)( aa_pair_type) / dist_count
            << ' '
            << m_SolubleSolubleCounts( distance)( aa_pair_type) / count
            << ' '
            << m_StartsSolubleCounts( distance)( aa_pair_type) / dist_count
            << ' '
            << m_EndsSolubleCounts( distance)( aa_pair_type) / dist_count
            << '\n';
        }
      }
      io::File::CloseClearFStream( output);
      // end
      return 0;
    }

  } // namespace app
} // namespace bcl

