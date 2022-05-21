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
#include "fold/bcl_fold_setup.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_protein_storage_file.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "command/bcl_command_flag_separator.h"
#include "fold/bcl_fold_default_flags.h"
#include "fold/bcl_fold_protocol_ensemble.h"
#include "fold/bcl_fold_protocols.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "sspred/bcl_sspred_method_handler.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

    //! @brief flag for returning the fold setup used
    //! @return the fold setup used
    const Setup &GetSetup()
    {
      return Setup::GetStaticInstance();
    }

  //////////
  // data //
  //////////

    //! @brief construct and return a static instance of the class
    //! @return a static instance of the class
    const Setup &Setup::GetStaticInstance()
    {
      return GetStaticInstanceNonConst();
    }

    //! @brief reset the static instance
    void Setup::InitializeStaticInstance()
    {
      GetStaticInstanceNonConst() = Setup();
    }

    //! @brief return vector of all flags applicable
    //! @return vector of all flags applicable
    util::ShPtrVector< command::FlagInterface> Setup::GetAllFlags()
    {
      // initialize static data to store flags
      util::ShPtrVector< command::FlagInterface> flags;

      // iterate over all the protocols possible
      for
      (
        std::vector< Protocol>::const_iterator
          protocol_itr( GetProtocols().Begin()), protocol_itr_end( GetProtocols().End());
        protocol_itr != protocol_itr_end; ++protocol_itr
      )
      {
        // insert all the flags in the vector
        const util::ShPtr< command::FlagInterface> sp_separator
        (
          new command::FlagSeparator
          (
            "\n==============================================================================================\n" +
            protocol_itr->GetName() + ":\t" + ( **protocol_itr)->GetDescription() + "\n\n"
          )
        );
        flags.PushBack( sp_separator);

        // pushback the flags for this protocol
        flags.Append( ( **protocol_itr)->GetAllFlags());
      }

      // end
      return flags;
    }

    //! @brief return string containing readme information
    //! @return string containing readme information
    const std::string &Setup::GetReadme()
    {
      // initialize readme string
      static std::string s_readme;

      // if readme is empty
      if( s_readme.empty())
      {
        // iterate over all the protocols possible
        for
        (
          std::vector< Protocol>::const_iterator
            protocol_itr( GetProtocols().Begin()), protocol_itr_end( GetProtocols().End());
          protocol_itr != protocol_itr_end; ++protocol_itr
        )
        {
          s_readme +=
            "===============================================\n" +
            protocol_itr->GetName() + "\n" +
            "===============================================\n" +
            ( **protocol_itr)->GetReadMe() + "\n\n";
        }
      }

      // end
      return s_readme;
    }

    //! @brief construct and return a non-const static instance of the class
    //! @return a non-const static instance of the class
    Setup &Setup::GetStaticInstanceNonConst()
    {
      // initialize static instance
      static Setup s_static_instance;

      // end
      return s_static_instance;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief private default constructor
    Setup::Setup() :
      m_EmptyModel(),
      m_NativeModel(),
      m_StartModel(),
      m_NumberRounds( DefaultFlags::GetFlagNumberModels()->GetFirstParameter()->GetNumericalValue< size_t>()),
      m_QualityMeasures( quality::Measures::GetCommandLineQualityMeasures()),
      m_SuperimposeMeasure
      (
        quality::SuperimposeMeasure
        (
          quality::SuperimposeMeasures::GetFlagSuperimposeMeasure()->GetFirstParameter()->GetValue()
        )
      ),
      m_Prefix( DefaultFlags::GetFlagPrefix()->GetFirstParameter()->GetValue()),
      m_StepStatuses(),
      m_Storage( assemble::ProteinStorageFile::GetDefaultStorage())
    {
      // if no native was given but quality measures were specified
      if( !DefaultFlags::GetFlagNativeModel()->GetFirstParameter()->GetWasSetInCommandLine())
      {
        // if quality measures were specified
        BCL_Assert
        (
          m_QualityMeasures.IsEmpty(),
          "No native model was given, but quality measures were specified!"
        );

        // if superimpose measure was specified
        BCL_Assert
        (
          m_SuperimposeMeasure == quality::GetSuperimposeMeasures().e_NoSuperimpose,
          "No native model was given, but a superimpose measure was specified!"
        );
      }

      // iterate over given step statuses
      for
      (
        util::ShPtrVector< command::ParameterInterface>::const_iterator
          step_status_itr( DefaultFlags::GetFlagPrintMinimization()->GetParameterList().Begin()),
          step_status_itr_end( DefaultFlags::GetFlagPrintMinimization()->GetParameterList().End());
        step_status_itr != step_status_itr_end;
        ++step_status_itr
      )
      {
        // insert into the set
        m_StepStatuses.Insert( opti::StepStatusEnum( ( *step_status_itr)->GetValue()));
      }
    }

    //! @brief Clone function
    //! @return pointer to new Setup
    Setup *Setup::Clone() const
    {
      return new Setup( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &Setup::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the native/template model if one was provided
    //! @return the native/template model if one was provided
    const util::ShPtr< assemble::ProteinModel> &Setup::GetNativeModel() const
    {
      // if the native model is not initialized and a pdb filename was given at the commandline
      if( !m_NativeModel.IsDefined() && DefaultFlags::GetFlagNativeModel()->GetFlag())
      {
        // read the model
        m_NativeModel =
          util::ShPtr< assemble::ProteinModel>
          (
            pdb::Factory().ProteinModelFromPDBFilename
            (
              DefaultFlags::GetFlagNativeModel()->GetFirstParameter()->GetValue()
            ).Clone()
          );

        BCL_Assert( m_NativeModel.IsDefined(), "native model is not defined");
      }

      // return native model
      return m_NativeModel;
    }

    //! @brief return the protein model constructed from aa sequences (given as fastas) that are supposed to be folded
    //! @return the protein model constructed from aa sequences (given as fastas) that are supposed to be folded
    const util::ShPtr< assemble::ProteinModel> &Setup::GetEmptyModel() const
    {
      // empty model is not initialized
      if( !m_EmptyModel.IsDefined())
      {
        // make sure the flags are correctly set
        BCL_Assert
        (
          DefaultFlags::GetFlagNativeModel()->GetFlag() ||
          DefaultFlags::GetFlagFastaRead()->GetFlag() ||
          DefaultFlags::GetFlagStartModel()->GetFlag(),
          "One of the flags for native, fasta, or start model should be set to construct fold setup!!"
        );

        // if native pool flag was given
        if( DefaultFlags::GetFlagUseNativeSSEsAsPool()->GetFlag())
        {
          BCL_Assert
          (
            DefaultFlags::GetFlagNativeModel()->GetFlag(),
            "A native structure has to be provided for using native pool flag!!"
          );
        }

        // construct the ProteinModelData
        util::ShPtr< assemble::ProteinModelData> sp_model_data( new assemble::ProteinModelData());

        // set the identification
        std::string identification( m_Prefix);

        // if a native was given
        if( DefaultFlags::GetFlagNativeModel()->GetFlag())
        {
          // set the identification to the filename minus path and extension
          identification = io::File::RemoveFullExtension
          (
            io::File::SplitToPathAndFileName
            (
              DefaultFlags::GetFlagNativeModel()->GetFirstParameter()->GetValue()
            ).Second()
          );
        }

        // update the protein model data w/ the identification
        sp_model_data->Insert
        (
          assemble::ProteinModelData::e_Identification,
          util::ShPtr< util::Wrapper< std::string> >( new util::Wrapper< std::string>( identification))
        );

        // initialize map to hold the min pool lengths
        storage::Map< biol::SSType, size_t> min_pool_sse_lengths( assemble::SSEPool::GetCommandLineMinSSELengths());

        // initialize empty pool
        util::ShPtr< assemble::SSEPool> sp_pool( new assemble::SSEPool());

        // if a native pdb was provided
        if( DefaultFlags::GetFlagNativeModel()->GetFirstParameter()->GetWasSetInCommandLine())
        {
          // if fasta flag was given give a warning
          if( DefaultFlags::GetFlagFastaRead()->GetFlag())
          {
            BCL_MessageCrt( "Given fasta flag is skipped since a native structure is provided");
          }

          // update empty model to native model
          m_EmptyModel =
            util::ShPtr< assemble::ProteinModel>( new assemble::ProteinModel( GetNativeModel()->GetEmptyChains()));

          // insert the native model into model data
          sp_model_data->Insert( assemble::ProteinModelData::e_NativeModel, GetNativeModel());

          // make a copy of the native model
          util::ShPtr< assemble::ProteinModel> sp_native_model_copy( GetNativeModel()->Clone());

          // filter the SSEs by given pool min sse sizes
          sp_native_model_copy->FilterByMinSSESizes( min_pool_sse_lengths);

          // insert the filtered native model to the PM data
          sp_model_data->Insert( assemble::ProteinModelData::e_NativeFilteredModel, sp_native_model_copy);

          // if the flag for using native SSE definitions were given
          if( DefaultFlags::GetFlagUseNativeSSEsAsPool()->GetFlag())
          {
            // if fasta flag was given give a warning
            if( assemble::SSEPool::GetFlagPoolRead()->GetFlag())
            {
              BCL_MessageCrt( "Given pool file is skipped because use native pool flag is given!");
            }

            // update the pool to have the corresponding SSE definitions from the native model
            sp_pool = util::ShPtr< assemble::SSEPool>
            (
              new assemble::SSEPool
              (
                sp_native_model_copy->GetSSEs(),
                true,
                DefaultFlags::GetFlagUseNativeSSEsAsPool()->GetFirstParameter()->GetValue() == "ideal"
              )
            );
          }
          else
          {
            // open pool file
            io::IFStream read;
            const std::string pool_file( assemble::SSEPool::GetFlagPoolRead()->GetFirstParameter()->GetValue());
            BCL_MessageCrt( "Reading pool from file " + pool_file);
            io::File::MustOpenIFStream( read, pool_file);

            // set the identification to the filename minus path and extension
            identification = io::File::RemoveFullExtension
            (
              io::File::SplitToPathAndFileName
              (
                DefaultFlags::GetFlagNativeModel()->GetFirstParameter()->GetValue()
              ).Second()
            );

            // read pool
            sp_pool->ReadSSEPool
            (
              read,
              *m_EmptyModel,
              min_pool_sse_lengths[ biol::GetSSTypes().HELIX],
              min_pool_sse_lengths[ biol::GetSSTypes().STRAND]
            );
          }
        }
        // if no native pdb is given, but just fasta files
        else
        {
          // initialize empty model
          m_EmptyModel = util::ShPtr< assemble::ProteinModel>( new assemble::ProteinModel());

          // create set for all chain ids represented in the pool
          const storage::Set< char> pool_chain_ids
          (
            assemble::SSEPool::GetChainsRepresented( assemble::SSEPool::GetFlagPoolRead()->GetFirstParameter()->GetValue())
          );

          // fasta chain ids stored
          storage::Set< char> fasta_chain_ids;

          // create a itr for looping over chain ids
          util::ShPtrVector< command::ParameterInterface>::const_iterator chain_id_itr;

          // if chain id flag was given
          if( DefaultFlags::GetFlagChainIdRead()->GetFlag())
          {
            // make sure the number of fastas and chain ids is the same
            BCL_Assert
            (
              DefaultFlags::GetFlagFastaRead()->GetParameterList().GetSize()
                == DefaultFlags::GetFlagChainIdRead()->GetParameterList().GetSize(),
              "Same number of fasta files and chain ids required!"
            );

            // initialize itr for the following loop
            chain_id_itr = DefaultFlags::GetFlagChainIdRead()->GetParameterList().Begin();
          }

          // iterate through the list of fastas and chain ids
          for
          (
            util::ShPtrVector< command::ParameterInterface>::const_iterator
              fasta_itr( DefaultFlags::GetFlagFastaRead()->GetParameterList().Begin()),
              fasta_itr_end( DefaultFlags::GetFlagFastaRead()->GetParameterList().End());
            fasta_itr != fasta_itr_end;
            ++fasta_itr
          )
          {
            // remove fast extension
            const std::string fasta_id( io::File::RemoveFullExtension( ( *fasta_itr)->GetValue()));

            // get the chain id, either the last char before the ".fasta" extension or in the list of chain ids
            const char chain_id
            (
              DefaultFlags::GetFlagChainIdRead()->GetFlag() ?
                ( *chain_id_itr)->GetValue()[ 0] :
                fasta_id[ fasta_id.length() - 1]
            );

            // make sure this chain is also found in the pool chain ids
            BCL_Assert
            (
              pool_chain_ids.Find( chain_id) != pool_chain_ids.End(),
              "The given chain id \"" + util::Format()( chain_id) + "\" is not found in the pool"
            );

            // check if chain id was already inserted
            BCL_Assert
            (
              fasta_chain_ids.Find( chain_id) == fasta_chain_ids.End(),
              "fasta files with identical chain identifiers have been found"
            );

            // insert chain id
            fasta_chain_ids.Insert( chain_id);

            // open fasta file
            io::IFStream read;

            // set the identification to the filename minus path and extension
            identification = io::File::RemoveFullExtension
            (
              io::File::SplitToPathAndFileName
              (
                DefaultFlags::GetFlagNativeModel()->GetFirstParameter()->GetValue()
              ).Second()
            );
               io::File::MustOpenIFStream( read, ( *fasta_itr)->GetValue());

            // insert chain for fasta file
            m_EmptyModel->Insert
            (
              pdb::Factory( biol::GetAAClasses().e_AABackBone).ChainFromFastaStream( chain_id, read)
            );

            // if chain ids were given on the command line increase the chain id itr
            if( DefaultFlags::GetFlagChainIdRead()->GetFlag())
            {
              ++chain_id_itr;
            }
          }

          // open pool file
          io::IFStream read;
          io::File::MustOpenIFStream( read, assemble::SSEPool::GetFlagPoolRead()->GetFirstParameter()->GetValue());

          // read pool
          sp_pool->ReadSSEPool
          (
            read,
            *m_EmptyModel,
            min_pool_sse_lengths[ biol::GetSSTypes().HELIX],
            min_pool_sse_lengths[ biol::GetSSTypes().STRAND]
          );

          // close stream
          io::File::CloseClearFStream( read);
        }

        // make sure pool has SSEs
        BCL_Assert( !sp_pool->IsEmpty(), "Pool must contain SSEs for assembly!");

        // if separate pool flag was provided and native pool is being used
        if( DefaultFlags::GetFlagPoolSeparate()->GetFlag() && !sp_pool->IsOverlapping())
        {
          // separate pools
          const bool success
          (
            sp_pool->Separate
            (
              assemble::SSEPool::GetCommandLineMinSSELengths(),
              DefaultFlags::GetFlagPoolSeparate()->GetFirstParameter()->GetNumericalValue< size_t>()
            )
          );
          if( !success)
          {
            BCL_MessageCrt( "Separate pools operation unsuccessful")
          }

          // make sure to prune it to remove loops
          sp_pool->Prune( assemble::SSEPool::GetCommandLineMinSSELengths());

        }

        // if the list is not empty
        storage::Set< sspred::Method> ss_pred_methods( sspred::Methods::GetCommandLineMethods());

        if( !ss_pred_methods.IsEmpty())
        {

          std::string prefix( DefaultFlags::GetFlagReadSequenceDataPath()->GetParameterList()( 1)->GetValue());
          std::string path( DefaultFlags::GetFlagReadSequenceDataPath()->GetFirstParameter()->GetValue());

          //add ss predictions to m_StartModel
          sspred::MethodHandler::ReadPredictionsForProteinModel
          (
            ss_pred_methods,
            *m_EmptyModel,
            DefaultFlags::GetFlagReadSequenceDataPath()->GetParameterList()( 1)->GetValue(),
            DefaultFlags::GetFlagReadSequenceDataPath()->GetFirstParameter()->GetValue()
          );
        }

        sp_model_data->Insert( assemble::ProteinModelData::e_Pool, sp_pool);
        // set data on the model
        m_EmptyModel->SetProteinModelData( sp_model_data);

        // true if ensemble of models is going to be folded
        if( ProtocolEnsemble::GetFlagEnsembleSize()->GetFlag())
        {
          BCL_MessageCrt( "setting ensemble");
          util::ShPtr< assemble::ProteinEnsemble> ensemble
          (
            new assemble::ProteinEnsemble( util::ShPtr< assemble::ProteinModel>( m_EmptyModel->HardCopy()))
          );

          m_EmptyModel = ensemble;
        }
      }

      // return
      return m_EmptyModel;
    }

    //! @brief return the starting model after adding sses from -start_model file to it
    //!        (this is used to refine an existing starting model)
    //! @return return the starting model after adding sses from -start_model file to it
    const util::ShPtr< assemble::ProteinModel> Setup::GetStartModel() const
    {
      // start model is not initialized
      if( !m_StartModel.IsDefined())
      {
        // use the empty model to get the sequences
        m_StartModel = util::ShPtr< assemble::ProteinModel>( GetEmptyModel()->Clone());

        // if a pdb was provided
        if( DefaultFlags::GetFlagStartModel()->GetFlag())
        {
          util::ShPtr< assemble::ProteinModel> sp_model_provided
          (
            new assemble::ProteinModel
            (
              pdb::Factory().ProteinModelFromPDBFilename( DefaultFlags::GetFlagStartModel()->GetFirstParameter()->GetValue())
            )
          );

          // iterate over the chains in the start model
          for
          (
            util::ShPtrVector< assemble::Chain>::iterator
              chain_itr( sp_model_provided->GetChains().Begin()),
              chain_itr_end( sp_model_provided->GetChains().End());
            chain_itr != chain_itr_end; ++chain_itr
          )
          {
            // iterate over each SSE in that given chain
            for
            (
              storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
                sse_itr( ( *chain_itr)->GetData().Begin()), sse_itr_end( ( *chain_itr)->GetData().End());
              sse_itr != sse_itr_end; ++sse_itr
            )

            // if this SSE has defined coordinates
            if( ( *sse_itr)->HasDefinedCoordinates())
            {
              // insert the SSEs from the chain
              m_StartModel->Insert( *sse_itr);
            }
          }

          // connect the sequences of the SSEs from the start model to the native model
          m_StartModel->ConnectSSEToChainData();

          // update any protein model data
          util::ShPtr< assemble::ProteinModelData> sp_data( m_StartModel->GetProteinModelData());
          sp_data->Insert( *( sp_model_provided->GetProteinModelData()));
          m_StartModel->SetProteinModelData( sp_data);
        }
      }

      // return model
      return m_StartModel;
    }

    //! @brief return quality measures to be calculated
    //! @return quality measures to be calculated
    const storage::Set< quality::Measure> &Setup::GetQualityMeasures() const
    {
      return m_QualityMeasures;
    }

    //! @brief return superimpose measure
    //! @return quality superimpose measure
    const quality::SuperimposeMeasure &Setup::GetSuperimposeMeasure() const
    {
      return m_SuperimposeMeasure;
    }

    //! @brief gives the prefix object
    //! @return the prefix that is prepended to output files
    const std::string &Setup::GetPrefix() const
    {
      return m_Prefix;
    }

    //! @brief gives the set of step statuses that will be printed by printers
    //! @return set of step statuses that will be printed by printers
    const storage::Set< opti::StepStatusEnum> &Setup::GetStepStatuses() const
    {
      return m_StepStatuses;
    }

    //! @brief get the protein storage
    //! @return protein storage
    const util::ShPtr< assemble::ProteinStorageFile> &Setup::GetStorage() const
    {
      return m_Storage;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Setup::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Setup::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
