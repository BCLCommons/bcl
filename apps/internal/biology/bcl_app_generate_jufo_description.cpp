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
#include "app/bcl_app_apps.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa.h"
#include "biol/bcl_biol_blast_profile_handler.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "descriptor/bcl_descriptor_base.h"
#include "descriptor/bcl_descriptor_iterator.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "sspred/bcl_sspred_jufo9d.h"
#include "storage/bcl_storage_template_instantiations.h"
#include "util/bcl_util_implementation.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class GenerateJufoDescription
    //! @brief This application is used for generating ANN training data for JUFO9D
    //!
    //! @author karakam
    //! @date 12/10/2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API GenerateJufoDescription :
      public Interface
    {
    private:

    //////////
    // data //
    //////////

      //! path for input files
      util::ShPtr< command::FlagStatic> m_PathFlag;
      util::ShPtr< command::ParameterInterface> m_PathParam;

      //! flag for reading predictions from different SSMethods
      util::ShPtr< command::FlagStatic> m_PdbListFlag;
      util::ShPtr< command::ParameterInterface> m_PdbListParam;

      //! flag for using pdb path hierarchy for all input files
      util::ShPtr< command::FlagInterface> m_PdbPathHierarchy;

      //! flag to define the output filename
      util::ShPtr< command::FlagInterface> m_OutputPrefixFlag;

      //! flag for providing a oligomeric state dictionary
      util::ShPtr< command::FlagStatic> m_OligomericStateDictionaryFlag;
      util::ShPtr< command::ParameterInterface> m_OligomericStateDictionaryParam;

      //! flag to define the number of requested entries for each of the states
      util::ShPtr< command::FlagStatic> m_NumberEntriesPerStateFlag;
      util::ShPtr< command::ParameterInterface> m_NumberEntriesPerStateParam;

      //! flag to create unbalanced input
      util::ShPtr< command::FlagInterface> m_UnbalancedInput;

      //! flag to create unbalanced input
      util::ShPtr< command::FlagInterface> m_NonrandomInput;

      //! flag to create input for testing, includes all residues in sequence, none skipped
      util::ShPtr< command::FlagInterface> m_TestingInput;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! @brief default constructor
      GenerateJufoDescription();

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      GenerateJufoDescription *Clone() const
      {
        return new GenerateJufoDescription( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief return first layer descriptors for the given amino acid sequence
      //! @param SEQUENCE AASequence of interest
      //! @param OLIGOMERIC_STATE oligomeric state of this sequence
      //! @return map that stores the descriptions for each
      storage::Map< util::SiPtr< const biol::AABase>, linal::Vector< double> >
        CalculateFirstLayerDescriptors( util::ShPtr< biol::AASequence> &SEQUENCE, const size_t OLIGOMERIC_STATE) const;

      //! @brief initializes the command object for that executable
      //! @return initalized command object
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief Main function
      //! @return return value of the application
      int Main() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // end
        return OSTREAM;
      }

    private:

      //! static instance of this class
      static const ApplicationType GenerateJufoDescription_Instance;

    }; // class GenerateJufoDescription

    //! @brief Main function
    //! @return return value of the application
    int GenerateJufoDescription::Main() const
    {

    /////////////////////
    // initializations //
    /////////////////////

      const size_t number_states( 9);
      const size_t number_output_nodes( 9);

      // initialize the datasets
      storage::Vector< storage::Vector< storage::VectorND< 2, linal::Vector< double> > > >
        complete_dataset( number_states);
      storage::Vector< storage::VectorND< 2, linal::Vector< double> > > nonrandom_complete_dataset;

      // initialize pdb factory and handler
      pdb::Factory pdb_factory( biol::GetAAClasses().e_AABackBone);

      // initialize sse_min_sizes
      storage::Map< biol::SSType, size_t> sse_min_sizes;
      sse_min_sizes[ biol::GetSSTypes().HELIX] = 5;
      sse_min_sizes[ biol::GetSSTypes().STRAND] = 3;

      // initialize a pdb id list
      storage::List< storage::Pair< std::string, bool> > pdb_id_list;

      // read the pdb list
      io::IFStream read;
      io::File::MustOpenIFStream( read, m_PdbListParam->GetValue());

      // initialize temp storage for pdbid and the soluble state
      std::string pdb_id, soluble_state;

      BCL_MessageStd( "Reading pdb list");
      // while reading the file
      while( read >> pdb_id && read >> soluble_state && !read.eof())
      {
        // determine is_membrane;
        bool is_membrane;
        if( soluble_state == "Membrane")
        {
          is_membrane = true;
        }
        else if( soluble_state == "Soluble")
        {
          is_membrane = false;
        }
        else
        {
          BCL_MessageCrt( "The soluble state is undefined: " + soluble_state);
          return 0;
        }
        // form the pair and pushback
        pdb_id_list.PushBack( storage::Pair< std::string, bool>( pdb_id, is_membrane));
      }

      io::File::CloseClearFStream( read);

      // initialize the oligo dictionary
      storage::Map< std::string, size_t> oligo_dictionary;

      // if a dictionary file is given
      if( m_OligomericStateDictionaryFlag->GetFlag())
      {
        BCL_MessageStd( "Reading oligomeric state dictionary");

        io::IFStream oligo_read;
        // open the stream and start reading
        io::File::MustOpenIFStream( oligo_read, m_OligomericStateDictionaryParam->GetValue());
        while( !oligo_read.eof())
        {
          // read pdb tag and oligo state (0/1) value and put them into the map
          std::string tag;
          size_t oligo_state, number_subunits;
          oligo_read >> tag >> oligo_state >> number_subunits;
          oligo_dictionary[ tag] = oligo_state;
        }
        io::File::CloseClearFStream( oligo_read);
      }

      // construct the membrane
      const biol::Membrane membrane( biol::Membrane::GetCommandLineMembrane());

    ///////////////
    // iteration //
    ///////////////

      // initialize a pdb counter
      size_t pdb_ctr( 1);

      BCL_MessageStd( "Starting to read pdbs");

      // iterate over the pdb list
      for
      (
        storage::List< storage::Pair< std::string, bool> >::const_iterator
          pdb_itr( pdb_id_list.Begin()), pdb_itr_end( pdb_id_list.End());
        pdb_itr != pdb_itr_end; ++pdb_itr, ++pdb_ctr
      )
      {
        BCL_MessageCrt
        (
          "Processing pdb with tag: " + pdb_itr->First() + " " +
          util::Format()( pdb_ctr) + "/" + util::Format()( pdb_id_list.GetSize())
        );

        // initialize path
        std::string input_path( m_PathParam->GetValue() + PATH_SEPARATOR);
        // if pdb hierarchy is used
        if( m_PdbPathHierarchy->GetFlag())
        {
          input_path += pdb_itr->First().substr( 1, 2) + PATH_SEPARATOR;
        }

        // read pdb
        const std::string prefix_without_chainid( pdb_itr->First().substr( 0, 4));
        const char chain_id( pdb_itr->First()[ 4]);
        io::IFStream read;
        io::File::MustOpenIFStream( read, input_path + prefix_without_chainid + "biobcl.pdb");
        pdb::Handler pdb( read);
        io::File::CloseClearFStream( read);
        assemble::ProteinModel this_model( pdb_factory.ProteinModelFromPDB( pdb));

        // if testing input, use all residues for generating descriptors, even the unresolved ones
        if( m_TestingInput->GetFlag())
        {
          this_model.AddLoops( true, false);
        }

        // set the oligo state to 0
        size_t this_oligo_state( 0);

        // if the tag is in the map look up the oligostate for this pdb and store it
        if( oligo_dictionary.Find( prefix_without_chainid) != oligo_dictionary.End())
        {
          this_oligo_state = oligo_dictionary[ prefix_without_chainid];
        }

        // make a reference for the chain
        util::ShPtr< assemble::Chain> this_chain( this_model.GetChain( chain_id));

        BCL_Assert
        (
          this_chain.IsDefined(), "The chain with chainid " + util::Format()( chain_id) + " cannot be found!"
        );

        // make a reference of the sequence
        util::ShPtr< biol::AASequence> this_sequence( this_chain->GetSequence());

        // read the blast profile
        BCL_MessageVrb( "Reading blast profile");
        const std::string blast_filename( input_path + pdb_itr->First() + ".ascii6");

        // read blast profile
        io::File::MustOpenIFStream( read, blast_filename);
        biol::BlastProfileHandler::ReadProfileForAASequence( read, *this_sequence);
        io::File::CloseClearFStream( read);

        // instantiate an undefined amino acid with 0 blast profile and set JUFO2 to default
        biol::AA undefined_aa( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().XXX)));
        undefined_aa.SetBlastProfile( biol::BlastProfile( linal::Vector< double>( 20, 0.0)));

        BCL_MessageVrb( "Calculating descriptors");

        const storage::Map< util::SiPtr< const biol::AABase>, linal::Vector< double> > description_map
        (
          CalculateFirstLayerDescriptors( this_sequence, this_oligo_state)
        );

        // get all SSEs
        const storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap> all_sses( this_chain->GetData());

      /////////////////////////////
      // iteration over residues //
      /////////////////////////////

        // initialize iterators
        biol::AASequence::const_iterator seq_aa_itr( this_sequence->Begin());
        const biol::AASequence::const_iterator seq_aa_itr_end( this_sequence->End());

        BCL_MessageVrb( "iterating over SSEs");

        // iterate over all SSEs
        for
        (
          storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
            sse_itr( all_sses.Begin()), sse_itr_end( all_sses.End());
          sse_itr != sse_itr_end; ++sse_itr
        )
        {
          BCL_MessageVrb( "\tSSE: " + ( *sse_itr)->GetIdentification());
          // create iterator to the first aa
          biol::AASequence::const_iterator sse_aa_itr( ( *sse_itr)->Begin());
          const biol::AASequence::const_iterator sse_aa_itr_end( ( *sse_itr)->End());

          // while the seq_aa_itr and sse_aa_itr do not match
          // this can happen if a residue is not resolved
          while
          (
            ( *seq_aa_itr)->GetSeqID() < ( *sse_aa_itr)->GetSeqID() &&
            seq_aa_itr != seq_aa_itr_end && sse_aa_itr != sse_aa_itr_end
          )
          {
            // increment seq_aa_itr
            ++seq_aa_itr;
          }

          // store sstype
          biol::SSType this_sstype( ( *sse_itr)->GetType());

          // if the size of this SSE is less than the requested size
          if( !m_TestingInput->GetFlag() && ( *sse_itr)->GetSize() < sse_min_sizes[ this_sstype])
          {
            BCL_MessageDbg( "skipping SSE " + ( *sse_itr)->GetIdentification());
            continue;
          }

          BCL_MessageDbg( "processing residues of SSE " + ( *sse_itr)->GetIdentification());
          // iterate over the residues
          for
          (
            ; seq_aa_itr != seq_aa_itr_end && sse_aa_itr != sse_aa_itr_end;
            ++seq_aa_itr, ++sse_aa_itr
          )
          {

            BCL_MessageVrb( "\t\tAA: " + ( *sse_aa_itr)->GetIdentification());

          //////////////////////////
          // finalize description //
          //////////////////////////

            // if CA coordinates are not defined
            if( !m_TestingInput->GetFlag() && !( *sse_aa_itr)->GetCA().GetCoordinates().IsDefined())
            {
              // then skip this amino acid
              continue;
            }

            biol::EnvironmentType biol_environment( biol::GetEnvironmentTypes().e_Solution);
            // define default environment, which is true if working with JUFO9D data
            size_t merged_index( biol_environment->GetReducedIndex() * 3 + this_sstype);

            // if this pdb has the membrane flag set
            if( pdb_itr->Second())
            {
              // then we determine the environment type using the membrane object created
              biol_environment = membrane.DetermineEnvironmentType( ( *sse_aa_itr)->GetCA().GetCoordinates());

              // if testing input
              if( m_TestingInput->GetFlag() && biol_environment == biol::GetEnvironmentTypes().e_Undefined)
              {
                biol_environment = biol::GetEnvironmentTypes().e_Solution;
              }

              // true if the environment type of this amino acid is Gap
              if( biol_environment->IsGap())
              {
                // then skip this amino acid
                BCL_MessageDbg( "skipping gap AA: " + ( *seq_aa_itr)->GetIdentification());
                continue;
              }

              BCL_MessageDbg( "\t\tAA: " + biol_environment.GetName() + util::Format()( merged_index));
            } // end membrane protein

            // find the description
            storage::Map< util::SiPtr< const biol::AABase>, linal::Vector< double> >::const_iterator
              find_itr( description_map.Find( util::SiPtr< const biol::AABase>( *seq_aa_itr)));

            // make sure the description exists
            BCL_Assert
            (
              find_itr != description_map.End(),
              "can't locate description for residue" + ( *seq_aa_itr)->GetIdentification()
            );

            // create prediction vector and set the corresponding cell to 1
            linal::Vector< double> prediction_vector( number_output_nodes, double( 0.0));
            prediction_vector( merged_index) = 1.0;

            // insert it into complete dataset
            const storage::VectorND< 2, linal::Vector< double> > this_data( find_itr->second, prediction_vector);

            complete_dataset( merged_index).PushBack( this_data);
            nonrandom_complete_dataset.PushBack( this_data);
            BCL_MessageCrt( "Entered data from: " + prefix_without_chainid);
          }
        } // end itr sses
      } // end itr pdb

    ////////////////////
    // normalizations //
    ////////////////////

      BCL_MessageCrt( "Finished with descriptor generation");
      BCL_MessageCrt( "Starting normalizations");

      // initialize number of entries per state
      size_t nr_entries_per_state( 0);
      // if defined in the commandline then use that
      if( m_NumberEntriesPerStateFlag->GetFlag())
      {
        nr_entries_per_state = m_NumberEntriesPerStateParam->GetNumericalValue< size_t>();
      }

      // iterate over the vector
      for( size_t this_state( 0); this_state < number_states; ++this_state)
      {
        // and get maximum number of states
        const size_t this_size( complete_dataset( this_state).GetSize());

        BCL_MessageTop
        (
          biol::GetEnvironmentTypes().GetReducedTypes()( this_state / 3).GetName() + "-" +
          biol::SSType( this_state % 3)->GetOneLetterCode() + " ==> " + util::Format()( this_size)
        );

        // if there are no entries for this state
        if( this_size == 0)
        {
          BCL_MessageVrb
          (
            "There are no entries for this state therefore skipping! " +
            biol::GetEnvironmentTypes().GetReducedTypes()( this_state / 3).GetName() + "-" +
            biol::SSType( this_state % 3)->GetOneLetterCode()
          );
        }

        // if nr_entries_per_state flag was not given
        // find the state with the highest count and set it as the nr_entries_per_state
        if( !m_NumberEntriesPerStateFlag->GetFlag() && this_size > nr_entries_per_state)
        {
          nr_entries_per_state = this_size;
        }
      }

      BCL_MessageCrt( "number entries per state: " + util::Format()( nr_entries_per_state));

    ////////////
    // output //
    ////////////

      BCL_MessageCrt( "constructing the final list");
      // initialize final list
      storage::Vector< storage::VectorND< 2, linal::Vector< double> > > final_list;

      // for unbalanced input, no iteration of number of datasets
      if( m_UnbalancedInput->GetFlag())
      {
        // for nonrandom input, just the way the sequence is
        if( m_NonrandomInput->GetFlag())
        {
          final_list = nonrandom_complete_dataset;
        }
        else
        {
          // initialize flag to keep track if we run out of data to print
          bool elements_left( true);

          // go until we deplete all the elements
          while( elements_left)
          {
            // assume no elements are left
            elements_left = false;

            // for all states
            for( size_t nr_state( 0); nr_state < number_states; ++nr_state)
            {
              // create a reference
              storage::Vector< storage::VectorND< 2, linal::Vector< double> > > &complete_ref( complete_dataset( nr_state));

              // if the size is 0
              if( complete_ref.GetSize() == 0)
              {
                continue;
              }

              // update the elements_left flag
              elements_left = true;

              // get a random iterator
              storage::Vector< storage::VectorND< 2, linal::Vector< double> > >::iterator rand_itr
              (
                random::GetGlobalRandom().Iterator( complete_ref.Begin(), complete_ref.End(), complete_ref.GetSize())
              );

              // push the data
              final_list.PushBack( *rand_itr);

              // remove this element from the list
              complete_ref.Remove( rand_itr);
            }
          }
        }
      }
      // if balanced
      else
      {
        // runs over number of datasets per state
        for( size_t nr_entry( 0); nr_entry < nr_entries_per_state; ++nr_entry)
        {
          // for all states
          for( size_t nr_state( 0); nr_state < number_states; ++nr_state)
          {
            // create a reference
            storage::Vector< storage::VectorND< 2, linal::Vector< double> > > &complete_ref( complete_dataset( nr_state));

            // if the size is 0
            if( complete_ref.GetSize() == 0)
            {
              continue;
            }
            storage::Vector< storage::VectorND< 2, linal::Vector< double> > >::const_iterator rand_itr
            (
              random::GetGlobalRandom().Iterator( complete_ref.Begin(), complete_ref.End(), complete_ref.GetSize())
            );
            final_list.PushBack( *rand_itr);
          }
        }
      }

      BCL_MessageCrt( "outputting final list");
      // open output file
      io::OFStream write;
      io::File::MustOpenOFStream( write, m_OutputPrefixFlag->GetFirstParameter()->GetValue());
      write << final_list;
      io::File::CloseClearFStream( write);

      return 0;
    }

    //! @brief return first layer descriptors for the given amino acid sequence
    //! @param SEQUENCE AASequence of interest
    //! @param WINDOWS_MAP map of windows
    //! @return map that stores the descriptions for each
    storage::Map< util::SiPtr< const biol::AABase>, linal::Vector< double> >
    GenerateJufoDescription::CalculateFirstLayerDescriptors
    (
      util::ShPtr< biol::AASequence> &SEQUENCE,
      const size_t OLIGOMERIC_STATE
    ) const
    {
      // create the descriptor to generate the jufo dataset
      util::Implementation< descriptor::Base< biol::AABase, float> > aa_descriptor
      (
        sspred::JUFO9D::GetJufo9DANNDescriptors( OLIGOMERIC_STATE)
      );

      // set the object up
      assemble::ProteinModelWithCache model
      (
        assemble::ProteinModel( util::ShPtr< assemble::Chain>( new assemble::Chain( SEQUENCE))),
        false
      );
      aa_descriptor->SetObject( model);

      // set the dimension (1 because we operate on elements of the sequence)
      aa_descriptor->SetDimension( 1);

      // create a descriptor iterator
      descriptor::Iterator< biol::AABase> itr( model.GetIterator());

      // initialize descriptor map
      storage::Map< util::SiPtr< const biol::AABase>, linal::Vector< double> > descriptor_map;

      // iterate over the amino acids
      for( ; itr.NotAtEnd(); ++itr)
      {
        linal::Vector< float> descriptor( aa_descriptor->GetSizeOfFeatures(), ( *aa_descriptor)( itr).Begin());
        descriptor_map[ *itr( 0)] = linal::Vector< double>( descriptor.Begin(), descriptor.End());
      }

      // end
      return descriptor_map;
    }

    //! @brief initializes the command object for that executable
    //! @return initalized command object
    util::ShPtr< command::Command> GenerateJufoDescription::InitializeCommand() const
    {
      // initialize a new command object
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // path for files
      sp_cmd->AddFlag( m_PathFlag);
      // flag for reading predictions from different SSMethods
      sp_cmd->AddFlag( m_PdbListFlag);
      // flag for using pdb path hierarchy for all input files
      sp_cmd->AddFlag( m_PdbPathHierarchy);
      // path for output files
      sp_cmd->AddFlag( m_OutputPrefixFlag);
      // flag for providing a oligomeric state dictionary
      sp_cmd->AddFlag( m_OligomericStateDictionaryFlag);
      // flag afor adjusting the membrane
      sp_cmd->AddFlag( biol::Membrane::GetFlagMembrane());
      // flag to define the number of requested entries for each of the states
      sp_cmd->AddFlag( m_NumberEntriesPerStateFlag);
      //flag to create unbalanced input
      sp_cmd->AddFlag( m_UnbalancedInput);
      //flag to create nonrandom, successive input
      sp_cmd->AddFlag( m_NonrandomInput);
      //flag to create input for testing, no residues skipped
      sp_cmd->AddFlag( m_TestingInput);
      // add the possibility to convert to natural aa types
      sp_cmd->AddFlag( pdb::Factory::GetFlagConvertToNaturalAAType());
      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief default constructor
    GenerateJufoDescription::GenerateJufoDescription() :
      m_PathFlag( new command::FlagStatic( "path", "flag for setting path where input files can be found")),
      m_PathParam( new command::Parameter( "path_param", "path where input files can be found", ".")),
      m_PdbListFlag( new command::FlagStatic( "pdb_list", "flag to provide tags of pdb files to be read")),
      m_PdbListParam( new command::Parameter( "pdb_list_filename", "name of the file that has list of tags of pdbs to be read", "pdbs.ls")),
      m_PdbPathHierarchy
      (
        new command::FlagStatic
        (
          "pdb_hierarchy",
          "boolean to indicate whether a pdb hiearchy is used in input paths so for pdbtag 1ABC.pdb it looks at {path}/AB/1ABC.pdb"
        )
      ),
      m_OutputPrefixFlag
      (
        new command::FlagStatic
        (
          "output_prefix", "Flag to change the prefix to be used for the output filename",
          command::Parameter
          (
            "output_prefix_value", "\tPrefix to be used for the output filename", "descriptor.out"
          )
        )
      ),
      m_OligomericStateDictionaryFlag( new command::FlagStatic( "oligo_dict", "flag for providing a dictionary for oligomeric states( pdbtag 0/1")),
      m_OligomericStateDictionaryParam( new command::Parameter( "oligo_dict_param", "file that contains the dictionary for oligomeric states", "")),
      m_NumberEntriesPerStateFlag
      (
        new command::FlagStatic
        (
          "nr_entries_per_state", "flag to set the requested number of entries per state"
        )
      ),
      m_NumberEntriesPerStateParam
      (
        new command::Parameter( "nr_entries_per_state_param", "\tthe requested number of entries per state", "0")
      ),
      m_UnbalancedInput( new command::FlagStatic( "creating_unbalanced_input", "flag to create unbalanced input")),
      m_NonrandomInput( new command::FlagStatic( "creating_nonrandom_input", "flag to create nonrandom input")),
      m_TestingInput( new command::FlagStatic( "creating_input_for_testing", "flag to create input for testing"))
    {
      // attach parameters to flags
      m_PathFlag->PushBack( m_PathParam);
      m_PdbListFlag->PushBack( m_PdbListParam);
      m_OligomericStateDictionaryFlag->PushBack( m_OligomericStateDictionaryParam);
      m_NumberEntriesPerStateFlag->PushBack( m_NumberEntriesPerStateParam);
    }

    const ApplicationType GenerateJufoDescription::GenerateJufoDescription_Instance
    (
      GetAppGroups().AddAppToGroup( new GenerateJufoDescription(), GetAppGroups().e_InternalBiol)
    );

  } // namespace app
} // namespace bcl
