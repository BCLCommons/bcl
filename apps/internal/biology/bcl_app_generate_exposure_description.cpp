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
#include "assemble/bcl_assemble_aa_neighbor_count.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container_generator_protein_model.h"
#include "biol/bcl_biol_aa.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "biol/bcl_biol_blast_profile_handler.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "descriptor/bcl_descriptor_base.h"
#include "descriptor/bcl_descriptor_iterator.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "sspred/bcl_sspred_jufo9d.h"
#include "util/bcl_util_implementation.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class GenerateExposureDescription
    //! @brief generates training data for predicting AA exposure (SASA) from sequence
    //!
    //! @author weinerbe
    //! @date Aug 30, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API GenerateExposureDescription :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      //! flag for input file
      util::ShPtr< command::FlagInterface> m_InputFileFlag;

      //! flag for output file
      util::ShPtr< command::FlagInterface> m_OutputFileFlag;

      //! flag for providing a oligomeric state dictionary
      util::ShPtr< command::FlagStatic> m_OligomericStateDictionaryFlag;
      util::ShPtr< command::ParameterInterface> m_OligomericStateDictionaryParam;

      //! flag for printing out NC and rSASA values instead of descriptors
      util::ShPtr< command::FlagInterface> m_NCFlag;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      GenerateExposureDescription();

      //! @brief Clone function
      //! @return pointer to new GenerateExposureDescription
      GenerateExposureDescription *Clone() const
      {
        return new GenerateExposureDescription( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief initializes the command object for that executable
      //! @return initalized command object
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief Main function
      //! @return return value of the application
      int Main() const;

      //! @brief Reads in SASA file
      //! @return map of seq id to SASA
      static storage::Map< size_t, double> ReadSASAFile( const std::string &FILENAME);

      //! @brief compute rSASA
      //! @param SASA_VALUE un-normalized sasa value
      //! @param AA_TYPE amino acid type
      //! @return rSASA
      static double ComputeRSASA( const double &SASA_VALUE, const biol::AAType &AA_TYPE);

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
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // end
        return OSTREAM;
      }

    private:

      //! static instance of this class
      static const ApplicationType GenerateExposureDescription_Instance;

    }; // class GenerateExposureDescription

    //! @brief Main function
    //! @return return value of the application
    int GenerateExposureDescription::Main() const
    {
      // initialize list of descriptors and features to train
      storage::Vector< storage::VectorND< 2, linal::Vector< double> > > final_list;

      // initialize neighbor count stats
      storage::Vector< storage::VectorND< 2, double> > nc_rsasa_values;

      // instantiate an undefined amino acid with 0 blast profile
      biol::AA undefined_aa( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().XXX)));
      undefined_aa.SetBlastProfile( biol::BlastProfile( linal::Vector< double>( 20, 0.0)));

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
          oligo_dictionary[ tag.substr( 0, 4)] = oligo_state;
        }
        io::File::CloseClearFStream( oligo_read);
      }

      // read in the input file
      io::IFStream read;
      io::File::MustOpenIFStream( read, m_InputFileFlag->GetFirstParameter()->GetValue());

      // iterate over the file
      while( !read.eof())
      {
        // read the line
        std::string this_line;
        std::getline( read, this_line);
        util::TrimString( this_line);
        storage::Vector< std::string> line_entries( util::SplitString( this_line));

        // check the size
        if( line_entries.IsEmpty())
        {
          continue;
        }
        BCL_Assert( line_entries.GetSize() == 3, "Each line in input file should have 3 entries");

        // set up fasta
        const std::string &fasta_filename( line_entries( 0));

        // read in per residue sasa file
        const storage::Map< size_t, double> sasa_map( ReadSASAFile( line_entries( 1)));

        // if the nc flag is passed
        if( m_NCFlag->GetFlag())
        {
          // get pdb file
          std::string pdb_filename( fasta_filename);
          pdb_filename.replace( fasta_filename.find( ".fasta"), 7, ".pdb");
          BCL_MessageStd( "Reading PDB file: " + pdb_filename);

          // open pdb
          const pdb::Factory factory;
          const assemble::ProteinModel protein_model( factory.ProteinModelFromPDBFilename( pdb_filename));
          BCL_Assert( protein_model.GetChains().GetSize() == 1, "Protein model should only have one chain");

          // create AA neighbor list generator
          const assemble::AANeighborCount neighbor_count;

          const util::ShPtr< math::FunctionInterfaceSerializable< assemble::ProteinModel, assemble::AANeighborListContainer> >
          sp_generator
          (
            assemble::AANeighborListContainerGeneratorProteinModel::AANeighborListGenerator
            (
              neighbor_count.GetThresholdRange().GetMax(),
              neighbor_count.GetMinimalSequenceSeparation(),
              false,
              false
            )
          );

          // get the container
          const assemble::AANeighborListContainer neighbor_container( sp_generator->operator ()( protein_model));

          // iterate over the residues
          for
          (
            assemble::AANeighborListContainer::const_iterator itr( neighbor_container.Begin()),
              itr_end( neighbor_container.End());
            itr != itr_end; ++itr
          )
          {
            // see if this residue has a SASA value
            storage::Map< size_t, double>::const_iterator find_itr( sasa_map.Find( itr->first->GetSeqID()));
            if( find_itr != sasa_map.End())
            {
              nc_rsasa_values.PushBack
              (
                storage::VectorND< 2, double>
                (
                  neighbor_count( itr->second),
                  ComputeRSASA( find_itr->second, itr->first->GetType())
                )
              );
            }
          }

          continue;
        }

        // read in the fasta
        io::IFStream input_stream;
        io::File::MustOpenIFStream( input_stream, fasta_filename);
        BCL_MessageStd( "Reading fasta " + fasta_filename);
        biol::AASequence this_sequence( biol::AASequenceFactory::BuildSequenceFromFASTA( input_stream));
        io::File::CloseClearFStream( input_stream);

        // set the oligo state to 0
        size_t this_oligo_state( 0);

        // if the tag is in the map look up the oligostate for this pdb and store it
        const std::string pdb_tag( fasta_filename.substr( fasta_filename.find( ".fasta") - 5, 4));
        if( oligo_dictionary.Find( pdb_tag) != oligo_dictionary.End())
        {
          this_oligo_state = oligo_dictionary[ pdb_tag];
        }

        // read in the blast profile
        std::string blast_filename( fasta_filename);
        blast_filename.replace( fasta_filename.find( ".fasta"), 7, ".ascii6");
        io::File::MustOpenIFStream( input_stream, blast_filename);
        biol::BlastProfileHandler::ReadProfileForAASequence( input_stream, this_sequence);
        io::File::CloseClearFStream( input_stream);

        // calculate the windows map
        const storage::Map< util::SiPtr< const biol::AABase>, util::SiPtrVector< const biol::AABase> > windows_map
        (
          biol::CreateWindowsFromAminoAcids( this_sequence, 15, undefined_aa)
        );

        // create the descriptor to generate the jufo dataset
        util::Implementation< descriptor::Base< biol::AABase, float> > aa_descriptor
        (
          sspred::JUFO9D::GetJufo9DANNDescriptors( this_oligo_state)
        );

        // set the object up
        assemble::ProteinModelWithCache model
        (
          assemble::ProteinModel
          (
            util::ShPtr< assemble::Chain>( new assemble::Chain( util::CloneToShPtr( this_sequence)))
          ),
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

        // iterate over the sequence
        for
        (
          biol::AASequence::const_iterator aa_itr( this_sequence.Begin()), aa_itr_end( this_sequence.End());
          aa_itr != aa_itr_end; ++aa_itr
        )
        {
          // find the description
          storage::Map< util::SiPtr< const biol::AABase>, linal::Vector< double> >::const_iterator
            find_itr( descriptor_map.Find( util::SiPtr< const biol::AABase>( *aa_itr)));

          // make sure the description exists
          BCL_Assert
          (
            find_itr != descriptor_map.End(),
            "can't locate description for residue" + ( *aa_itr)->GetIdentification()
          );

          // find the sasa value for this residue
          const storage::Map< size_t, double>::const_iterator sasa_find_itr( sasa_map.Find( ( *aa_itr)->GetSeqID()));

          // if it was found
          if( sasa_find_itr != sasa_map.End())
          {
            // insert it into complete dataset
            const linal::Vector< double> prediction_vector
            (
              1,
              ComputeRSASA( sasa_find_itr->second, ( *aa_itr)->GetType())
            );
            const storage::VectorND< 2, linal::Vector< double> > this_data( find_itr->second, prediction_vector);
            final_list.PushBack( this_data);
          }
        }
      }

      // close the stream
      io::File::CloseClearFStream( read);

      // write the data
      io::OFStream write;
      io::File::MustOpenOFStream( write, m_OutputFileFlag->GetFirstParameter()->GetValue());
      if( m_NCFlag->GetFlag())
      {
        for
        (
          storage::Vector< storage::VectorND< 2, double> >::const_iterator itr( nc_rsasa_values.Begin()),
            itr_end( nc_rsasa_values.End());
          itr != itr_end; ++itr
        )
        {
          write << itr->First() << '\t' << itr->Second() << '\n';
        }
      }
      else
      {
        write << final_list;
      }
      io::File::CloseClearFStream( write);

      // end
      return 0;
    }

    //! @brief Reads in SASA file
    //! @return map of seq id to SASA
    storage::Map< size_t, double> GenerateExposureDescription::ReadSASAFile( const std::string &FILENAME)
    {
      // initialize map
      storage::Map< size_t, double> sasa_map;

      // open file
      io::IFStream read;
      io::File::MustOpenIFStream( read, FILENAME);

      // initialize variables
      size_t seq_id;
      double sasa;

      // iterate over the file
      while( read >> seq_id >> sasa && !read.eof())
      {
        // update the map
        sasa_map[ seq_id] = sasa;
      }

      // end
      return sasa_map;
    }

    //! @brief compute rSASA
    //! @param SASA_VALUE un-normalized sasa value
    //! @param AA_TYPE amino acid type
    //! @return rSASA
    double GenerateExposureDescription::ComputeRSASA( const double &SASA_VALUE, const biol::AAType &AA_TYPE)
    {
      // calculate and return
      return SASA_VALUE / AA_TYPE->GetAAProperty( biol::AATypeData::e_SASA);
    }

    //! @brief initializes the command object for that executable
    //! @return initalized command object
    util::ShPtr< command::Command> GenerateExposureDescription::InitializeCommand() const
    {
      // initialize a new command object
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // add members
      sp_cmd->AddFlag( m_InputFileFlag);
      sp_cmd->AddFlag( m_OutputFileFlag);
      sp_cmd->AddFlag( m_OligomericStateDictionaryFlag);
      sp_cmd->AddFlag( m_NCFlag);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief default constructor
    GenerateExposureDescription::GenerateExposureDescription() :
      m_InputFileFlag
      (
        new command::FlagStatic
        (
          "input",
          "\tinput file to be used. Each line has the form, "
            "[fasta file] [msms.area file] [Membrane/Soluble]",
          command::Parameter( "input_file", "path and name of input file", "")
        )
      ),
      m_OutputFileFlag
      (
        new command::FlagStatic
        (
          "output",
          "\toutput file to be created containing training data",
          command::Parameter( "output_file", "path and name of output file", "exposure.data")
        )
      ),
      m_OligomericStateDictionaryFlag
      (
        new command::FlagStatic( "oligo_dict", "flag for providing a dictionary for oligomeric states( pdbtag 0/1")
      ),
      m_OligomericStateDictionaryParam
      (
        new command::Parameter( "oligo_dict_param", "file that contains the dictionary for oligomeric states", "")
      ),
      m_NCFlag( new command::FlagStatic( "nc", "flag for printing correlation between neighbor count and rSASA"))
    {
      m_OligomericStateDictionaryFlag->PushBack( m_OligomericStateDictionaryParam);
    }

    const ApplicationType GenerateExposureDescription::GenerateExposureDescription_Instance
    (
      GetAppGroups().AddAppToGroup( new GenerateExposureDescription(), GetAppGroups().e_InternalBiol)
    );

  } // namespace app
} // namespace bcl
