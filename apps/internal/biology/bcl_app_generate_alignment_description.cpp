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
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "storage/bcl_storage_template_instantiations.h"

namespace bcl
{
  namespace app
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class GenerateAlignmentDescription
    //! @brief TODO: add brief class comment
    //!
    //! @author karakam
    //! @date 03/01/2010
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API GenerateAlignmentDescription :
      public Interface
    {
    private:

    //////////
    // data //
    //////////

      //! name of the file that has the list of alignment files
      util::ShPtr< command::ParameterInterface> m_AlignmentListParam;

      //! window size
      util::ShPtr< command::FlagStatic> m_WindowSizeFlag;
      util::ShPtr< command::ParameterInterface> m_WindowSizeParam;

      //! path to the directory where the blast profiles reside
      util::ShPtr< command::FlagStatic> m_BlastPathFlag;
      util::ShPtr< command::ParameterInterface> m_BlastPathParam;

      //! path to the directory where the SS predictions reside
      util::ShPtr< command::FlagStatic> m_SSPredPathFlag;
      util::ShPtr< command::ParameterInterface> m_SSPredPathParam;

      //! flag to identify which SSPredictions to use
      util::ShPtr< command::FlagInterface> m_SSPredMethodsFlag;

      //! flag to identify the ratio of positive pairs to be used
      util::ShPtr< command::FlagInterface> m_PositiveRatio;

      //! flag for setting the output filename
      util::ShPtr< command::FlagInterface> m_OutputFilenameFlag;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      GenerateAlignmentDescription();

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      GenerateAlignmentDescription *Clone() const
      {
        return new GenerateAlignmentDescription( *this);
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
      static const ApplicationType GenerateAlignmentDescription_Instance;

    }; // class GenerateAlignmentDescription

    //! @brief Main function
    //! @return return value of the application
    int GenerateAlignmentDescription::Main() const
    {

//      // initialize stream
//      io::IFStream read;
//      io::File::MustOpenIFStream( read, m_AlignmentListParam->GetValue());
//
//      // get the sspred methods set
//      storage::Set< sspred::Method> ss_methods;
//      storage::Vector< sspred::Method> methods_commandline( m_SSPredMethodsFlag->GetObjectList< sspred::Method>());
//      ss_methods.InsertElements( methods_commandline.Begin(), methods_commandline.End());
//
//      // introduce typedef
//      typedef storage::Pair< util::ShPtrList< biol::AASequence>, util::ShPtr< align::AlignmentSimple< biol::AABase> > >
//       SequencesAndAlignmentPair;
//
//      // store the provided window size
//      const size_t window_size( m_WindowSizeParam->GetNumericalValue< size_t>());
//      BCL_Assert( window_size % 2 == 1, "The window size should be an odd number not " + util::Format()( window_size));
//      // calculate the radius from the windows size
//      const size_t window_radius( ( window_size - 1) / 2);
//
//      // initialize the list to store all alignments
//      storage::List< SequencesAndAlignmentPair> alignments_list;
//
//      // list of sequences
//      storage::Map< std::string, util::ShPtr< biol::AASequence> > sequences_list;
//
//      // create a container for storing the residue pairs that are aligned for every sequence pair
//      storage::Map< std::string, storage::Map< std::string, storage::List< storage::VectorND< 2, size_t> > > >
//        aligned_aa_pairs;
//
//    ///////////////////////////////////////////////////
//    // reading alignments, blast profiles and sspred //
//    ///////////////////////////////////////////////////
//
//      // read the name of string files
//      const storage::Vector< std::string> list_of_alignment_filenames( util::StringListFromIStream( read));
//      io::File::CloseClearFStream( read);
//
//      // iterate over the each filename
//      for
//      (
//        storage::Vector< std::string>::const_iterator filename_itr( list_of_alignment_filenames.Begin()),
//          filename_itr_end( list_of_alignment_filenames.End());
//        filename_itr != filename_itr_end; ++filename_itr
//      )
//      {
//        // attach the ifstream to the filename
//        io::File::MustOpenIFStream( read, *filename_itr);
//
//        // read the alignment from file and pushback it into the alignment list
//        SequencesAndAlignmentPair this_sequence_alignment_pair
//        (
//          align::HandlerPIR< biol::AABase>().ReadAlignment( read)
//        );
//
//        // clear the stream
//        io::File::CloseClearFStream( read);
//
//        // iterate over the sequence in this alignment
//        for
//        (
//          util::ShPtrList< biol::AASequence>::iterator seq_itr( this_sequence_alignment_pair.First().Begin()),
//            seq_itr_end( this_sequence_alignment_pair.First().End());
//          seq_itr != seq_itr_end; ++seq_itr
//        )
//        {
//          // get the fasta name for this sequence
//          const std::string fasta_header( ( *seq_itr)->GetFastaHeader());
//          const std::string fasta_name( fasta_header.substr( 1, fasta_header.length() - 1));
//
//          // open the stream to blast profile
//          io::File::MustOpenIFStream( read, m_BlastPathParam->GetValue() + PATH_SEPARATOR + fasta_name + ".ascii6");
//
//          // read blast profile and clear the stream
//          biol::BlastProfileHandler::ReadProfileForAASequence( read, **seq_itr);
//          io::File::CloseClearFStream( read);
//
//          // now read the secondary structure predictions
//          sspred::MethodHandler::ReadPredictionsForAASequence
//          (
//            ss_methods, **seq_itr, fasta_name, m_SSPredPathParam->GetValue()
//          );
//
//          // insert this sequence into sequences_list
//          sequences_list[ fasta_header] = *seq_itr;
//        }
//
//        // push this sequence alignment pair to the list
//        alignments_list.PushBack( this_sequence_alignment_pair);
//      }
//
//    //////////////////////////////
//    // finding aligned residues //
//    //////////////////////////////
//
//      // initialize counter to store the number of aligned pairs
//      size_t nr_aligned_pairs( 0);
//
//      BCL_MessageStd( "Looking for aligned residues");
//      // start iterating over the read alignments
//      for
//      (
//        storage::List< SequencesAndAlignmentPair>::const_iterator
//          seq_align_itr( alignments_list.Begin()), seq_align_itr_end( alignments_list.End());
//        seq_align_itr != seq_align_itr_end; ++seq_align_itr
//      )
//      {
//        // store the number of sequences
//        const size_t number_sequences( seq_align_itr->First().GetSize());
//
//        // iterate over the sequences and store the fasta names as a vector, as well as the sizes
//        storage::Vector< std::string> fasta_names_vector;
//        storage::Vector< size_t> sequence_length_vector;
//        for
//        (
//          util::ShPtrList< biol::AASequence>::const_iterator seq_itr( seq_align_itr->First().Begin()),
//          seq_itr_end( seq_align_itr->First().End());
//          seq_itr != seq_itr_end; ++seq_itr
//        )
//        {
//          // insert this fasta name into the vector
//          fasta_names_vector.PushBack( ( *seq_itr)->GetFastaHeader());
//          sequence_length_vector.PushBack( ( *seq_itr)->GetSize());
//        }
//
//        // iterate over the assignments in the alignment
//        for
//        (
//          util::ShPtrList< align::Assignment< biol::AABase> >::const_iterator
//            alignment_itr( seq_align_itr->Second()->GetAssignments().Begin()),
//            alignment_itr_end( seq_align_itr->Second()->GetAssignments().End());
//          alignment_itr != alignment_itr_end;
//          ++alignment_itr
//        )
//        {
//          // create a reference on the assignment
//          const align::Assignment< biol::AABase> &assignment( **alignment_itr);
//
//          // iterate over the members of the assignment
//          for( size_t index_a( 0); index_a < number_sequences; ++index_a)
//          {
//            // iterate over the members of the assignment
//            for( size_t index_b( index_a + 1); index_b < number_sequences; ++index_b)
//            {
//              // create references on the members
//              const util::SiPtr< const biol::AABase> sp_amino_acid_a( *assignment.GetNthMember( index_a));
//              const util::SiPtr< const biol::AABase> sp_amino_acid_b( *assignment.GetNthMember( index_b));
//
//              // if either of amino acids are not defined
//              if( !sp_amino_acid_a.IsDefined() || !sp_amino_acid_b.IsDefined())
//              {
//                // skip
//                continue;
//              }
//
//              // if either of amino acids are in the beginning or ending part of sequence
//              if
//              (
//                size_t( sp_amino_acid_a->GetSeqID()) <= window_radius ||
//                sp_amino_acid_a->GetSeqID() + window_radius > sequence_length_vector( index_a) ||
//                size_t( sp_amino_acid_b->GetSeqID()) <= window_radius ||
//                sp_amino_acid_b->GetSeqID() + window_radius > sequence_length_vector( index_b)
//              )
//              {
//                // skipm_OutputFilenameFlag->GetFirstParameter()->GetValue()
//                continue;
//              }
//
//              BCL_Message
//              (
//                util::Message::e_Verbose,
//                "aligned residue from sequences: : " + fasta_names_vector( index_a) + " and " +
//                fasta_names_vector( index_b) + sp_amino_acid_a->GetIdentification() + " vs " +
//                sp_amino_acid_b->GetIdentification()
//              );
//
//              // insert the aligned aa information into corresponding list
//              aligned_aa_pairs[ fasta_names_vector( index_a)][ fasta_names_vector( index_b)].PushBack
//              (
//                storage::VectorND< 2, size_t>( sp_amino_acid_a->GetSeqID(), sp_amino_acid_b->GetSeqID())
//              );
//              // increment
//              ++nr_aligned_pairs;
//            }
//          }
//        } // iterate over assignments
//      } // iterate over all sequence and alignments
//
//    ////////////////////////////
//    // generate negative pairs//
//    ////////////////////////////
//
//      // get the positive ratio
//      const double positive_ratio( m_PositiveRatio->GetFirstParameter()->GetNumericalValue< double>());
//      const double negative_ratio( 1.0 - positive_ratio);
//      const size_t nr_positives( nr_aligned_pairs);
//      const size_t nr_total( nr_aligned_pairs / positive_ratio);
//      const size_t nr_negatives( nr_total - nr_positives);
//      BCL_Message
//      (
//        util::Message::e_Standard,
//        "positive ratio: " + util::Format()( positive_ratio) +
//        "\nnegative ratio: " + util::Format()( negative_ratio) +
//        "\nnr total: " + util::Format()( nr_total) +
//        "\nnr positives: " + util::Format()( nr_positives) +
//        "\nnr negatives: " + util::Format()( nr_negatives)
//      );
//
//      // create a container for storing the residue pairs that are negative examples
//      storage::Map< std::string, storage::Map< std::string, storage::List< storage::VectorND< 2, size_t> > > >
//        unaligned_aa_pairs;
//
//      BCL_MessageStd( "Creating negative results");
//      // iterate until there are enough negatives
//      for( size_t ctr_negative( 0); ctr_negative < nr_negatives; ++ctr_negative)
//      {
//        // get a random element on the alignments list
//        storage::List< SequencesAndAlignmentPair>::const_iterator align_itr
//        (
//          random::GetGlobalRandom().Iterator( alignments_list.Begin(), alignments_list.End(), alignments_list.GetSize())
//        );
//        // now pick two random sequences from the alignment
//        util::SiPtrVector< const biol::AASequence> seq_vector( align_itr->First().Begin(), align_itr->First().End());
//
//        // pick two random sequence indices
//        const size_t seq_index_a( random::GetGlobalRandom().Random< size_t>( seq_vector.GetSize() - 1));
//        size_t seq_index_b( seq_index_a);
//        while( seq_index_a == seq_index_b)
//        {
//          seq_index_b = random::GetGlobalRandom().Random< size_t>( seq_vector.GetSize() - 1);
//        }
//
//        // create references to the strings for the sequence fasta headers a and b
//        util::SiPtr< const biol::AASequence> seq_a( seq_vector( seq_index_a));
//        util::SiPtr< const biol::AASequence> seq_b( seq_vector( seq_index_b));
//        const std::string &fasta_a( seq_a->GetFastaHeader());
//        const std::string &fasta_b( seq_b->GetFastaHeader());
//
//        // aligned pairs in these sequences
//        storage::List< storage::VectorND< 2, size_t> >& this_aligned_pairs
//        (
//          aligned_aa_pairs[ fasta_a][ fasta_b]
//        );
//
//        // initiate bool for success
//        bool success( false);
//
//        // while no success
//        while( !success)
//        {
//          // generate two random numbers
//          storage::VectorND< 2, size_t> random_seqid_pair
//          (
//            random::GetGlobalRandom().Random< size_t>( window_radius + 1, seq_a->GetSize() - window_radius),
//            random::GetGlobalRandom().Random< size_t>( window_radius + 1, seq_b->GetSize() - window_radius)
//          );
//
//          // check that these two seqids are not already defined as aligned
//          storage::List< storage::VectorND< 2, size_t> >::const_iterator pair_itr
//          (
//            std::find( this_aligned_pairs.Begin(), this_aligned_pairs.End(), random_seqid_pair)
//          );
//
//          // if the pair itr was not found
//          if( pair_itr == this_aligned_pairs.End())
//          {
//            // insert it into the list
//            unaligned_aa_pairs[ fasta_a][ fasta_b].PushBack( random_seqid_pair);
//
//            BCL_Message
//            (
//              util::Message::e_Verbose,
//              "Negative example from " + fasta_a + " and " + fasta_b + " " +
//              util::Format()( random_seqid_pair.First()) + " vs " + util::Format()( random_seqid_pair.Second())
//            )
//
//            // set success to true
//            success = true;
//          }
//        }
//      }
//
//    ///////////////////////////////////
//    // initialize descriptor objects //
//    ///////////////////////////////////
//
//      BCL_MessageStd( "constructing descriptors");
//
//    ////////////////////
//    // AA descriptors //
//    ////////////////////
//
//      BCL_MessageStd( "Creating aa descriptors");
//
//      // form the property list
//      storage::List< biol::AATypeData::PropertyType> aa_properties_list;
//      aa_properties_list.PushBack( biol::AATypeData::e_StericalParameter);
//      aa_properties_list.PushBack( biol::AATypeData::e_Polarizability);
//      aa_properties_list.PushBack( biol::AATypeData::e_Volume);
//      aa_properties_list.PushBack( biol::AATypeData::e_IsoelectricPoint);
//      aa_properties_list.PushBack( biol::AATypeData::e_SASA);
//      aa_properties_list.PushBack( biol::AATypeData::e_FreeEnergyHelix);
//      aa_properties_list.PushBack( biol::AATypeData::e_FreeEnergyStrand);
//      aa_properties_list.PushBack( biol::AATypeData::e_FreeEnergyCoil);
//      aa_properties_list.PushBack( biol::AATypeData::e_FreeEnergyCore);
//      aa_properties_list.PushBack( biol::AATypeData::e_FreeEnergyTransition);
//      aa_properties_list.PushBack( biol::AATypeData::e_FreeEnergySolution);
//      aa_properties_list.PushBack( biol::AATypeData::e_FreeEnergyCoreHelix);
//      aa_properties_list.PushBack( biol::AATypeData::e_FreeEnergyTransitionHelix);
//      aa_properties_list.PushBack( biol::AATypeData::e_FreeEnergySolutionHelix);
//      aa_properties_list.PushBack( biol::AATypeData::e_FreeEnergyCoreStrand);
//      aa_properties_list.PushBack( biol::AATypeData::e_FreeEnergyTransitionStrand);
//      aa_properties_list.PushBack( biol::AATypeData::e_FreeEnergySolutionStrand);
//      aa_properties_list.PushBack( biol::AATypeData::e_FreeEnergyCoreCoil);
//      aa_properties_list.PushBack( biol::AATypeData::e_FreeEnergyTransitionCoil);
//      aa_properties_list.PushBack( biol::AATypeData::e_FreeEnergySolutionCoil);
//
//      // create the AA descriptors subset that will also be used in global descriptors
//      util::ShPtr< descriptor::List< biol::AABase> >  descriptors_aa_subset( new descriptor::List< biol::AABase>());
//      // properties
//      descriptors_aa_subset->PushBack
//      (
//        util::ShPtr< descriptor::Interface< biol::AABase> >( new descriptor::AAPropertyOld( aa_properties_list))
//      );
//      // blast profile
//      descriptors_aa_subset->PushBack
//      (
//        util::ShPtr< descriptor::Interface< biol::AABase> >( new descriptor::AABlastProfileOld())
//      );
//
//      // initialize aa descriptor list that will contain all aa descriptor not only global ones
//      util::ShPtr< descriptor::List< biol::AABase> > descriptors_aa( descriptors_aa_subset->Clone());
//      // iterate over the ssmethods
//      for
//      (
//        storage::Set< sspred::Method>::const_iterator
//          method_itr( ss_methods.Begin()), method_itr_end( ss_methods.End());
//        method_itr != method_itr_end; ++method_itr
//      )
//      {
//        // create three_state boolean
//        const bool three_state( *method_itr != sspred::GetMethods().e_JUFO16D);
//
//        // create the descriptor
//        descriptors_aa->PushBack
//        (
//          util::ShPtr< descriptor::Interface< biol::AABase> >
//          (
//            new descriptor::AASSTMPrediction( *method_itr, three_state)
//          )
//        );
//      }
//
//    ////////////////////////
//    // window descriptors //
//    ////////////////////////
//
//      BCL_MessageStd( "Creating window descriptors");
//
//      // create the window descriptor
//      descriptor::WindowOld< biol::AABase> descriptor_window( descriptors_aa);
//
//      // create the window paired descriptor
//      util::ShPtr< descriptor::List< storage::VectorND< 2, util::SiPtr< const biol::AABase> > > > paired_descriptors
//      (
//        new descriptor::List< storage::VectorND< 2, util::SiPtr< const biol::AABase> > >()
//      );
//      // BLOSUM
//      paired_descriptors->PushBack
//      (
//        util::ShPtr< descriptor::Interface< storage::VectorND< 2, util::SiPtr< const biol::AABase> > > >
//        (
//          new descriptor::AAPairProbability
//          (
//            score::AAAssignmentBLOSUM::GetBLOSUMMatrix( score::AAAssignmentBLOSUM::e_BLOSUM_90)
//          )
//        )
//      );
//      // PAM
//      paired_descriptors->PushBack
//      (
//        util::ShPtr< descriptor::Interface< storage::VectorND< 2, util::SiPtr< const biol::AABase> > > >
//        (
//          new descriptor::AAPairProbability
//          (
//            score::AAAssignmentPAM::GetPAMMatrix( score::AAAssignmentPAM::e_PAM_250)
//          )
//        )
//      );
//      // PHAT
//      paired_descriptors->PushBack
//      (
//        util::ShPtr< descriptor::Interface< storage::VectorND< 2, util::SiPtr< const biol::AABase> > > >
//        (
//          new descriptor::AAPairProbability
//          (
//            score::AAAssignmentPHAT::GetPHATMatrix( score::AAAssignmentPHAT::e_PHAT_85)
//          )
//        )
//      );
//      // paired window descriptor
//      descriptor::WindowPaired< biol::AABase> descriptor_window_paired( paired_descriptors);
//
//    //////////////////////////////////////
//    // precalculate global descriptions //
//    //////////////////////////////////////
//
//      BCL_MessageStd( "Pre-calculating global descriptors");
//
//      // initialize global description list
//      descriptor::List< biol::AASequence> global_sequence_descriptor;
//      // insert global sequence average descriptor
//      global_sequence_descriptor.PushBack
//      (
//        util::ShPtr< descriptor::Interface< biol::AASequence> >
//        (
//          new descriptor::AASequenceAverage( descriptors_aa_subset)
//        )
//      );
//      // insert chain length
//      global_sequence_descriptor.PushBack
//      (
//        util::ShPtr< descriptor::Interface< biol::AASequence> >
//        (
//          new descriptor::AASequenceLength()
//        )
//      );
//
//      // initialize place to store the average descriptions for each sequence
//      storage::Map< std::string, storage::List< double> > sequence_average_descriptions;
//
//      // iterate over the sequences
//      for
//      (
//        storage::Map< std::string, util::ShPtr< biol::AASequence> >::const_iterator
//          seq_itr( sequences_list.Begin()), seq_itr_end( sequences_list.End());
//        seq_itr != seq_itr_end; ++seq_itr
//      )
//      {
//        BCL_MessageStd( "Pre-calculating global descriptors for " + seq_itr->first);
//
//        // generate description and insert it into the map
//        sequence_average_descriptions[ seq_itr->first] = global_sequence_descriptor( *seq_itr->second);
//      }
//
//    ///////////////////////////
//    // generate descriptions //
//    ///////////////////////////
//
//      // initialize storage for the description of true positives
//      storage::List< storage::VectorND< 2, linal::Vector< double> > > all_descriptions;
//
//      BCL_MessageStd( "generating descriptions");
//
//      // iterate over the first fasta name aligned residues map
//      for
//      (
//        storage::Map< std::string, storage::Map< std::string, storage::List< storage::VectorND< 2, size_t> > > >::const_iterator
//          seq_a_itr( aligned_aa_pairs.Begin()), seq_a_itr_end( aligned_aa_pairs.End());
//        seq_a_itr != seq_a_itr_end; ++seq_a_itr
//      )
//      {
//        // store the name of the fasta
//        const std::string fasta_name_a( seq_a_itr->first);
//        // create a reference on the corresponding AASequence
//        const biol::AASequence &sequence_a( *sequences_list[ fasta_name_a]);
//        // store the global descriptions for this sequence
//        storage::List< double> sequence_description_a( sequence_average_descriptions[ fasta_name_a]);
//
//        // iterate over the second fasta name
//        for
//        (
//          storage::Map< std::string, storage::List< storage::VectorND< 2, size_t> > >::const_iterator
//            seq_b_itr( seq_a_itr->second.Begin()), seq_b_itr_end( seq_a_itr->second.End());
//          seq_b_itr != seq_b_itr_end; ++seq_b_itr
//        )
//        {
//          // store the name of the fasta
//          const std::string fasta_name_b( seq_b_itr->first);
//          // create a reference on the corresponding AASequence
//          const biol::AASequence &sequence_b( *sequences_list[ fasta_name_b]);
//          // store the global descriptions for this sequence
//          storage::List< double> sequence_description_b( sequence_average_descriptions[ fasta_name_b]);
//
//          BCL_MessageStd( "generating descriptions for " + fasta_name_a + " and " + fasta_name_b);
//          // form a list to store both aligned and unaligned residues
//          storage::List< storage::Pair< storage::VectorND< 2, size_t>, bool> > all_pairs;
//
//          // create a reference on the unaligned pairs
//          storage::List< storage::VectorND< 2, size_t> > &
//            this_unaligned_pairs( unaligned_aa_pairs[ fasta_name_a][ fasta_name_b]);
//
//          // iterate over the aligned pairs
//          for
//          (
//            storage::List< storage::VectorND< 2, size_t> >::const_iterator
//              pair_itr( seq_b_itr->second.Begin()), pair_itr_end( seq_b_itr->second.End());
//            pair_itr != pair_itr_end; ++pair_itr
//          )
//          {
//            // add them to all pairs
//            all_pairs.PushBack( storage::Pair< storage::VectorND< 2, size_t>, bool>( *pair_itr, true));
//          }
//          // iterate over the unaligned residues
//          for
//          (
//            storage::List< storage::VectorND< 2, size_t> >::const_iterator
//              pair_itr( this_unaligned_pairs.Begin()), pair_itr_end( this_unaligned_pairs.End());
//            pair_itr != pair_itr_end; ++pair_itr
//          )
//          {
//            // add them to all pairs
//            all_pairs.PushBack( storage::Pair< storage::VectorND< 2, size_t>, bool>( *pair_itr, false));
//          }
//
//          // now iterate over all pairs
//          for
//          (
//            storage::List< storage::Pair< storage::VectorND< 2, size_t>, bool> >::const_iterator
//              pair_itr( all_pairs.Begin()), pair_itr_end( all_pairs.End());
//            pair_itr != pair_itr_end; ++pair_itr
//          )
//          {
//            // store the seq_ids
//            const size_t seqid_a( pair_itr->First().First());
//            const size_t seqid_b( pair_itr->First().Second());
//
//            BCL_Message
//            (
//              util::Message::e_Verbose,
//              sequence_a.GetAA( seqid_a)->GetIdentification() + " vs" + sequence_b.GetAA( seqid_b)->GetIdentification()
//            );
//
//            // generate windows
//            util::SiPtrVector< const biol::AABase>
//              window_a( sequence_a.GetData().SubShPtrVector( seqid_a - window_radius - 1, window_size));
//            util::SiPtrVector< const biol::AABase>
//              window_b( sequence_b.GetData().SubShPtrVector( seqid_b - window_radius - 1, window_size));
//            // generate window pair
//            const storage::VectorND< 2, util::SiPtrVector< const biol::AABase> >
//              window_pair( window_a, window_b);
//
//            // generate descriptions
//            storage::List< double> this_description;
//            this_description.Append( descriptor_window( window_a));
//            this_description.Append( double( seqid_a));
//            this_description.Append( descriptor_window( window_b));
//            this_description.Append( double( seqid_b));
//            this_description.Append( descriptor_window_paired( window_pair));
//            this_description.Append( sequence_description_a);
//            this_description.Append( sequence_description_b);
//            BCL_Message
//            (
//              util::Message::e_Debug,
//              "Description \n" + util::Format()( this_description) + " ==> " + util::Format()( pair_itr->Second())
//            );
//
//            // insert this descriptions into map of list of all descriptions
//            all_descriptions.PushBack
//            (
//              storage::VectorND< 2, linal::Vector< double> >
//              (
//                linal::Vector< double>( this_description.Begin(), this_description.End()),
//                linal::Vector< double>( 1, double( pair_itr->Second()))
//              )
//            );
//          }
//        }
//      }
//
//      // now output all the descriptions to a file
//      BCL_Message
//      (
//        util::Message::e_Standard,
//        "outputting the description to file " + m_OutputFilenameFlag->GetFirstParameter()->GetValue()
//      );
//      io::OFStream write;
//      io::File::MustOpenOFStream( write, m_OutputFilenameFlag->GetFirstParameter()->GetValue());
//      write << all_descriptions;
//      io::File::CloseClearFStream( write);

      // end
      return 0;
    }

    //! @brief initializes the command object for that executable
    //! @return initalized command object
    util::ShPtr< command::Command> GenerateAlignmentDescription::InitializeCommand() const
    {
      // initialize a new command object
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // name of the file that has the list of alignment files
      sp_cmd->AddParameter( m_AlignmentListParam);

      // window size
      sp_cmd->AddFlag( m_WindowSizeFlag);

      // blast path
      sp_cmd->AddFlag( m_BlastPathFlag);

      // sspred path
      sp_cmd->AddFlag( m_SSPredPathFlag);

      // ssmethods
      sp_cmd->AddFlag( m_SSPredMethodsFlag);

      // flag to identify the ratio of positive pairs to be used
      sp_cmd->AddFlag( m_PositiveRatio);

      // flag for setting the output filename
      sp_cmd->AddFlag( m_OutputFilenameFlag);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief default constructor
    GenerateAlignmentDescription::GenerateAlignmentDescription() :
      m_AlignmentListParam
      (
        new command::Parameter
        (
          "alignment_list", "name of the file that has the list of alignment files", command::ParameterCheckFileExistence()
        )
      ),
      m_WindowSizeFlag( new command::FlagStatic( "window_size", "flag to change the size of the window to be used")),
      m_WindowSizeParam
      (
        new command::Parameter
        (
          "window_size_param", "size of the window", command::ParameterCheckRanged< size_t>( 0, 100), "5"
        )
      ),
      m_BlastPathFlag( new command::FlagStatic( "blast_path", "flag to give the path where blast profiles residue")),
      m_BlastPathParam( new command::Parameter( "blast_path_param", "path where blast profiles reside", ".")),
      m_SSPredPathFlag( new command::FlagStatic( "sspred_path", "flag to give the path where ss predictions reside")),
      m_SSPredPathParam( new command::Parameter( "sspred_path_param", "path where ss predictions reside", ".")),
      m_SSPredMethodsFlag
      (
        new command::FlagDynamic
        (
          "ss_methods", "flag for identifying the SSPred methods to use",
          command::Parameter
          (
            "ss_methods_param", "SSMethods of interest",
            command::ParameterCheckEnumerate< sspred::Methods>(),
            sspred::GetMethods().e_JUFO9D.GetName()
          ), 0, sspred::GetMethods().GetEnumCount()
        )
      ),
      m_PositiveRatio
      (
        new command::FlagStatic
        (
          "positive_ratio", "flag for setting the ratio of positive inputs to be use",
          command::Parameter
          (
            "positive_ratio_value", "ratio of positive to the total number of inputs",
            command::ParameterCheckRanged< double>( 0.0, 1.0),
            "0.5"
          )
        )
      ),
      m_OutputFilenameFlag
      (
        new command::FlagStatic
        (
          "output", "flag for setting the output filename",
          command::Parameter
          (
            "output_filename", "filename of the outputfile", "output.out"
          )
        )
      )
    {
      // attach parameters to flags
      m_WindowSizeFlag->PushBack( m_WindowSizeParam);
      m_BlastPathFlag->PushBack( m_BlastPathParam);
      m_SSPredPathFlag->PushBack( m_SSPredPathParam);
    }

    //! flag to define the number of requested entries for each of the states
    util::ShPtr< command::FlagInterface> m_NumberEntriesPerStateFlag;

    const ApplicationType GenerateAlignmentDescription::GenerateAlignmentDescription_Instance
    (
      GetAppGroups().AddAppToGroup( new GenerateAlignmentDescription(), GetAppGroups().e_InternalBiol)
    );

  } // namespace app
} // namespace bcl
