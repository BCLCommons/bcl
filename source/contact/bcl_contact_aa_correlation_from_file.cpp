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
#include "contact/bcl_contact_aa_correlation_from_file.h"

// includes from bcl - sorted alphabetically
#include "align/bcl_align_aligner_dp.h"
#include "align/bcl_align_alignment_leaf.h"
#include "align/bcl_align_handler_fasta.h"
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "function/bcl_function_binary_sum.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_z_score.h"
#include "score/bcl_score_aa_assignments.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AACorrelationFromFile::s_Instance
    (
      util::Enumerated< descriptor::Base< biol::AABase, float> >::AddInstance
      (
        new AACorrelationFromFile()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AACorrelationFromFile::AACorrelationFromFile( const std::string &SUFFIX) :
      m_FileSuffix( SUFFIX),
      m_FileExtension( ".corr_mat_bcl")
    {
      function::BinarySum< const biol::AABase, const biol::AABase, double> assign_score;
      AddScore( score::GetAAAssignments().e_PAM250, 2.5, assign_score);
      AddScore( score::GetAAAssignments().e_BLOSUM45, 2.5, assign_score);

      // Based on Elizabeth Dong's Weight tables -0.672 -0.907 -0.635  0.00
      // fold recognition scores ENCLOSED_SINGLE_GAP, ENCLOSED_MULTIPLE_GAP, BOUNDARY_SINGLE_GAP, BOUNDARY_MULTIPLE_GAP
      score::AssignmentWithGap< biol::AABase> assign_gap_score
      (
        util::CloneToShPtr( assign_score),
        -0.672, // ENCLOSED_SINGLE_GAP
        -0.907, // ENCLOSED_MULTIPLE_GAP
        -0.635, // BOUNDARY_SINGLE_GAP
        0.0 // BOUNDARY_MULTIPLE_GAP
      );
      m_Assignment = assign_gap_score;
    }

    //! @brief Clone function
    //! @return pointer to new BaseElement
    AACorrelationFromFile *AACorrelationFromFile::Clone() const
    {
      return new AACorrelationFromFile( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AACorrelationFromFile::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &AACorrelationFromFile::GetAlias() const
    {
      static const std::string s_name( "AACorrelationFromFile");
      return s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t AACorrelationFromFile::GetNormalSizeOfFeatures() const
    {
      return 1;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT_A, ELEMENT_B: the element pair of interest
    //! @param STORAGE storage for the descriptor
    void AACorrelationFromFile::Calculate
    (
      const iterate::Generic< const biol::AABase> &ELEMENT_A,
      const iterate::Generic< const biol::AABase> &ELEMENT_B,
      linal::VectorReference< float> &STORAGE
    )
    {
      // if the sequence maps are empty, this indicates that this object is now operating over a new sequence,
      // so it is necessary to reload the files
      if( m_SeqToMatrixColMap.IsEmpty())
      {
        LoadFiles();
      }

      storage::Map< util::SiPtr< const biol::AABase>, size_t, biol::AALessThanSeqID>::const_iterator
        itr_a( m_SeqToMatrixColMap.Find( *ELEMENT_A)), itr_b( m_SeqToMatrixColMap.Find( *ELEMENT_B));

      if( itr_a == m_SeqToMatrixColMap.End() || itr_b == m_SeqToMatrixColMap.End())
      {
        STORAGE( 0) = m_CorrelationMatrix.GetValue();
      }
      else
      {
        STORAGE( 0) = m_CorrelationMatrix( itr_a->second, itr_b->second);
      }
//      BCL_Message
//      (
//        util::Message::e_Standard,
//        util::Format()( ELEMENT_A->GetSeqID()) + " " + util::Format()( ELEMENT_B->GetSeqID()) + " "
//        + std::string( size_t( 1), ELEMENT_A->GetData()->GetType()->GetOneLetterCode())
//        + std::string( size_t( 1), ELEMENT_B->GetData()->GetType()->GetOneLetterCode())
//        + " " + util::Format()( STORAGE( 0))
//      );
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void AACorrelationFromFile::SetObjectHook()
    {
      m_SeqToMatrixColMap.Reset();
      m_CorrelationMatrix.Reset();
    }

    //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
    //! @return the cache preference, assuming this feature has its normal dimension setting
    descriptor::CachePreference AACorrelationFromFile::GetNormalCachePreference() const
    {
      return descriptor::e_PreferCache;
    }

    //! @brief function to load files; should only be called the first time Calculate is called with a new sequence
    //! since the results are often in the cache
    void AACorrelationFromFile::LoadFiles()
    {
      m_SeqToMatrixColMap.Reset();
      m_CorrelationMatrix.Reset();

      util::SiPtr< const assemble::ProteinModelWithCache> sp_protein_model( this->GetCurrentObject());

      // Get the filename
      util::ShPtr< util::Wrapper< std::string> > sp_filename_wrapper
      (
        sp_protein_model->GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile)
      );

      std::string pdb_filename( sp_filename_wrapper->GetData());

      // Remove the last extension
      std::string basename( io::File::RemoveLastExtension( io::File::RemoveCompressionExtension( pdb_filename)) + m_FileSuffix);
      // Add my extension
      std::string sym_matrix_filename( basename + m_FileExtension);

      // Read in the file's data
      io::IFStream input;
      io::File::MustOpenIFStream( input, sym_matrix_filename);
      // Read and close file stream
      io::Serialize::Read( m_CorrelationMatrix, input);
      io::File::CloseClearFStream( input);

      // Set Correlation Matrix default return value to average correlation
      m_CorrelationMatrix.SetValue( m_CorrelationMatrix.GetAverage());

      // Postfix string for files of precomputed alignments without correlation matrices of type mfasta
      const std::string msa_postfix_mfasta( ".mfasta");

      // Create alignment node
      util::ShPtr< align::AlignmentNode< biol::AABase> > alignment_node_ptr;

      BCL_MessageVrb( "Reading mfasta");
      // Create potential filenames
      std::string target_alignment_file_mfasta( basename + msa_postfix_mfasta);

      // Check for mfasta file first
      if( io::DirectoryEntry( target_alignment_file_mfasta).DoesExist())
      {
        // Create necessary stream and open it
        io::IFStream read_stream;
        io::File::MustOpenIFStream( read_stream, target_alignment_file_mfasta);

        // Create necessary PIR handler
        align::HandlerFasta< biol::AABase> mfasta_handler;
        alignment_node_ptr = mfasta_handler.ReadAlignment( read_stream, biol::AASequence());
        BCL_MessageVrb( "Read mfasta");
      }
      if( !alignment_node_ptr.IsDefined())
      {
        BCL_MessageStd( "Missing target alignment: " + target_alignment_file_mfasta);
        return;
      }

      // Get the AA sequence for the current protein
      util::ShPtr< align::AlignmentInterface< biol::AABase> > sp_seq_to_align
      (
        new align::AlignmentLeaf< biol::AABase>( sp_protein_model->GetChains().FirstElement()->GetSequence())
      );

      // Get the target sequence from the alignment
      util::ShPtr< align::AlignmentInterface< biol::AABase> > sp_alignment_target_seq
      (
        new align::AlignmentLeaf< biol::AABase>( alignment_node_ptr->GetSequences().FirstElement().HardCopy())
      );

      // Get rid of alignment as soon as first element is copied to save memory
      alignment_node_ptr = util::ShPtr< align::AlignmentNode< biol::AABase> >();

      // Set up aligner DP
      align::AlignerDP< biol::AABase> aligner;
      aligner.SetScoringFunction( m_Assignment);
      // Align the current protein's sequence to the target sequence
      BCL_MessageVrb( "Generate Alignment");
      const align::AlignmentNode< biol::AABase> align_with_target_seq
      (
        aligner.AlignPair( sp_seq_to_align, sp_alignment_target_seq).First()
      );
      BCL_MessageStd( "Aligned");
      BCL_MessageDbg
      (
        "score test is: " + util::Format()( aligner.AlignPair( sp_seq_to_align, sp_alignment_target_seq).Second())
      );
      BCL_MessageDbg
      (
        "Alignment size afterwards is: " + util::Format()( align_with_target_seq.GetSequences().FirstElement()->GetSize())
      );
      BCL_MessageDbg( "Generate Alignment done");

      // Iterate through the alignment and set the map values such that the current protein's amino acids point to the
      // sequence IDs of the aligned target alignment
      const util::ShPtrList< align::Assignment< biol::AABase> > &assignments( align_with_target_seq.GetAssignments());

      BCL_MessageDbg( "Assignments size:" + util::Format()( assignments.GetSize()));

      // Iterate through alignment and first check for gap in first (target) sequence
      for
      (
        util::ShPtrList< align::Assignment< biol::AABase> >::const_iterator
          itr( assignments.Begin()), itr_end( assignments.End());
        itr != itr_end;
        ++itr
      )
      {
        util::SiPtrList< const biol::AABase>::const_iterator itr_members( ( *itr)->GetMembers().Begin());

        // If not a gap iterate through second assignment through last checking for non-gap
        if( !itr_members->IsDefined())
        {
          continue;
        }
        // Store the current protein's AABase
        const biol::AABase &current_base( **itr_members);
        ++itr_members;

        for
        (
          util::SiPtrList< const biol::AABase>::const_iterator itr_assignment_end( ( *itr)->GetMembers().End());
          itr_members != itr_assignment_end;
          ++itr_members
        )
        {
          // pushback AA if non gap position found
          if( itr_members->IsDefined())
          {
            BCL_Assert
            (
              ( *itr_members)->GetSeqID() <= int( m_CorrelationMatrix.GetSize()),
              "Sequence ID is larger than given correlation matrix."
            );
            m_SeqToMatrixColMap[ current_base] = ( *itr_members)->GetSeqID() - 1;
//            BCL_MessageStd( "MP: " + util::Format()( ( *itr_members)->GetSeqID()));
            break;
          }
        }
      }
      BCL_MessageDbg( "Map size:" + util::Format()( m_SeqToMatrixColMap.GetSize()));
      BCL_MessageVrb( "Done with LoadFiles");
    }

    //! @brief adds the given pair score with the given weight to the scoring function
    //! @param SCORE the pair score to add
    //! @param SCORE_WEIGHT the weight for the pair score
    //! @param SCORE_FCT the scoring function (will be modified)
    void AACorrelationFromFile::AddScore
    (
      const score::AAAssignment &SCORE,
      double SCORE_WEIGHT,
      function::BinarySum< const biol::AABase, const biol::AABase, double> &SCORE_FCT
    )
    {
      util::ShPtr< math::ZScore> zscore( score::GetAAAssignments().GetZScore( SCORE));
      double offset( zscore->GetMean());
      double stddev( zscore->GetSigma());

      // correct weight to be a Z-score
      double Zscore_weight( SCORE_WEIGHT / stddev);

      SCORE_FCT += ( Zscore_weight * ( ( **SCORE) - offset));
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AACorrelationFromFile::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "data for amino acid pairs based on their sequence index and the given correlation matrix file with suffix " + m_FileExtension
      );
      parameters.AddInitializer
      (
        "suffix",
        "suffix of the file to be used for the correlation values (MSA must have the same basename_suffix) - default is just PDBFileBasename" + m_FileExtension,
        io::Serialization::GetAgent( &m_FileSuffix)
      );

      return parameters;

      return parameters;
    }

  } // namespace contact
} // namespace bcl
