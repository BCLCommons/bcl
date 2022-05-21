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

//////////////////////////////////////////////////////////////////////////////////
// automatically built code for simulating artificial neural network            //
//////////////////////////////////////////////////////////////////////////////////

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "contact/bcl_contact_ann.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector.h"

// external includes - sorted alphabetically
#include <cmath>

namespace bcl
{
  namespace contact
  {
    //! normalization values A*x+b
     const double ANN_NORMALIZE_CONTACT_HELIX_HELIX_A[ 544] =
    {
       0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.000833333, 0.000833333, 0.000833333, 1.2
    };

    //! normalization values a*x+B
     const double ANN_NORMALIZE_CONTACT_HELIX_HELIX_B[ 544] =
    {
       0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, -0.1
    };

    //! weights for layer 1
    const double ANN_WEIGHTS_CONTACT_HELIX_HELIX_LAYER_0[ 8704] =
    {
         0.092417,   0.063652,  -0.250940,   0.037502,  -0.143599,   0.203308,   0.137931,   0.157879,   0.239902,  -0.980965,  -0.527698,   0.107550,   0.298287,  -0.185739,   0.040926,   0.280935,   0.175741,   0.359968,  -0.756405,   0.260418,   0.285450,   0.361010,   0.625310,   0.358165,   0.080417,   0.068143,  -0.092924,   0.126960,  -0.662853,   0.075495,   0.157454,  -0.094829,
        -0.083609,  -0.013628,  -0.049112,  -0.011374,   0.044771,  -0.023007,   0.225046,  -0.551397,  -0.908100,   0.352592,   0.168542,   0.210194,  -0.006737,   0.206544,   0.077005,   0.106860,   0.108839,  -0.475650,   0.224929,  -0.001023,   0.186924,  -0.113542,  -0.412661,   0.023323,   0.095220,   0.185293,  -0.238625,   0.082033,   0.154339,  -0.238330,   0.131518,   0.053800,
         0.254620,   0.551053,   0.093016,  -0.191312,   0.111308,  -0.221015,  -0.326193,   0.120600,   0.111881,   0.197477,  -0.333557,  -0.192120,   0.020128,  -0.018659,  -0.041061,   0.096507,  -0.225787,  -0.128096,   0.308734,   0.128217,  -0.217549,   0.562777,  -0.023451,   0.060654,   0.410473,   0.187506,  -0.280354,  -0.092157,   0.105532,   0.207984,   0.015539,   0.007084,
         0.001245,   0.116281,   0.071678,   0.750332,  -0.655515,   0.155686,   0.051955,  -0.112436,   0.336343,   0.294829,   0.345635,   0.200539,  -0.183235,   0.274611,  -0.230887,  -0.192331,   0.325088,   0.102540,  -0.198242,  -0.351230,   0.219988,   0.181565,   0.209856,  -0.076948,  -0.228629,   0.271372,  -0.428046,   0.054287,   0.046348,   0.262583,   0.135738,   0.309954,
         0.064355,   0.716487,  -0.137707,   0.489550,   0.715721,  -0.583715,  -0.270295,   0.394844,  -0.025997,   0.080001,  -1.689140,  -0.853717,   0.258957,   0.594712,   0.909563,   0.006411,  -0.110465,   0.500837,  -0.157588,   0.139269,  -0.066579,   0.004603,   0.218267,   0.116667,   0.190966,   0.325156,  -0.015856,  -0.043438,   0.174067,   0.003430,  -0.041629,   0.517416,
        -0.070088,   0.158658,   0.212113,   0.381927,   0.184490,  -0.660205,   0.356214,  -0.015515,  -0.580944,  -0.415761,  -0.017866,   0.110834,   0.123887,   0.152970,   0.031879,  -0.322150,  -0.056466,   0.153489,   0.178223,  -0.403077,   0.070970,  -0.337699,  -0.198852,   0.113940,  -0.103076,   0.129729,  -0.148589,  -0.227241,  -0.117656,  -0.604884,   0.848101,  -0.274485,
         0.332908,  -0.168280,  -0.115805,  -0.626330,  -0.195461,  -0.135807,  -0.133067,   0.041390,  -0.734338,  -0.413934,  -0.450784,  -0.357415,  -0.104596,   0.081944,  -0.258033,   0.237241,  -0.310416,  -0.238567,  -0.734246,  -0.109813,   0.165001,   0.436753,  -0.082722,   0.286784,  -0.099372,  -0.055057,  -0.119545,   0.028327,   0.287157,  -0.149588,   0.055202,  -0.207417,
         0.811645,  -0.096836,   0.222582,   0.163670,  -0.407382,   0.400883,  -0.347715,   0.006527,   0.386766,   0.208074,   0.002256,   0.159907,  -0.196878,  -0.459771,   0.410498,   0.165473,  -0.436372,   0.406265,   0.180119,   0.255205,  -0.122568,   0.417992,  -0.150142,   0.192228,   0.157552,  -0.461412,  -0.274358,   0.196008,   0.293086,   0.375738,   0.331867,  -0.203781,
         0.408661,   0.058886,  -0.475376,  -0.147273,   0.099506,   0.296888,   0.586393,   0.117667,   0.117755,  -0.204708,   0.165455,  -0.167999,  -0.449077,  -0.311748,  -0.131598,  -0.114957,  -0.196305,  -0.015573,  -0.108986,   0.110729,   0.186779,  -0.135485,   0.111610,  -1.155700,  -0.877745,   0.461706,  -0.020110,  -0.095775,  -0.109713,   0.354545,  -0.057642,   0.240913,
         0.244263,   0.039324,   0.572714,   0.120702,   0.209203,   0.344993,   0.165568,   0.118755,  -0.274766,   0.094814,  -0.270605,  -0.288347,   0.242960,   0.051701,  -0.224632,  -0.095121,  -0.046578,   0.209557,   0.161124,   0.059490,  -0.040148,  -0.246578,  -0.638301,   0.293527,  -0.117897,  -0.322575,  -0.057775,  -0.074342,  -0.139661,   0.012304,  -0.330826,  -0.124563,
         0.101918,   0.452481,  -0.028132,  -0.038975,  -0.043866,   0.302553,   0.356304,  -0.049816,   0.357823,  -0.068480,   0.013432,   0.033825,   0.108763,   0.144499,  -0.243501,  -0.058275,   0.043212,   0.043718,   0.012883,   0.074736,  -0.256005,  -0.208094,  -0.000870,   0.196718,   0.250419,   0.070373,   0.085442,  -0.124056,   0.066919,  -0.032568,  -0.371626,   0.187030,
         0.068672,  -0.318821,   0.149139,   0.146611,   0.029381,   0.256072,  -0.302862,  -0.084388,   0.047899,   0.028452,  -0.018838,   0.079690,  -0.037355,   0.192702,   0.076212,  -0.090393,   0.154332,  -0.051788,  -0.403985,   0.082936,  -0.200634,   0.189129,   0.330158,  -0.065316,  -0.035525,  -0.359325,   0.080794,  -0.003606,  -0.277930,  -0.155861,   0.043089,   0.189123,
        -0.158608,   0.167768,  -0.302295,   0.456354,  -0.267559,  -0.310907,  -0.073708,   0.078525,   0.088445,   0.120493,  -0.063392,   0.270391,   0.359689,  -0.151253,   0.217030,  -0.282347,  -0.194921,   0.456527,   0.480011,  -0.040681,  -0.443584,   0.029392,   0.269489,   0.096359,  -2.480370,  -0.122358,   0.251067,   0.739252,   0.415947,   0.343461,  -0.205110,   0.207688,
        -0.132621,  -0.147735,  -0.248226,  -0.654831,   0.244156,   0.081923,   0.176266,   0.305907,   0.017287,  -0.213875,   0.194037,   0.028173,  -0.066390,  -0.311526,   0.280736,   0.376291,   0.007522,   0.041198,   0.318554,   0.090351,   0.174615,   0.635569,  -0.761956,  -0.154442,   0.277059,   0.391588,   0.072300,   0.080482,   0.500284,   0.019738,   0.212167,   0.320832,
         0.324815,  -0.027932,   0.176900,  -0.160182,  -0.494667,  -0.176839,  -0.117452,   0.066858,  -0.135982,  -0.237341,  -0.173742,  -0.013415,   0.594733,   0.157391,   0.139552,   0.344218,   0.025273,  -0.231792,  -0.011361,  -0.256434,   0.387133,   0.252646,  -0.813497,  -0.276147,  -0.158369,   0.110125,  -0.509233,   0.154137,  -0.064656,   0.146461,  -0.097948,  -0.257161,
        -0.575791,   0.116774,  -0.018104,   0.119864,   0.056982,   0.159422,  -0.016923,  -0.043453,  -0.313769,   0.417209,   0.562848,   0.146610,   0.167783,   0.497111,   0.304005,  -0.300601,  -0.065808,  -0.269683,  -0.280783,  -0.122642,   0.040012,   0.274650,   0.175732,  -0.070550,  -0.019210,   0.254412,  -0.236758,   0.338602,  -0.430897,  -0.357006,   0.098269,   0.138387,
         0.236403,   0.287032,  -0.081233,  -0.010970,   0.059889,  -0.094171,  -0.181962,  -0.311434,   0.689136,   0.333925,   0.540008,   0.160843,   1.009540,  -0.336731,  -0.050503,   0.705199,  -1.005180,   0.587382,   0.141151,   0.173147,   0.642069,  -0.036682,  -0.278765,   0.509165,  -0.324512,  -0.094163,   0.161713,  -0.030241,   0.206919,  -1.099180,   1.990020,  -0.092721,
         0.179962,   0.120846,  -0.103739,  -0.021018,   0.040262,   0.013016,   0.160182,   0.056110,   0.062530,  -0.583807,  -0.096073,   0.036638,  -0.193789,   0.004921,   0.141442,   0.250166,   0.037346,   0.047418,  -0.139376,  -0.001227,  -0.030934,   0.078680,  -0.149498,  -0.110602,   0.075405,  -0.300073,   0.031653,   0.117808,  -0.128325,  -0.268603,  -0.101099,   0.165863,
         0.012127,   0.097221,   0.127817,   0.145079,  -0.059090,   0.319890,  -0.039792,   0.155811,  -0.010886,   0.135983,  -0.158219,   0.010965,  -0.396450,   0.525783,  -0.150720,  -0.289030,  -0.109340,   0.013661,  -0.009700,  -0.013648,  -0.387428,   0.336615,   0.209139,  -0.420497,   0.156622,  -0.007757,  -0.104885,   0.052511,  -0.136848,   0.073256,   0.031843,   0.123636,
         0.282542,   0.233121,   0.226921,   0.078098,  -0.021277,  -0.133398,   0.346918,   0.197645,   0.338796,  -0.046678,  -0.164703,   0.030444,  -0.036756,   0.070512,   0.047246,  -0.144643,   0.108392,   0.036001,   0.028640,  -0.255933,  -0.382377,   0.284463,  -0.105645,   0.127686,  -0.413336,  -0.066336,   0.079506,   0.050764,   0.155593,   0.187835,   0.094281,   0.204212,
         0.147596,   0.095707,  -0.070765,  -0.044061,   0.875920,  -0.170469,   0.302588,   0.191597,   0.322038,  -0.393627,   0.006930,   0.241027,  -0.032330,  -0.117461,  -0.084865,   0.143207,   0.348504,  -0.200308,   0.087454,   0.144823,   0.094579,   0.080717,   0.219745,  -0.116307,  -0.141312,   0.193813,   0.113142,   0.297513,  -0.216482,   0.135949,   0.111025,   0.151773,
         0.087889,  -0.008850,   0.887854,  -0.246355,  -0.172400,   0.465683,  -0.047565,  -0.294657,   0.346146,  -0.046952,  -0.268471,   0.432654,  -0.182175,  -0.355027,  -0.180881,   0.066293,   0.286554,  -0.456277,   0.005068,   0.032305,   0.556650,  -0.146691,  -0.309970,   0.089276,  -0.187340,  -0.020460,   0.384285,   0.045055,   0.200617,   0.241240,  -0.117170,   0.476943,
         0.834025,   0.336243,  -0.463196,  -0.189161,  -0.321005,   0.510895,  -0.295698,  -0.221768,   0.180470,  -0.214261,   0.262572,   0.103987,  -0.481927,   0.144217,   0.303952,  -0.696920,   0.027094,   0.162849,  -0.519066,   0.048240,   0.151621,  -0.024475,  -0.034154,  -0.005266,   0.022406,   0.016796,   0.275194,  -0.113816,   0.013995,   0.611936,   0.177139,   0.133660,
         0.344538,   0.037754,   0.105355,  -0.256904,   0.233840,   0.226534,   0.393230,   0.036047,   0.030655,  -0.099965,   0.168050,   0.086734,  -0.061575,  -0.010106,   0.020647,  -0.106346,   0.035562,  -0.050805,  -0.113234,   0.097691,   0.269089,   0.297841,   0.071969,   0.069306,   0.067410,   0.026098,   0.095403,   0.402097,  -0.097500,  -0.336939,   0.392786,   0.499263,
         0.200082,  -0.354791,   0.490112,   0.649453,  -0.199823,   0.135101,  -0.107576,   0.062890,   0.171752,  -0.053751,   0.045775,  -0.248459,  -0.010312,  -0.028661,   0.109480,   0.016322,  -0.294413,  -0.006089,  -0.192833,   0.037185,   0.068400,   0.015855,  -0.039543,   0.054992,   0.109771,  -0.172997,   0.566183,   0.174223,  -0.558580,   0.381235,  -0.127281,   0.244395,
        -0.010306,   0.110804,   0.512040,   0.236999,  -0.129278,  -0.055181,  -0.667023,   0.331880,   0.071229,  -0.377358,   0.274369,   0.210571,   0.392584,   0.061745,  -0.204845,   0.117659,   0.140922,   0.185684,   0.047371,   0.037040,   0.014289,   0.136123,  -0.037438,   0.568823,   0.046431,   0.012486,  -0.245408,  -0.082350,  -0.025976,   0.395597,  -0.078602,  -0.012522,
        -0.019660,  -0.167878,  -0.095417,  -0.065003,   0.134025,  -0.100557,  -0.131603,  -0.465307,   0.026407,  -0.077889,  -0.180321,  -0.225019,  -0.180412,  -0.036065,   0.160208,   0.267834,   0.086740,   0.018944,   0.211896,  -0.027805,   0.085679,   0.016524,   0.119357,  -0.245767,   0.355789,   0.205518,   0.149309,   0.048555,   0.077073,   0.178808,  -0.071990,   0.122217,
        -0.084613,   0.042204,   0.467164,  -0.007442,  -0.038884,  -0.169672,   0.106198,  -0.116259,   0.054047,  -0.179949,  -0.119649,  -0.073265,   0.140343,   0.134489,   0.072504,   0.170151,   0.201581,  -0.104841,  -0.028180,   0.272161,   0.606851,  -0.242039,   0.026991,   0.123488,   0.402444,  -0.115909,   0.032272,   0.473850,  -0.017237,  -0.143729,  -0.041287,   0.011227,
        -0.133708,  -0.153758,  -0.107972,  -0.288416,  -0.041301,  -0.018174,  -0.102735,  -0.137427,  -0.257812,  -0.024980,   0.127536,  -0.024082,   0.246666,   0.187805,   0.104193,   0.058077,   0.017054,   0.482602,   0.597571,  -0.077427,  -0.350651,  -0.110867,  -0.166160,   0.457309,  -0.341864,  -0.370724,  -0.132865,  -0.183668,   0.064898,  -0.086403,  -0.446209,   0.088171,
         0.168514,  -0.593234,   0.082364,  -0.084436,   0.262088,  -0.156037,   0.188314,   0.015348,   0.096892,   0.301628,  -0.161378,   0.021607,  -0.144809,   0.065255,   0.011394,   0.474296,   0.628369,  -0.545293,  -0.291862,   0.275546,  -0.247938,   0.242630,   0.063560,  -0.239323,  -0.116325,   0.362516,  -0.098443,  -0.360136,   0.055437,  -0.053750,   0.186423,  -0.502647,
         0.248235,  -0.018719,   0.114275,  -0.079797,  -0.153261,   0.015383,  -0.081511,   0.085523,   0.064661,  -0.025358,   0.012542,   0.046266,  -0.005351,   0.303962,   0.446928,  -0.433507,   0.211334,   0.469305,   0.670604,  -0.370708,   0.215631,   0.458301,   0.285872,   0.268276,  -0.018288,  -0.155265,   0.365330,   0.077120,  -0.192074,   0.019844,   0.011191,  -0.135340,
        -0.021492,   0.033131,  -0.200278,   0.216281,   0.101146,   0.293545,   0.221442,   0.023688,   0.248315,   0.295259,   0.327916,   0.517761,  -0.348358,   0.098841,  -0.049721,   0.071487,  -0.046009,   0.238586,  -0.005285,   0.160617,  -0.385488,  -0.126940,   0.175932,   0.181922,  -0.455390,  -0.147028,  -0.048999,  -0.419007,  -0.131111,   0.317715,  -0.068470,   0.154461,
         0.147879,   0.066249,  -0.025944,   0.177840,   0.434635,   0.223550,  -0.070878,   0.225201,   0.199038,  -0.490785,   0.298623,  -0.094928,  -0.383241,   0.039361,  -0.111023,   0.104310,  -0.298084,  -0.258006,   0.007257,  -0.245042,   0.063364,   0.039089,  -0.345285,   0.095701,   0.337575,  -0.149863,  -0.077044,   0.184629,   1.183580,   0.220399,   0.157250,   0.092706,
         0.197211,   0.245791,   0.224096,   0.085519,   0.114405,   0.237111,   0.728451,  -0.408796,  -0.504201,  -0.208883,  -0.190530,  -0.132295,  -0.063748,  -0.070466,   0.019490,  -0.145069,   0.093648,  -0.437242,   0.317474,  -0.184581,   0.033076,   0.048032,  -0.014305,  -0.430588,  -0.194392,  -0.191132,   0.149081,  -0.224907,   0.175578,  -0.882268,   2.302120,  -3.000950,
         0.059163,   0.059564,  -0.137891,  -0.197507,   0.276141,   0.051400,   0.026116,   0.125860,  -0.483990,   0.527050,   0.286875,  -0.190441,  -0.251049,  -0.101355,  -0.191369,  -0.035114,  -0.297995,  -0.314594,  -0.118283,  -0.329596,   0.290016,   0.188665,  -0.608149,  -0.064815,   0.318414,  -0.006532,   0.008358,   0.042440,  -0.036409,   0.072922,   0.089428,  -0.013476,
        -0.186882,  -0.167315,   0.210361,   0.065923,   0.014474,   0.176882,  -0.054620,  -0.404756,   0.263770,  -0.131571,  -0.389140,  -0.187567,  -0.401367,   0.023237,  -0.178896,  -0.320896,  -0.156229,  -0.051780,   0.130816,   0.152954,  -0.346145,   0.147189,   0.012084,  -0.136733,   0.028240,   0.037808,   0.179752,  -0.181901,   0.375513,   0.169560,   0.160948,   0.040307,
        -0.170595,  -0.096321,   0.067464,   0.020303,   0.214155,  -0.449309,  -0.120909,  -0.087331,   0.132337,   0.144694,  -0.033307,  -0.028464,   0.182417,  -0.003709,   0.185565,   0.012526,   0.045248,  -0.239923,   0.160220,   0.020285,   0.162134,   0.198724,   0.331805,   0.320118,   0.227844,   0.151089,   0.144435,   0.079628,   0.111772,   0.138426,   0.107397,  -0.116730,
         0.128636,   0.135521,   0.309089,  -1.432480,  -0.161493,   0.096166,  -0.205776,   0.066633,  -0.214570,  -0.177940,   0.169076,   0.040554,   0.226668,  -0.106754,   0.214592,  -0.007877,  -0.432897,   0.112065,   0.038054,  -0.012313,  -0.092834,  -0.219912,  -0.036183,   0.133710,  -0.047271,  -0.069950,  -0.010492,  -0.167830,   0.346749,   0.014381,   0.161854,  -0.083291,
         0.290616,  -1.379590,  -1.117010,   0.108416,  -0.424705,  -0.755055,  -1.190700,   0.067200,  -0.482816,  -0.787635,  -0.136051,  -0.424385,   0.156971,   0.349535,  -0.876160,   0.067457,  -0.081153,  -0.606673,  -0.043424,   0.079169,  -0.001357,   0.050908,  -0.083523,  -0.081361,   0.056320,  -0.270400,   0.111708,  -0.041373,   0.083513,   0.006763,   0.288497,  -1.645560,
        -1.008450,   0.081175,  -0.032350,  -0.249177,  -0.477029,   0.251634,  -0.032812,  -0.365801,  -0.202759,  -0.102809,   0.050625,   0.185963,  -0.350181,   0.286776,   0.158314,  -0.756123,   0.118985,   0.312907,   0.014462,   0.157303,   0.099027,  -0.104013,   0.342287,  -0.072620,  -0.027736,   0.058424,  -0.023265,  -0.038350,   0.092987,  -0.795892,  -0.429002,   0.091548,
         0.235805,   0.183544,  -0.008576,   0.072329,   0.428101,   0.210383,   0.095201,   0.114041,   0.007226,   0.060437,   0.409883,   0.280612,   0.053979,  -1.005530,   0.168129,   0.256977,   0.060611,   0.118273,   0.089259,   0.070278,   0.128497,  -0.137139,   0.143344,  -0.112286,   0.264212,   0.163525,   0.071561,  -0.628181,   0.040982,   0.188589,  -0.445188,  -0.284041,
        -0.675272,   0.147822,  -0.126049,  -0.307534,  -0.207380,  -0.030037,   0.196925,   0.154345,  -0.405646,   0.078170,   0.220912,  -0.440275,   0.029292,   0.135919,  -0.122346,   0.278264,   0.324378,   0.055269,   0.019492,  -0.232122,   0.419211,   0.055343,   0.015825,   0.058420,  -0.142093,  -0.067221,   0.151128,  -0.145254,  -0.177154,  -0.219605,  -0.617931,   0.474680,
        -0.614448,  -0.447175,   0.108220,  -0.232452,   0.009592,   0.104242,  -0.463267,  -0.056075,  -0.068943,  -0.266148,  -0.092288,  -0.033487,  -0.311704,  -0.092684,   0.099546,   0.057794,   0.131668,  -0.054803,  -0.029272,   0.181868,   0.161964,  -0.003698,   0.157183,  -0.249839,  -0.437330,   0.055380,   0.365426,  -0.012450,   0.240724,   0.142683,   0.228819,   0.371926,
        -0.045925,   0.085745,   0.101886,   0.075041,   0.306552,  -0.041007,  -0.023374,  -0.100821,   0.286393,  -0.054856,   0.170795,   0.170407,   0.100664,   0.145582,   0.088783,   0.089467,   0.019097,  -0.090749,   0.226019,   0.168726,   0.104154,  -0.250802,  -0.599934,   0.063273,   0.050086,   0.195163,   0.345959,  -0.058406,   0.094364,   0.256976,  -0.178977,  -0.004723,
        -0.051037,  -0.208148,   0.303343,   0.033946,   0.043491,   0.126950,   0.129062,   0.162768,  -0.014687,   0.074072,   0.008037,   0.015911,   0.216995,   0.138228,   0.068017,   0.045841,   0.215432,   0.071460,   0.019977,  -0.629443,  -0.220951,   0.068315,  -0.021155,   0.067956,  -0.055639,  -0.113451,  -0.056789,  -0.085268,  -0.156824,  -0.114282,   0.101170,   0.108677,
         0.080631,   0.154027,   0.013617,  -0.113913,  -0.036853,  -0.025652,   0.062467,  -0.121815,  -0.039349,  -0.013924,   0.112861,  -0.013829,   0.046478,   0.049780,  -0.012209,  -0.074933,   0.170661,  -1.120470,   0.071538,  -0.248712,  -0.044169,   0.099786,  -0.075029,   0.123597,  -0.071476,   0.113869,  -0.049237,  -0.310202,  -0.021666,   0.129240,   0.096738,   0.123369,
        -0.161029,  -0.021735,   0.076307,  -0.173815,  -0.035468,   0.015278,  -0.021960,  -0.043029,   0.094122,  -0.096341,   0.132496,   0.175994,  -0.041192,   0.022541,   0.116249,  -0.912103,   0.454028,  -0.147989,   0.615167,   0.822510,   0.555590,  -0.068008,   0.363549,   0.533821,   0.122588,   0.477315,  -0.193787,  -0.364907,   0.831914,  -0.217605,  -0.270904,   0.055355,
         0.391914,   0.151092,   0.141408,   0.261171,  -0.268698,   0.079346,   0.135711,  -0.047844,   0.173579,  -0.018656,   0.005294,   0.123092,   0.052875,  -1.054590,   0.957983,  -0.107928,   0.009747,   0.152188,   0.551964,  -0.301389,   0.320131,   0.197557,  -0.033888,  -0.007895,  -0.195694,  -0.163175,   0.214862,   0.003178,  -0.093443,   0.731744,   0.184962,   0.074007,
        -0.077505,   0.008987,  -0.087261,   0.119197,   0.157729,   0.019951,  -0.004430,   0.107930,   0.163301,   0.125645,  -0.019008,  -1.031700,   0.977586,  -0.004034,  -0.008572,  -0.041200,   0.088318,   0.063942,   0.200301,  -0.009890,   0.009019,  -0.175993,   0.015730,   0.067787,   0.276107,  -0.095049,   0.078610,   0.089319,   0.019501,   0.013728,   0.399399,   0.022157,
         0.206585,   0.020961,   0.211403,   0.115034,  -0.038173,   0.045843,   0.040669,  -0.038774,  -0.226964,  -0.200345,   0.125053,  -0.483231,   0.021833,  -0.002465,   0.339412,  -0.190988,   0.151752,   0.390463,   0.075680,  -0.109582,  -0.027936,  -0.196798,   0.225106,  -0.037158,   0.001574,   0.084656,   0.199731,  -0.121369,   0.325432,  -0.007713,   0.022615,  -0.055645,
         0.162005,   0.117655,   0.012297,   0.120299,   0.137002,   0.097738,  -0.168689,   0.022841,   0.045241,  -0.080594,   0.314611,  -0.058061,   0.310773,  -0.112782,   0.248179,   0.149648,   0.126562,   0.079467,   0.028209,  -0.107544,   0.377049,   0.006623,   0.044603,  -0.020407,   0.163089,   0.115812,   0.084816,   0.294075,   0.066817,   0.958572,   0.794224,  -0.570295,
         0.035257,   0.026315,   0.137767,   0.064901,   0.027294,   0.138359,  -0.073164,   0.003994,  -0.150785,  -0.155851,  -0.359962,  -0.313452,   0.238368,  -0.048786,  -0.103406,   0.137942,   0.239508,   0.169139,   0.231653,   0.024152,  -0.105279,   0.087178,   0.138627,  -0.128265,  -0.234093,  -0.000604,  -0.185103,  -0.001747,   0.051837,   0.232247,   0.001085,   0.005520,
        -0.030657,  -0.108989,  -0.024843,   0.016928,  -0.117583,   0.147612,  -0.041227,  -0.221365,  -0.261796,  -0.088298,  -0.000320,   0.215166,   0.194163,   0.178340,   0.331806,  -0.123846,   0.035658,   0.004597,  -0.321446,  -0.277264,   0.013211,  -0.015548,  -0.025268,  -0.020287,   0.154706,   0.228188,   0.418498,   0.240156,  -0.246933,   0.044029,   0.159262,  -0.036658,
         0.005951,  -0.011719,  -0.121286,   0.032794,  -0.032603,  -0.051846,  -0.274372,  -0.228304,  -0.242871,   0.082543,  -0.242936,  -0.258528,   0.073025,  -0.273579,   0.148600,  -0.195088,  -0.014740,   0.132285,  -0.250592,   0.127368,   0.266668,  -0.215858,   0.125815,   0.018611,  -0.014379,   0.046782,   0.099953,  -0.108272,   0.187158,   0.037313,  -0.164152,  -0.182725,
         0.000017,  -0.026207,   0.049427,  -1.071170,   0.165974,  -0.093883,   0.034672,   0.389951,   0.013095,   0.022845,  -0.083906,   0.047217,   0.416825,  -0.269698,  -0.033412,   0.019773,   0.091086,  -0.087703,   0.125760,   0.045194,   0.110368,  -0.028105,  -0.283249,   0.079726,  -0.031570,  -0.015793,  -0.005043,   0.017806,   0.197642,  -0.108762,  -0.128377,   0.009568,
        -0.066829,  -0.829888,   0.521437,  -0.045779,   0.625651,   0.629959,   1.232710,  -0.068390,   0.271418,   0.402808,   0.239855,   0.323324,  -0.162254,  -0.105155,   0.378103,  -0.228400,   0.032444,   0.434100,   0.169890,   0.105700,   0.440569,   0.209239,  -0.094948,  -0.077367,   0.063679,  -0.179435,  -0.024463,  -0.015551,   0.017740,  -0.048773,   0.072176,  -1.017380,
        -0.032164,  -0.182558,  -0.115513,   0.121450,   0.309387,   0.064239,   0.062604,  -0.151035,   0.071220,   0.009493,  -0.006719,  -0.169796,   0.210291,  -0.050394,  -0.049776,   0.299892,   0.192011,   0.294106,   0.246128,  -0.142134,  -0.126080,  -0.106696,  -0.125535,  -0.101153,  -0.183842,  -0.116060,  -0.163980,   0.063824,   0.124357,  -0.489440,   0.187344,  -0.105736,
        -0.177746,   0.160630,  -0.007145,  -0.131497,  -0.012193,  -0.275928,   0.412554,  -0.130741,   0.068528,   0.004771,  -0.013467,   0.018836,  -0.050474,   0.070763,   0.093987,   0.069297,   0.053895,   0.161351,  -0.023032,   0.014130,  -0.006823,   0.024969,  -0.024947,  -0.223986,   0.155386,   0.205591,   0.122435,  -0.500683,  -0.125102,  -0.320949,   0.242059,   0.260354,
         0.216290,  -0.380922,   0.203534,   0.074603,   0.241675,  -0.133027,   0.150976,   0.006214,   0.537121,  -0.310626,   0.276958,  -0.089099,   0.285087,   0.064244,   0.242972,   0.146566,   0.251434,  -0.066323,   0.191250,   0.032990,   0.099938,   0.021903,   0.066765,  -0.022000,   0.103000,  -0.133629,  -0.315186,  -0.107188,   0.123356,   0.107557,   0.100517,  -0.031505,
         0.082651,   0.195402,  -0.171151,  -0.064845,  -0.005343,  -0.012804,   0.189412,  -0.055382,   0.238688,  -0.108441,   0.194923,   0.278655,   0.210497,   0.026886,   0.126048,   0.070599,   0.062724,   0.015452,   0.174489,   0.012906,  -0.004673,   0.042171,  -0.677556,   0.554260,   0.675956,  -0.199806,  -0.454487,   0.060018,   0.023369,   0.044574,  -0.233533,  -0.407873,
         0.171912,  -0.005968,   0.183328,  -0.047983,  -0.522675,  -0.144062,  -0.071507,  -0.086218,   0.258980,   0.284533,  -0.113488,   0.159965,   0.220837,   0.220258,   0.154441,  -0.163423,   0.210027,   0.203133,  -0.097239,   0.094251,  -0.131485,   0.181430,   0.328855,  -0.032378,  -0.136855,  -0.099557,  -0.203669,  -0.361445,   0.004924,  -0.395666,   0.181564,  -0.077965,
        -0.038867,  -0.019247,  -0.369714,  -0.019938,  -0.051021,   0.223318,   0.151165,   0.075033,   0.051458,  -0.044866,  -0.009445,   0.071017,  -0.066021,  -0.006425,  -0.130563,  -0.214169,  -0.121773,   0.126430,   0.158856,  -0.871490,  -0.136577,  -0.042884,  -0.142809,   0.008008,  -0.067791,  -0.032808,  -0.020077,  -0.033753,   0.226657,   0.182903,   0.067353,  -0.057329,
         0.004642,   0.196370,   0.400069,   0.597799,   0.124445,   0.169304,   0.427582,  -0.022426,  -0.001767,   0.111612,   0.223045,   0.039625,   0.121750,   0.050071,   0.275894,   0.036865,   0.530632,  -1.392660,  -0.492539,   0.122008,  -0.417238,  -0.193634,   0.136951,   0.087879,   0.006704,  -0.107535,   0.161273,  -0.228743,   0.428032,  -0.012565,  -0.339169,   0.133576,
         0.234429,  -0.235914,  -0.026823,  -0.058014,   0.255100,   0.234622,   0.195378,   0.164266,   0.165343,  -0.086869,   0.342308,   0.136878,   0.267657,   0.015029,   0.803169,  -1.578160,  -1.150500,   0.364686,  -0.198045,  -0.627899,  -0.736961,   0.477845,  -0.329917,  -0.109796,  -0.407437,  -0.280199,   0.254872,   0.466052,  -0.815608,   0.220226,   0.251389,  -0.500506,
         0.145415,   0.146412,   0.056083,   0.190866,  -0.127193,   0.099366,   0.354833,   0.057602,   0.020049,   0.171152,   0.340571,   0.085946,   0.977958,  -2.994820,  -1.357120,   0.363920,   0.007590,  -0.189549,  -0.217828,   0.472837,  -0.034090,  -0.332917,  -0.337993,  -0.353969,   0.282478,   0.391288,  -0.154306,   0.020791,   0.189180,  -0.631109,   0.184624,   0.222970,
         0.189362,   0.360679,   0.354899,  -0.151260,   0.360468,  -0.092346,  -0.054480,  -0.059518,   0.006092,   0.083664,   0.491263,  -1.250620,  -0.541895,  -0.040509,   0.106621,   0.293761,  -0.175942,   0.073906,   0.290674,   0.071647,   0.356739,  -0.044671,   0.193601,  -0.093766,   0.529918,   0.250936,   0.076303,  -1.029240,   0.234716,   0.073990,  -0.076610,   0.240937,
         0.292939,   0.122870,   0.143048,  -0.182350,   0.062093,   0.025782,   0.276345,   0.129782,   0.494652,  -0.828795,  -0.375061,   0.320468,  -0.328037,  -0.006179,  -0.360858,   0.343846,  -0.140263,  -0.298265,   0.070098,   0.005568,   0.170344,   0.139064,  -0.113026,   0.323099,  -0.017641,  -0.833286,   0.152324,  -0.104997,  -0.047716,   0.167148,   0.206509,  -0.068303,
        -0.022778,  -0.095282,  -0.000801,   0.139491,  -0.085532,   0.011939,   0.412081,  -0.817220,  -0.258817,  -0.072160,  -0.182645,  -0.066278,  -0.265532,   0.289044,  -0.209487,  -0.348051,   0.140139,  -0.082227,   0.083225,   0.205766,  -0.390782,   0.064601,   0.054667,  -0.020601,   0.054541,  -0.031587,   0.142790,  -0.126500,  -0.078826,   0.238108,   1.528620,   0.066713,
        -0.010094,   0.049222,   0.077971,   0.058346,   0.001659,  -0.123260,  -0.037058,   0.011380,  -0.053692,  -0.346102,  -0.006124,  -0.188363,   0.124673,   0.158878,   0.119773,  -0.086172,   0.328015,   0.301386,   0.101666,  -0.193399,  -0.103068,  -0.082502,  -0.044959,  -0.159946,  -0.110509,   0.264333,   0.049508,   0.075543,  -0.036982,   0.085180,  -0.027330,   0.048965,
         0.066548,  -0.006227,   0.020352,   0.033154,  -0.025462,  -0.064158,  -0.028280,  -0.758832,   0.003398,  -0.130917,  -0.055720,   0.174669,   0.261902,   0.057858,   0.273928,   0.063880,   0.049209,   0.022280,  -0.132517,  -0.331250,  -0.042541,  -0.017011,   0.094089,  -0.037553,   0.231251,   0.223028,   0.219605,   0.050528,  -0.086873,  -0.023529,   0.155456,  -0.023806,
         0.059197,  -0.054268,   0.135093,   0.127502,   0.066004,  -0.535117,  -0.070508,  -0.168403,  -0.117619,  -0.027928,   0.069043,  -0.061391,   0.164449,  -0.121422,   0.037513,  -0.218821,  -0.000986,   0.010968,  -0.114009,   0.235445,   0.054113,  -0.035213,   0.186314,   0.081782,   0.096744,   0.059488,   0.043157,   0.108240,   0.056198,   0.002203,   0.037094,  -0.141755,
         0.059517,   0.056949,   0.204156,  -1.175620,  -0.249296,  -0.082598,  -0.179919,   0.130121,   0.105628,  -0.144410,   0.118211,   0.006347,  -0.043475,  -0.208650,   0.010457,   0.015649,   0.078848,   0.003314,  -0.004816,   0.183591,  -0.022760,  -0.045244,  -0.078430,   0.055386,  -0.007763,   0.002839,  -0.034198,  -0.139945,   0.164710,  -0.065376,   0.035878,  -0.079709,
         0.148513,  -1.150210,  -0.334258,   0.014664,   0.938796,   0.443281,   0.865818,  -0.020886,   0.486624,   0.887584,   0.337225,   0.635629,  -0.089466,   0.006238,   0.592125,  -0.350821,  -0.100647,   0.465075,   0.380065,   0.069993,   0.192227,   0.249964,  -0.191394,  -0.077075,   0.047664,  -0.029002,   0.095272,  -0.149844,   0.170574,   0.032982,   0.326997,  -1.470780,
        -0.731765,  -0.021129,   0.026020,   0.287475,   0.258645,   0.192308,   0.051536,   0.199203,  -0.080509,   0.034129,  -0.042861,  -0.051560,   0.004321,   0.163586,   0.069773,  -0.060923,   0.211678,   0.227613,   0.017700,   0.095841,  -0.020408,  -0.100722,   0.153416,  -0.186340,  -0.067337,   0.068391,  -0.082430,   0.016478,   0.130726,  -0.762674,  -0.041746,  -0.145027,
        -0.070579,   0.204403,  -0.199022,  -0.003490,   0.131300,  -0.252898,   0.264469,  -0.170510,   0.204907,   0.101020,   0.161929,   0.126930,   0.151556,  -0.254701,   0.091172,   0.138835,  -0.040751,   0.129972,   0.209399,   0.112534,   0.113130,  -0.117234,  -0.082197,  -0.224762,   0.002369,   0.118276,   0.009716,  -0.703971,   0.101310,  -0.106271,   0.326397,   0.158420,
         0.101059,  -0.335666,   0.248032,   0.121438,   0.068636,  -0.162056,   0.059040,  -0.126242,   0.429878,  -0.051061,   0.090381,   0.019855,   0.151077,  -0.024667,   0.166403,  -0.070099,   0.076800,  -0.053885,   0.086961,  -0.103675,   0.090404,   0.076564,   0.069792,   0.069561,   0.126716,  -0.503097,  -0.076070,  -0.053714,   0.209037,   0.215119,   0.358024,   0.005729,
         0.102297,   0.513729,  -0.129091,   0.162103,  -0.119329,  -0.131789,   0.432068,  -0.179629,  -0.005556,   0.141983,   0.236996,   0.192992,  -0.007886,  -0.004716,   0.193090,  -0.038391,   0.087419,  -0.133412,   0.267425,   0.101475,  -0.125637,   0.103332,  -0.487990,   0.528387,   0.462448,  -0.086632,  -0.319526,  -0.237553,  -0.236291,  -0.066361,  -0.185721,  -0.472221,
         0.016371,  -0.146917,   0.144592,  -0.027848,  -0.559965,  -0.003412,   0.096715,  -0.238001,   0.046756,   0.025732,  -0.089483,   0.027640,   0.175747,   0.108258,  -0.046183,  -0.085799,   0.173364,   0.101939,  -0.020201,   0.088030,  -0.367229,   0.337463,   0.130776,   0.082504,  -0.322224,  -0.277205,  -0.470872,  -0.108541,  -0.025236,  -0.636994,   0.027023,  -0.193502,
         0.100448,  -0.040170,  -0.381366,   0.041710,   0.029236,  -0.069328,   0.008343,  -0.138091,  -0.012278,  -0.041279,   0.129179,   0.132411,  -0.023000,  -0.019220,  -0.002689,  -0.207265,  -0.100402,   0.003095,  -0.128620,  -0.541029,   0.322734,  -0.000797,  -0.093057,   0.120237,   0.004566,  -0.009722,   0.014297,   0.043727,   0.179588,   0.120974,   0.049448,  -0.294775,
         0.125191,   0.204098,   0.338045,   0.436090,   0.142024,   0.064204,   0.184223,  -0.067618,   0.054020,  -0.052633,   0.188216,  -0.024651,   0.012606,  -0.063484,   0.069961,   0.085852,   0.123928,  -0.812386,   0.051447,   0.102560,  -0.379925,  -0.169747,   0.112078,   0.168984,   0.111206,  -0.045316,   0.314305,  -0.102624,   0.340318,   0.044393,  -0.219056,   0.159487,
         0.086144,  -0.252659,   0.123226,   0.012158,   0.227001,   0.046599,   0.074107,   0.048107,   0.086088,  -0.177718,   0.541078,   0.180765,   0.043240,   0.028227,   0.342735,  -1.068970,  -0.546939,   0.164305,  -0.354029,  -0.423264,  -0.752643,   0.301195,  -0.509282,  -0.345354,  -0.165953,  -0.344049,   0.260896,   0.353350,  -0.726943,   0.201749,   0.191780,  -0.454125,
        -0.068479,   0.027800,   0.021640,   0.171061,   0.060898,   0.025897,   0.150921,  -0.033527,  -0.056649,  -0.043741,   0.247430,  -0.021601,   0.290765,  -1.182170,  -0.765214,   0.078672,  -0.190406,  -0.125919,  -0.163166,   0.364528,  -0.202013,  -0.330922,   0.010589,  -0.295807,   0.094598,   0.220565,  -0.280087,   0.112509,   0.092245,  -0.472785,   0.226054,   0.187879,
        -0.066920,   0.159227,   0.218353,  -0.025795,   0.179401,  -0.090720,  -0.038038,  -0.129537,  -0.000739,   0.139302,   0.110873,  -0.846282,  -0.145510,  -0.246714,  -0.067737,   0.431003,   0.129113,   0.112840,   0.321184,   0.103057,   0.154699,  -0.080800,   0.132258,  -0.036841,   0.440651,  -0.027075,   0.017285,  -0.727089,   0.377113,   0.301232,   0.115568,   0.023120,
         0.177975,   0.145558,   0.031840,  -0.049433,   0.186699,   0.007667,   0.165271,   0.001485,   0.254732,  -0.546279,  -0.344300,   0.127337,  -0.459056,  -0.168899,  -0.285203,   0.323127,  -0.224617,  -0.418927,   0.020761,  -0.089850,   0.158159,   0.114339,  -0.278173,   0.041148,   0.067172,  -0.528401,   0.211464,  -0.141565,   0.023078,   0.258170,   0.142324,   0.092329,
         0.059713,  -0.047664,   0.266121,   0.085729,  -0.136236,   0.042230,   0.053386,  -0.558915,  -0.245940,   0.004048,  -0.423644,  -0.128970,  -0.453396,   0.218384,  -0.214005,  -0.477566,   0.016718,  -0.152439,   0.150266,   0.192578,  -0.580744,   0.043527,   0.002738,  -0.131270,   0.038968,  -0.084142,   0.054643,   0.200861,   0.179510,   0.233577,   1.548830,   0.181734,
         0.016966,   0.054617,  -0.022384,   0.002164,   0.195524,   0.155793,   0.025243,  -0.059538,  -0.364068,   0.528985,  -0.007166,  -0.168143,  -0.135160,  -0.160730,  -0.384366,   0.188011,  -0.245551,  -0.354575,  -0.044301,  -0.115767,   0.240275,   0.084340,  -0.340741,  -0.082800,   0.324468,   0.124372,  -0.157503,  -0.054627,  -0.092327,   0.109002,   0.018851,   0.075413,
         0.102630,   0.022517,   0.023604,   0.232416,   0.026903,   0.101882,  -0.200546,   0.316232,  -0.133030,  -0.112819,  -0.099555,  -0.101321,  -0.241590,   0.250644,  -0.011094,  -0.207201,  -0.237897,   0.061495,  -0.079469,  -0.058478,  -0.104030,   0.043882,  -0.100317,   0.030878,   0.105364,   0.082757,   0.240179,  -0.034073,   0.119486,   0.025605,   0.104521,   0.054342,
        -0.089257,  -0.075794,  -0.051401,   0.154282,  -0.084642,   0.372179,  -0.219631,  -0.050310,   0.078161,   0.127912,   0.057853,   0.017093,   0.252252,   0.024718,   0.248311,  -0.051394,   0.009961,  -0.106411,   0.145637,   0.239293,   0.022927,   0.129013,   0.210083,   0.125074,   0.211296,  -0.039804,   0.200013,   0.051765,   0.169244,   0.038945,  -0.051742,  -0.089857,
         0.205482,  -0.094068,   0.155599,  -0.251020,  -0.502382,   0.065905,  -0.340083,  -0.237069,  -0.531052,   0.112543,  -0.116562,  -0.294117,  -0.153680,  -0.413120,   0.178994,   0.228684,  -0.386763,   0.244437,   0.293891,  -0.095281,  -0.288974,  -0.206531,  -0.031226,   0.107345,   0.029098,   0.020362,   0.288043,  -0.026673,   0.504981,  -0.051098,   0.177172,   0.053770,
         0.053221,  -0.319702,  -1.049090,   0.021280,   0.120419,  -0.246500,  -0.600259,   0.400896,  -0.247225,  -0.426341,  -0.155976,  -0.214305,   0.257942,   0.336035,  -0.449370,  -0.034065,   0.137088,  -0.825027,  -0.117449,   0.000123,   0.180311,   0.092747,   0.188079,  -0.069196,  -0.076085,  -0.177381,  -0.040731,  -0.032031,   0.088643,   0.028296,   0.022040,  -0.640052,
        -0.862209,  -0.020721,  -0.106858,   0.065019,  -0.425065,   0.209190,   0.058622,  -0.228063,   0.232323,  -0.101928,  -0.139884,   0.017627,  -0.030690,   0.140130,  -0.057408,  -0.496712,   0.222636,   0.044568,  -0.165926,  -0.093125,  -0.033179,  -0.077593,   0.076249,  -0.086232,  -0.102628,   0.050812,   0.037154,  -0.010514,  -0.067069,  -0.041182,  -0.056392,  -0.181689,
        -0.111828,   0.369141,  -0.243170,  -0.022908,   0.100944,  -0.318280,   0.169013,  -0.117339,   0.280876,   0.179011,  -0.035647,   0.378143,   0.045309,  -0.809154,   0.104225,   0.193919,   0.106913,  -0.015575,   0.166825,   0.102960,  -0.130654,  -0.182977,   0.158378,  -0.136272,  -0.029755,   0.027668,  -0.182005,  -0.209163,   0.324302,  -0.084029,  -0.177773,   0.051233,
        -0.333582,  -0.015661,  -0.169935,  -0.309601,   0.265468,   0.085611,  -0.179201,  -0.129953,  -0.148341,   0.000662,   0.187290,  -0.004495,   0.116026,   0.075124,   0.111313,   0.208722,  -0.029680,  -0.050787,  -0.113551,  -0.129491,   0.164082,   0.018499,   0.033885,   0.188592,  -0.154255,   0.083516,   0.278527,  -0.008101,   0.017083,  -0.124184,  -0.239902,   0.229375,
        -0.262629,  -0.004743,   0.077069,   0.098214,  -0.033417,  -0.039581,  -0.290358,   0.013351,   0.018097,   0.037290,   0.056622,   0.073340,   0.309844,   0.180760,   0.108196,  -0.043188,   0.008357,  -0.049731,   0.004365,   0.078436,   0.005329,  -0.013555,  -0.505842,   0.277601,   0.497942,  -0.226890,  -0.089848,   0.139015,   0.183348,  -0.101223,  -0.037871,   0.321943,
         0.080781,  -0.085431,   0.172758,  -0.043702,   0.053805,  -0.003884,  -0.190542,   0.156344,   0.318823,   0.248877,   0.268272,  -0.010373,   0.108523,   0.129310,   0.052310,  -0.016892,   0.153808,   0.092237,  -0.053444,   0.113376,  -0.327365,  -0.016573,   0.440423,  -0.196745,  -0.330085,   0.019896,  -0.073896,  -0.024285,  -0.089132,  -0.330399,   0.075174,  -0.051466,
        -0.027460,  -0.068141,  -0.394629,   0.062737,   0.130924,   0.158468,   0.091988,   0.167625,   0.090968,  -0.054305,   0.084779,   0.101766,  -0.016928,   0.040517,  -0.031040,  -0.029908,  -0.025014,  -0.026175,  -0.072439,  -0.618377,   0.467885,  -0.075484,  -0.231943,  -0.100524,  -0.035510,   0.281078,  -0.220444,  -0.240074,   0.169352,  -0.118323,  -0.107198,  -0.242074,
        -0.103877,   0.158014,   0.105918,   0.387594,   0.217262,   0.114943,   0.112264,  -0.078237,  -0.053030,  -0.137640,   0.200952,   0.007199,  -0.216833,   0.006775,   0.050210,  -0.024092,   0.239384,  -0.910612,  -0.022028,   0.009204,  -0.076005,   0.137272,   0.224448,  -0.000303,   0.172983,   0.277307,   0.415230,   0.047131,   0.039584,  -0.158678,   0.028411,   0.007183,
        -0.125916,   0.340647,   0.227052,  -0.094929,   0.065579,   0.001636,   0.053387,   0.053821,   0.029285,  -0.123696,   0.128562,  -0.131193,   0.130464,   0.148458,   0.242171,  -1.146900,  -0.714014,   0.250250,   0.196694,   0.359855,   0.464251,  -0.171165,   0.092912,   0.358505,  -0.080511,   0.147493,   0.089297,  -0.171412,   0.028496,   0.009952,  -0.138509,   0.164020,
         0.390510,   0.205805,   0.063617,   0.387656,  -0.042111,   0.202572,   0.198096,  -0.085382,   0.021703,  -0.170122,   0.201429,   0.074120,   0.484069,  -1.192130,  -0.945049,   0.099140,  -0.201405,   0.011228,   0.038281,   0.238834,  -0.043390,  -0.162580,  -0.304643,  -0.228105,   0.209062,   0.268129,  -0.228396,   0.206162,   0.159167,  -0.392359,   0.181356,   0.167538,
         0.523356,   0.328641,   0.167438,  -0.182939,   0.099964,  -0.250335,   0.042812,  -0.208525,   0.116297,   0.118962,   0.359359,  -1.131100,  -0.210485,  -0.077235,  -0.251503,   0.137852,  -0.135115,   0.361355,   0.192801,  -0.141072,   0.415478,  -0.231715,   0.102410,   0.026464,   0.058618,   0.141665,   0.109702,  -0.469768,   0.117373,  -0.127035,   0.293610,   0.047401,
         0.025111,  -0.033011,   0.162951,  -0.006351,   0.002836,   0.106993,   0.192554,   0.067662,   0.372628,  -0.891525,  -0.547270,  -0.061428,   0.021922,   0.224707,  -0.049708,   0.014645,   0.247587,   0.255086,   0.209069,  -0.097116,   0.144644,  -0.088736,   0.176307,   0.137705,  -0.233938,  -0.376736,   0.201573,  -0.046440,  -0.139509,  -0.067359,   0.138351,   0.042376,
         0.256799,   0.061145,   0.001575,   0.147116,   0.081503,   0.007919,   0.540712,  -0.805117,  -0.487145,   0.048372,   0.124048,  -0.025615,   0.332549,   0.164383,   0.382270,   0.322622,  -0.132471,   0.009209,  -0.183536,   0.101957,   0.280802,   0.139345,   0.073234,  -0.241349,   0.161119,   0.221829,   0.387062,   0.172071,  -0.001229,  -0.056395,   3.429960,  -1.300920,
         0.086571,   0.019875,  -0.154624,  -0.002354,   0.298267,   0.114839,  -0.130431,   0.273594,  -0.063658,  -0.575575,  -0.141951,   0.268340,  -0.446204,  -0.301693,  -0.043482,   0.350204,  -0.330122,  -0.225846,   0.444734,  -0.202832,   0.098849,   0.073394,  -0.296021,   0.123890,   0.237938,  -0.358943,  -0.072355,   0.103934,  -0.104372,  -0.020022,   0.357316,  -0.138788,
        -0.013793,   0.232307,  -0.180919,   0.087615,  -0.082205,  -0.107705,  -0.104435,   0.423199,  -0.058360,  -0.091185,  -0.016083,   0.436134,   0.194606,  -0.394851,   0.010603,  -0.072681,   0.126109,   0.180867,  -0.095882,  -0.115852,   0.179213,   0.092546,  -0.011771,  -0.450631,   0.119907,   0.163721,  -0.184578,   0.092364,  -0.129089,  -0.045201,  -0.217893,  -0.030123,
         0.384917,  -0.000194,   0.116775,   0.034210,   0.053433,   0.084586,   0.297144,   0.416938,  -0.092410,  -0.298265,   0.188054,  -0.013473,  -0.141440,   0.042124,  -0.058869,   0.040174,   0.220704,   0.321223,  -0.274637,  -0.146086,   0.029372,  -0.009141,  -0.143382,   0.192880,   0.039411,   0.271911,  -0.016412,   0.043235,  -0.054706,  -0.007686,   0.499184,  -0.038997,
         0.076499,   0.171866,   0.045760,   0.512079,   0.381336,   0.035049,  -0.588350,  -0.412770,  -0.384835,   0.251961,  -0.370339,  -0.189837,   0.066825,  -0.345336,   0.273613,   0.483996,  -0.421302,  -0.000982,   0.290997,  -0.303256,  -0.184490,  -0.372231,  -0.015609,   0.024512,   0.420238,   0.096279,  -0.070751,  -0.032372,  -0.225114,   0.253980,  -0.144522,   0.026604,
         0.024604,  -0.117875,   0.535173,  -0.067331,   0.197700,   0.044603,  -0.351193,  -0.029070,   0.133257,  -0.308525,   0.227595,   0.370318,   0.154526,  -0.022298,  -0.070596,   0.190960,   0.169219,  -0.368445,   0.061882,  -0.009651,   0.030923,  -0.401431,   0.267893,   0.006577,  -0.121526,   0.180956,   0.016255,  -0.083576,  -0.034679,   0.120145,  -0.016905,  -0.066710,
         0.495349,   0.026966,   0.235234,   0.393750,   0.744885,  -0.302432,   0.276070,   0.290612,   0.186383,   0.097009,  -0.043126,  -0.214373,   0.334979,  -0.036186,  -0.089809,   0.173179,  -0.056109,  -0.142152,   0.095004,  -0.029954,  -0.209766,   0.149968,  -0.176658,   0.049098,   0.335353,   0.031476,   0.161916,   0.259400,   0.219501,  -0.104602,   0.017652,   0.147966,
        -0.089036,  -0.354557,  -0.242829,   0.267251,  -0.344817,  -0.135834,  -0.304478,  -0.209225,   0.151559,   0.121060,  -0.117978,  -0.205179,   0.184042,   0.277473,  -0.094737,  -0.109643,  -0.355469,   0.304714,   0.131146,   0.116052,  -0.158327,   0.125607,   0.258220,   0.231029,  -0.002887,   0.187470,   0.299466,  -0.345582,  -0.269480,  -0.000003,  -0.099222,  -0.198794,
        -0.517629,  -0.134260,  -0.086355,  -0.353178,   0.129330,  -0.165426,   0.285986,   0.251698,  -0.178620,   0.094990,   0.386349,   0.200088,  -0.202370,  -0.113522,   0.418547,  -0.077517,   0.064136,  -0.105092,   0.018746,   0.136026,  -0.008194,   0.183856,   0.128909,   0.141874,   0.187428,   0.467464,  -0.428294,   0.151646,  -0.096445,  -0.022775,   0.115106,   0.102897,
         0.166006,   0.100736,  -0.060089,   0.033856,   0.142827,   0.086553,   0.151340,  -0.031912,   0.358727,  -0.398569,   0.098690,   0.127279,   0.002607,   0.148975,  -0.010636,   0.113765,  -0.360740,  -0.076405,   0.340988,   0.003335,   0.031248,   0.122435,   0.024160,  -0.421999,  -0.691113,   0.180090,  -0.151596,  -0.156763,  -0.485649,   0.096304,  -0.420166,  -0.715114,
         0.538412,  -0.043429,   0.076950,   0.059781,  -0.300126,   0.024496,   0.659217,   0.084781,  -0.052455,  -0.250346,   0.143580,   0.173760,  -0.110000,   0.039857,   0.012545,   0.174082,   0.048056,  -0.041012,  -0.186489,   0.051090,   0.218987,  -0.623381,  -0.415442,   0.025350,  -0.296831,   0.107215,   0.046228,   0.013278,  -0.003454,  -0.074917,   0.594453,   0.198200,
        -0.124721,  -0.172822,  -0.132631,   0.078265,  -0.007106,  -0.142230,  -0.132188,   0.274872,   0.178599,   0.219989,  -0.127235,   0.061547,  -0.073771,  -0.006420,   0.075079,  -0.158506,   0.091453,   0.009995,   0.066003,  -0.136287,   0.044804,   0.356975,  -0.109511,  -0.209997,   0.084537,   0.115260,  -0.168600,   0.031462,  -0.058853,   0.032541,   0.206582,   0.370454,
        -0.056795,  -0.074469,  -0.108862,  -0.075614,  -0.162986,   0.061889,   0.010627,   0.265467,   0.064154,  -0.091380,   0.123685,   0.015429,   0.364838,  -0.090349,   0.076424,  -0.009133,  -0.003982,  -0.323761,   0.432641,   0.016711,  -0.217731,  -0.218548,  -0.333580,  -0.080537,  -0.263523,  -0.381507,  -0.028022,   0.263181,   0.266586,   0.074620,  -0.482661,   0.001123,
         0.170864,  -0.498082,  -0.236604,  -0.162579,   0.210140,   0.428123,   0.060917,   0.144158,   0.148127,   0.240471,  -0.066044,  -0.009317,   0.005691,   0.164491,   0.235363,  -0.117842,  -0.003519,  -0.313288,  -0.064009,  -0.188066,   0.050509,  -0.123804,   0.073501,   0.183271,  -0.081089,   0.188574,   0.105544,   0.005551,   0.057143,   0.408172,   0.104685,   0.183296,
        -0.062895,   0.156555,  -0.080706,  -0.091439,   0.023822,   0.177850,  -0.074857,   0.065143,   0.143667,  -0.116837,  -0.180854,   0.084085,  -0.061841,   0.167721,   0.473510,   0.027975,  -0.234323,  -0.021481,   0.563694,   0.223834,  -0.032556,   0.392942,   0.289175,   0.023345,  -0.037042,  -0.026906,  -0.007294,   0.001989,   0.179957,   0.089699,   0.173037,   0.099791,
         0.007875,  -0.069227,   0.051197,  -0.068888,  -0.171310,   0.046724,   0.205734,   0.051952,   0.281786,  -0.101667,   0.125373,  -0.040559,   0.202032,   0.197098,   0.071947,  -0.362351,   0.177154,  -0.057842,  -0.117529,   0.077052,  -0.052446,   0.153543,  -0.034341,   0.189935,  -0.236570,  -0.131101,  -0.023747,  -0.190193,  -0.212963,  -0.255250,  -0.247547,   0.048321,
         0.015757,   0.244349,   0.060117,   0.187273,   0.259960,  -0.327430,   0.109090,   0.220732,  -0.023859,   0.437591,   0.216048,   0.007966,   0.072232,  -0.118349,   0.089033,   0.087108,  -0.195095,  -0.120204,  -0.194517,   0.041999,  -0.002682,   0.173138,  -0.270349,   0.156856,   0.220184,  -0.172590,  -0.116882,  -0.159540,   0.191706,   0.046825,   0.063002,  -0.029138,
        -0.049148,   0.018834,  -0.099411,   0.131862,  -0.065560,  -0.084923,  -0.130268,   0.598666,   0.557968,   0.027735,   0.196060,   0.370507,   0.090279,   0.102154,   0.216714,   0.288994,   0.622124,   0.405530,  -0.065625,  -0.041398,   0.444939,   0.109808,  -0.001367,  -0.189266,   0.059985,   0.110328,   0.203107,   0.110595,  -0.225783,   1.216410,   1.079800,   0.904600,
        -0.099044,   0.158188,  -0.207601,  -0.085350,   0.365184,   0.029480,   0.118419,   0.164167,  -0.526349,   0.954017,   0.495157,  -0.189778,  -0.308038,  -0.155582,  -0.224502,  -0.029038,  -0.161026,  -0.057619,  -0.184525,  -0.242927,   0.375115,   0.039273,  -0.521480,  -0.092157,   0.393752,   0.304545,   0.099958,   0.060405,  -0.105000,   0.007695,   0.131806,   0.127288,
        -0.114712,  -0.145528,   0.074407,   0.151968,  -0.131453,   0.188734,  -0.192795,  -0.159572,   0.343629,  -0.182316,  -0.212164,  -0.200617,  -0.468326,  -0.198987,  -0.077260,  -0.306174,   0.146110,   0.142884,   0.094344,   0.182978,  -0.299317,   0.117224,   0.055000,   0.089779,  -0.089237,   0.094847,   0.093580,   0.038574,   0.320130,   0.039889,   0.351572,   0.190697,
        -0.064618,  -0.066157,  -0.120265,   0.083039,   0.021344,  -0.422793,   0.058168,  -0.089673,   0.199709,   0.247266,  -0.036755,  -0.058352,   0.311835,   0.014103,   0.165368,   0.021594,  -0.158624,  -0.267151,   0.258854,   0.093660,   0.233817,  -0.055330,   0.373575,   0.190537,   0.031165,   0.083068,   0.112329,   0.152579,  -0.059453,  -0.098788,   0.091813,  -0.071036,
         0.068381,   0.111064,   0.260017,  -1.144330,  -0.275918,   0.276784,  -0.419855,  -0.243827,  -0.359129,  -0.077423,   0.163297,  -0.055496,  -0.121265,  -0.308991,   0.266204,  -0.061636,  -0.343923,   0.119590,   0.091765,   0.226989,  -0.131489,  -0.141489,   0.117579,  -0.002311,   0.123602,  -0.065290,  -0.059434,  -0.195621,   0.590225,   0.033956,   0.327562,  -0.016585,
         0.447766,  -1.213940,  -1.058330,   0.093796,  -0.335441,  -0.634408,  -0.932216,   0.000686,  -0.414634,  -0.648319,  -0.267751,  -0.216252,   0.420671,   0.564813,  -0.825480,   0.223329,  -0.154946,  -0.654451,   0.016516,  -0.016400,   0.104504,   0.193553,   0.159499,  -0.218816,   0.032709,  -0.198834,  -0.015755,   0.006010,   0.291231,  -0.031001,   0.525212,  -1.423050,
        -1.092800,   0.256233,   0.064446,  -0.040150,  -0.536540,   0.296878,   0.010225,  -0.050445,  -0.110089,  -0.087779,   0.113342,   0.331995,  -0.094746,   0.306067,   0.201896,  -0.980385,   0.232968,   0.196158,  -0.200673,   0.119678,   0.193578,  -0.104507,   0.330593,  -0.113066,   0.078955,   0.106242,  -0.013870,   0.056856,   0.198234,  -0.647480,  -0.322625,   0.019548,
         0.161082,   0.283408,   0.010530,  -0.078770,   0.350207,   0.153447,   0.244990,   0.080200,  -0.150565,  -0.051572,   0.334252,   0.261447,   0.110622,  -0.949419,   0.214798,   0.259810,   0.163307,   0.000764,  -0.004982,  -0.038616,  -0.104984,  -0.285737,   0.024358,  -0.194696,   0.276923,   0.057424,  -0.028789,  -0.696401,   0.138166,   0.350639,  -0.501197,  -0.392070,
        -0.529480,   0.104572,  -0.237020,  -0.210679,  -0.020083,  -0.023063,   0.046103,  -0.068942,  -0.387266,   0.125236,   0.244214,  -0.479860,   0.110934,   0.051104,  -0.132062,   0.398866,   0.060065,  -0.008931,   0.080051,  -0.092763,   0.347397,  -0.049643,   0.038479,   0.007812,  -0.034040,  -0.391243,   0.233106,   0.008623,   0.027228,  -0.229499,  -0.326123,   0.352848,
        -0.519588,  -0.214987,  -0.002042,   0.131016,  -0.138328,   0.194189,  -0.281793,   0.015667,  -0.019211,  -0.319996,  -0.072579,  -0.122530,   0.177431,   0.171393,   0.063457,  -0.176982,   0.009941,  -0.119222,   0.100946,   0.085288,   0.010888,  -0.097377,  -0.258634,  -0.071980,   0.212184,   0.071682,   0.450509,  -0.107303,   0.058148,   0.006893,   0.010803,   0.351258,
         0.214356,   0.027234,   0.168501,  -0.110340,   0.525732,   0.009068,   0.022339,   0.055291,   0.142977,  -0.287428,   0.004908,   0.285255,  -0.007030,   0.141969,  -0.138841,   0.038823,   0.179134,  -0.033950,   0.071074,  -0.033830,  -0.262391,   0.354733,   0.090750,  -0.147357,   0.023776,   0.126844,   0.333891,   0.051529,   0.097749,   0.300674,   0.022639,   0.141903,
        -0.172843,  -0.393597,   0.414293,  -0.062690,   0.072778,  -0.082811,   0.119097,  -0.152905,  -0.190191,  -0.021014,  -0.110296,   0.071213,   0.036220,   0.075813,   0.139919,  -0.084117,   0.260316,   0.133950,  -0.243450,   0.697699,   0.440209,   0.057041,   0.097740,   0.108105,   0.075067,   0.127897,   0.043618,  -0.050797,   0.001925,  -0.040211,  -0.112048,  -0.104684,
         0.305655,   0.104572,  -0.049044,  -0.032997,   0.124408,  -0.140126,  -0.140053,  -0.329556,   0.169949,  -0.154560,   0.077546,  -0.015896,  -0.175963,   0.160967,   0.013978,  -0.078634,  -0.113269,   1.321110,   0.323086,  -0.142384,   0.152995,   0.095381,   0.102296,   0.067490,   0.079600,   0.251254,   0.109742,  -0.043870,  -0.035117,   0.026182,   0.173609,   0.123271,
        -0.036589,   0.252789,   0.029827,  -0.132394,   0.001400,  -0.068949,  -0.098258,  -0.076348,  -0.138732,  -0.158584,   0.162758,   0.180710,  -0.101777,   0.021557,  -0.166186,   1.889550,   0.658534,  -0.109885,   0.455559,   0.458824,   0.551741,  -0.026886,   0.156433,   0.197104,   0.339509,   0.293833,  -0.144688,  -0.286555,   0.572377,  -0.214328,  -0.232814,   0.409874,
         0.309159,  -0.071793,   0.008043,   0.112107,  -0.012594,   0.052058,   0.047458,   0.039978,   0.077358,  -0.126386,  -0.050946,   0.043997,  -0.129051,   1.856180,   0.969757,  -0.105122,  -0.256002,   0.143382,   0.433741,  -0.132678,   0.059368,  -0.044400,   0.089561,   0.009681,   0.027618,  -0.283372,  -0.146381,  -0.148449,   0.017804,   0.169109,   0.330445,   0.016173,
        -0.239474,  -0.167108,   0.031998,   0.175373,   0.024285,  -0.076385,  -0.094409,  -0.011648,   0.097475,   0.019301,  -0.115465,   1.519370,   1.114730,   0.105668,  -0.127131,   0.060338,   0.018298,  -0.073696,   0.036056,  -0.207383,  -0.066437,  -0.067276,  -0.029687,   0.012616,   0.124877,   0.031056,  -0.108517,   0.258701,   0.094590,   0.022835,  -0.441882,  -0.232774,
         0.267641,   0.063839,   0.100622,   0.100437,   0.004237,  -0.044061,  -0.080424,   0.094382,  -0.180539,   1.242940,   0.108465,  -0.390702,  -0.021495,   0.136590,   0.403392,  -0.409284,   0.027427,   0.298907,   0.258446,  -0.001258,  -0.096263,  -0.296579,   0.253514,  -0.110218,  -0.055233,   0.359922,   0.186580,   0.118533,  -0.126518,  -0.180360,   0.117423,  -0.017318,
         0.080697,  -0.048447,   0.065235,  -0.026436,   0.101252,  -0.054722,  -0.276338,   0.918105,  -0.060736,  -0.142655,   0.199860,   0.089748,   0.330203,  -0.017850,   0.343272,   0.275131,  -0.004146,  -0.007610,  -0.134477,  -0.183917,   0.198038,   0.040503,  -0.147550,   0.124066,   0.262228,   0.054286,   0.002294,   0.118382,  -0.074566,   0.702586,   1.575870,  -0.606881,
         0.358551,   0.214895,   0.036474,   0.117703,   0.008277,   0.006457,   0.065646,   0.171822,   0.022749,  -0.452205,  -1.148060,   0.157864,   0.141835,  -0.007714,   0.127907,   0.148849,   0.060680,   0.192059,  -0.565681,   0.224397,  -0.029756,   0.175063,   0.198360,   0.129722,   0.248365,  -0.250549,   0.025714,   0.098952,  -0.172054,   0.014157,   0.073585,   0.136700,
         0.071936,   0.092138,   0.177466,   0.048988,   0.037105,   0.004838,   0.146324,  -0.148385,  -0.965788,   0.119999,   0.074089,   0.026103,   0.032222,  -0.009517,   0.097676,   0.027462,   0.177183,   0.113128,   0.014533,  -0.051461,   0.140851,   0.121633,   0.114765,  -0.028162,   0.042651,  -0.038793,  -0.107694,   0.045601,   0.057616,   0.187515,  -0.002790,   0.012394,
         0.330063,   0.047088,   0.035707,   0.197113,   0.051967,   0.297450,  -0.492446,   0.173488,   0.098288,  -0.229981,  -0.042075,   0.142927,  -0.004302,  -0.134328,  -0.119363,   0.038786,   0.386428,   0.287241,  -0.030552,   0.136618,   0.007911,   0.003429,  -0.135332,   0.076948,  -0.055891,   0.091570,   0.240580,   0.097639,   0.034353,   0.119666,   0.126494,   0.056875,
         0.122544,   0.045304,  -0.048847,   0.507815,   0.264330,   0.132737,   0.045946,   0.075647,  -0.060220,  -0.083097,  -0.062002,  -0.028334,   0.016509,  -0.602123,   0.224479,   0.269084,  -0.124582,   0.084203,   0.101523,  -0.198679,  -0.014740,   0.109132,  -0.368121,   0.009486,   0.236632,   0.028540,   0.063440,   0.108002,  -0.165762,   0.089194,   0.091161,  -0.070226,
        -0.049856,   0.161081,   0.141862,   0.011836,   0.485823,   0.382665,   0.472584,  -0.390561,   0.426290,   0.500546,   0.182798,   0.539826,  -0.140959,  -0.327997,   0.430027,  -0.047449,   0.154783,   0.016969,   0.248951,   0.028709,   0.118040,  -0.114423,   0.004231,  -0.030068,  -0.158993,  -0.028208,   0.062194,   0.041102,  -0.056235,   0.160858,  -0.062932,   0.479584,
         0.062849,   0.114787,  -0.006142,   0.118942,   0.258694,  -0.012489,  -0.010437,   0.046010,   0.297837,   0.134498,   0.044172,  -0.108507,  -0.071671,  -0.009724,  -0.064145,   0.084987,   0.022699,   0.003878,  -0.106917,   0.048371,   0.148039,   0.289905,   0.080317,   0.147421,   0.357059,  -0.002670,   0.142831,   0.071059,   0.012657,   0.254881,   0.161769,   0.173605,
        -0.002217,  -0.286517,   0.089685,  -0.004279,  -0.173988,  -0.062442,  -0.000842,   0.074290,   0.114785,   0.029950,  -0.184564,  -0.011024,   0.126748,   0.052845,  -0.085575,  -0.119457,  -0.633122,   0.034584,   0.213821,   0.177871,  -0.001139,   0.132474,   0.201801,   0.063828,   0.095928,   0.107417,   0.187148,  -0.081500,  -0.174096,   0.025941,  -0.033209,   0.120164,
         0.041994,  -0.056438,   0.009432,   0.139680,   0.020358,   0.141236,   0.175852,   0.133044,   0.089474,   0.064653,   0.097131,   0.154905,   0.053033,   0.111773,   0.120473,   0.099022,   0.172028,   0.110524,   0.050556,   0.152172,  -0.039672,  -0.036021,   0.155666,  -0.022810,   0.054760,   0.292161,  -0.107945,   0.170935,  -0.001141,   0.137090,   0.092016,  -0.025217,
         0.097300,   0.114680,   0.265061,   0.329210,   0.036904,   0.080309,   0.164126,   0.006642,   0.216045,  -0.980280,   0.183406,  -0.056918,   0.159191,   0.093895,  -0.042594,   0.158966,   0.022056,   0.043042,   0.182057,   0.109764,  -0.033495,   0.110781,   0.074337,  -0.192392,  -0.708976,   0.110636,  -0.501573,   0.098125,   0.137469,   0.155083,  -0.042993,  -0.127090,
         0.349861,   0.213995,  -0.067230,   0.046030,   0.043195,   0.024163,   0.286185,   0.046325,  -0.108461,  -0.025586,   0.031437,   0.027920,  -0.110996,   0.136816,  -0.057567,   0.099004,   0.249913,   0.019641,   0.044644,   0.243120,   0.136744,  -0.352462,  -0.586221,   0.126032,  -0.662962,  -0.472954,   0.104236,   0.187073,  -0.542790,  -0.485269,   0.240630,  -0.597544,
         0.193247,   0.114866,  -0.678816,   0.076787,   0.097193,  -0.103505,   0.091027,   0.228613,   0.036398,  -0.007492,   0.143201,   0.105653,  -0.042630,   0.023432,   0.285123,   0.116592,   0.139749,   0.194620,   0.045805,  -0.058634,  -0.575280,   0.004174,  -0.066540,  -0.196571,   0.041536,   0.221267,  -0.058388,  -0.027897,  -0.032406,  -0.015899,   0.232714,   0.307940,
        -0.128951,  -0.083794,   0.003324,  -0.273788,  -0.109021,  -0.017306,   0.106416,   0.205290,   0.020408,   0.074834,   0.015172,   0.025021,   0.145031,   0.085240,   0.085687,   0.077049,   0.032965,   0.132244,  -0.022740,   0.031810,  -0.055266,   0.134088,  -0.079338,   0.169556,  -0.316779,  -0.082772,   0.070780,   0.099747,   0.127482,   0.124158,  -0.143299,  -0.109683,
         0.128728,  -0.117435,   0.068548,   0.065596,  -0.188199,   0.037684,   0.017356,  -0.000014,   0.216178,   0.134996,  -0.253817,   0.104751,   0.013423,   0.146911,   0.107142,   0.047700,  -0.503372,  -0.009530,   0.231194,   0.269092,   0.484874,  -0.223634,   0.226802,   0.330617,   0.096257,   0.510186,   0.025135,  -0.450376,   0.312415,   0.076853,  -0.088209,   0.117830,
         0.315275,   0.226949,  -0.050853,  -0.104331,  -0.015647,   0.163054,   0.065790,   0.067027,   0.091457,   0.059163,   0.136693,   0.030056,   0.063218,  -0.387588,  -0.072325,   0.017261,  -0.001351,  -0.004720,   0.133099,   0.009861,  -0.068843,   0.133892,   0.227801,   0.234457,   0.146180,   0.083981,  -0.009077,  -0.035671,   0.235827,  -0.733965,   0.054543,  -0.039740,
         0.193335,  -0.013810,  -0.456008,   0.223468,   0.066974,   0.138761,   0.242696,   0.055890,   0.080934,   0.199210,   0.102975,   0.026793,  -0.134232,   0.149200,  -0.057498,  -0.282762,  -0.074349,   0.199216,  -0.156293,  -0.105168,  -0.056616,  -0.136191,   0.377037,   0.342654,  -0.472349,   0.161655,   0.117675,  -0.084568,  -0.075371,   0.055931,  -0.005818,   0.025259,
         0.229197,   0.076162,   0.047851,   0.179621,   0.017887,  -0.018451,   0.088999,   0.022501,   0.236884,  -0.057499,  -0.677926,  -0.010923,   0.263401,   0.069559,   0.295519,  -0.753089,   0.094462,   0.068591,   0.113082,  -0.014971,   0.015922,   0.141284,  -0.014634,   0.201561,  -0.066613,   0.185361,  -0.102315,   0.103815,  -0.034353,  -0.043824,   0.120788,  -0.012244,
         0.161404,   0.051389,  -0.071197,   0.083197,   0.007355,   0.158418,   0.189356,  -0.053797,  -0.916006,   0.110774,   0.310659,   0.036305,  -0.163938,  -0.058095,   0.213327,   0.255729,   0.190602,   0.209314,   0.045000,  -0.170152,   0.159327,  -0.086118,   0.075924,  -0.131180,   0.091572,   0.109706,   0.025745,   0.005014,  -0.016266,   0.207850,   0.994792,   0.186421,
        -0.034020,   0.080029,   0.010338,  -0.038055,   0.178560,  -0.062249,   0.047428,   0.054838,  -0.456993,   0.059268,   0.119285,  -0.023470,   0.092065,  -0.117021,  -0.076019,   0.010607,  -0.031393,   0.018675,  -0.088851,  -0.285678,   0.036967,   0.091802,  -0.185518,  -0.020808,   0.046241,   0.399855,   0.064194,  -0.000557,  -0.060882,   0.076703,   0.016860,  -0.009059,
        -0.039073,  -0.113099,   0.090938,   0.024486,   0.024983,  -0.010231,  -0.296248,  -0.335423,   0.118262,  -0.083768,  -0.321173,  -0.146258,  -0.162422,   0.066100,  -0.119178,  -0.172190,  -0.168834,  -0.301594,   0.088166,   0.135847,  -0.240804,   0.018302,   0.266272,  -0.035341,  -0.064285,  -0.024047,   0.095390,  -0.003284,   0.217847,   0.009519,   0.123492,   0.126460,
        -0.242227,  -0.108391,  -0.008535,  -0.142741,  -0.015945,   0.057771,   0.065099,  -0.138723,   0.062030,   0.180910,   0.201530,   0.010199,   0.363532,   0.026441,   0.037502,  -0.051733,  -0.104438,  -0.123311,   0.253596,   0.001399,   0.043700,   0.161508,   0.240303,   0.144745,   0.182273,   0.129318,  -0.060275,  -0.162056,   0.058128,  -0.094978,  -0.035134,  -0.009029,
        -0.117740,   0.002453,   0.080159,  -0.680515,  -0.082212,  -0.150508,   0.131617,   0.248409,   0.148291,  -0.238765,   0.210120,   0.102262,   0.203303,  -0.009988,  -0.070767,  -0.306354,   0.064349,   0.061615,  -0.106436,   0.272446,   0.044539,   0.001994,   0.051871,   0.069398,  -0.222394,   0.069810,   0.088111,  -0.216149,   0.546976,   0.012516,   0.123226,   0.097877,
         0.035615,  -0.446032,  -1.040700,   0.065457,   0.084588,  -0.579768,  -0.693409,   0.319828,  -0.321641,  -0.195343,  -0.216088,  -0.368885,   0.237785,   0.268994,  -0.395595,  -0.060375,   0.091231,  -0.001587,  -0.020175,   0.147742,  -0.069110,   0.128946,   0.301845,   0.091417,   0.057583,  -0.169207,  -0.067476,  -0.064396,   0.091530,   0.042599,   0.241929,  -0.818330,
        -1.099630,   0.073313,  -0.117864,  -0.141023,  -0.199838,   0.291513,  -0.203003,  -0.091161,  -0.094972,  -0.155376,   0.001216,   0.153960,  -0.161211,   0.025800,   0.184589,  -0.355248,   0.028724,   0.140828,   0.254800,   0.087229,   0.164092,  -0.066921,   0.169383,  -0.077379,  -0.108954,   0.004767,  -0.131257,  -0.073712,   0.008981,  -0.122201,  -0.379068,  -0.162260,
        -0.018919,   0.257348,   0.081367,  -0.154277,   0.315619,   0.218735,   0.160294,  -0.043537,  -0.216492,  -0.166917,   0.433796,  -0.085007,  -0.171193,  -0.382894,   0.093095,   0.038933,  -0.222198,  -0.124990,  -0.120009,  -0.040624,  -0.052172,  -0.097558,   0.087135,  -0.046800,   0.039410,   0.063750,  -0.094378,   0.026645,  -0.088959,   0.153853,  -0.147966,  -0.160158,
        -0.401891,   0.096981,  -0.029799,  -0.044467,  -0.032462,  -0.093513,   0.082410,  -0.006948,  -0.171924,   0.142829,   0.128709,  -0.181223,   0.124708,   0.042806,  -0.046840,  -0.004383,   0.275896,  -0.040595,  -0.000457,  -0.092905,   0.235880,   0.011926,  -0.043354,   0.109842,  -0.255283,  -0.114416,  -0.063954,   0.068720,   0.067992,  -0.206278,  -0.270711,   0.244547,
        -0.335602,  -0.150033,  -0.022159,  -0.214798,   0.033022,   0.089439,  -0.157831,   0.082711,   0.098287,   0.011740,  -0.076813,   0.064051,  -0.017568,   0.185407,   0.115465,  -0.063336,   0.037869,  -0.190499,   0.039294,   0.103808,   0.162885,   0.044296,  -0.633611,   0.004713,   0.602069,  -0.041612,  -0.089597,  -0.005084,   0.147723,   0.057441,  -0.009420,   0.013128,
         0.113581,   0.067981,  -0.038686,  -0.051415,  -0.062951,  -0.147960,  -0.178561,   0.087323,   0.243706,   0.049239,   0.219739,   0.061268,   0.106324,   0.020507,   0.102137,  -0.072655,   0.144796,  -0.030989,   0.056888,   0.085742,  -0.354665,   0.193216,   0.241062,   0.058035,  -0.238371,   0.018076,  -0.199019,   0.093139,  -0.198992,  -0.155467,   0.041918,   0.095143,
         0.045088,   0.004794,  -0.215317,   0.087611,   0.076862,  -0.337058,   0.040031,  -0.071358,   0.198112,   0.032116,  -0.044594,   0.062277,   0.085266,   0.092269,  -0.046046,  -0.037229,  -0.042012,  -0.003528,  -0.232723,  -0.094667,   0.272559,   0.072048,   0.092389,   0.178428,   0.289495,  -0.064662,   0.173716,   0.003826,   0.117050,   0.059543,  -0.071323,  -0.286716,
         0.238805,   0.096127,  -0.094409,   0.384844,   0.167662,   0.228357,  -0.143601,  -0.190403,  -0.083176,   0.020792,   0.086180,   0.101689,  -0.185139,   0.057157,   0.000319,   0.032580,   0.103010,  -0.215864,  -0.203915,   0.128909,   0.195289,   0.147166,   0.125977,  -0.092928,   0.344124,   0.288102,   0.139317,   0.230935,  -0.084161,  -0.239892,   0.417493,  -0.118925,
        -0.319819,   0.374505,   0.121863,  -0.029044,  -0.037723,   0.029976,   0.037031,  -0.016409,   0.001237,  -0.247625,   0.326578,  -0.052009,   0.104217,  -0.007620,  -0.023102,  -0.505972,  -0.488522,   0.247219,   0.080295,  -0.225309,  -0.120432,   0.060021,  -0.176337,   0.227686,  -0.083834,   0.141938,   0.144229,   0.083715,  -0.108766,   0.113777,   0.010979,  -0.270728,
         0.058016,   0.058736,   0.360098,   0.392460,   0.054488,   0.029986,   0.062365,  -0.117495,   0.144499,   0.074754,   0.082668,   0.170431,   0.112450,  -0.585743,  -0.551698,  -0.080768,  -0.160113,  -0.133337,  -0.153016,   0.229823,  -0.148896,  -0.419411,  -0.152349,  -0.213556,  -0.028992,   0.020007,  -0.342681,   0.122703,   0.018314,  -0.339040,  -0.027751,  -0.041328,
         0.072443,   0.089896,   0.025684,  -0.176077,   0.224991,  -0.050606,  -0.200977,  -0.064562,  -0.057023,  -0.032081,  -0.048105,  -0.447256,  -0.022756,  -0.103693,   0.003206,   0.065684,  -0.139329,  -0.078096,   0.205643,  -0.079154,   0.246142,  -0.212282,  -0.117458,  -0.170404,   0.376480,   0.132340,  -0.191677,  -0.234065,   0.221530,  -0.024141,  -0.205911,  -0.235301,
        -0.166931,  -0.073754,   0.125464,  -0.091180,  -0.122059,   0.014731,   0.052855,   0.023800,   0.076170,  -0.006849,  -0.242255,   0.065390,   0.054815,   0.153494,   0.169366,  -0.071746,   0.049651,   0.138000,   0.020098,   0.082632,   0.086464,  -0.251322,   0.152888,   0.034893,  -0.249315,  -0.243587,   0.230714,   0.019330,   0.065741,   0.024123,   0.125786,   0.245155,
         0.243341,  -0.051624,   0.076944,   0.031830,   0.053479,   0.207822,  -0.011294,  -0.044405,  -0.237685,   0.111120,   0.055853,  -0.048452,   0.022475,   0.070749,  -0.015191,   0.090353,  -0.174712,   0.098721,   0.014163,   0.046250,  -0.050682,   0.083501,   0.034582,  -0.175151,   0.132338,   0.010587,   0.174932,   0.284138,   0.036801,  -0.107956,   1.654620,  -0.731285,
         0.103319,   0.174733,   0.215336,   0.179731,   0.258160,   0.112686,  -0.025678,   0.227395,   0.299564,  -1.050800,  -0.224667,   0.068455,  -0.141952,  -0.080393,  -0.193158,   0.044082,   0.034822,  -0.113408,  -0.099736,   0.039113,  -0.001040,   0.105504,  -0.168431,   0.219125,   0.147179,  -0.380020,  -0.162601,   0.030005,  -0.037666,  -0.131262,   0.138225,   0.015576,
         0.002306,   0.083356,  -0.005239,   0.162644,  -0.031422,  -0.049277,   0.316514,  -0.033161,  -0.183917,   0.143818,  -0.181300,   0.067167,   0.037639,  -0.086170,   0.054763,   0.037234,  -0.105788,  -0.011952,  -0.189502,  -0.094975,   0.084330,   0.064850,   0.168760,  -0.426741,   0.020643,  -0.002893,  -0.460307,  -0.199778,  -0.050334,   0.148121,  -0.138927,   0.077173,
         0.341092,   0.041222,   0.281761,   0.166601,   0.106331,  -0.091499,   0.197288,   0.376572,  -0.052346,  -0.236750,  -0.234215,   0.066369,  -0.218372,  -0.094341,  -0.181857,  -0.060490,   0.327944,   0.501348,  -0.254543,   0.123104,   0.048832,  -0.149743,  -0.269413,   0.046004,  -0.275204,  -0.073815,   0.356921,   0.172365,  -0.063338,   0.003051,   0.260158,   0.043625,
         0.019268,   0.113949,  -0.072763,   0.455030,   0.323889,  -0.082615,  -0.213742,  -0.066899,  -0.146332,  -0.090473,  -0.274104,  -0.133414,   0.105728,  -0.151724,   0.180137,   0.303731,  -0.276483,  -0.137931,   0.117454,  -0.333870,  -0.081012,   0.036295,  -0.169353,  -0.154389,   0.183892,   0.076001,   0.198713,   0.207230,  -0.311530,   0.235371,  -0.105732,  -0.091606,
         0.084120,   0.303404,   0.112908,  -0.248802,   0.277802,   0.192499,  -0.083546,  -0.472224,   0.195761,   0.111080,   0.004009,   0.313772,  -0.254864,  -0.367823,   0.136622,  -0.099564,   0.000950,  -0.303838,  -0.040076,  -0.003875,  -0.076784,  -0.113102,  -0.020950,  -0.074941,  -0.264343,   0.037120,   0.355385,  -0.110778,   0.152374,   0.046577,   0.037060,   0.556323,
         0.449091,   0.373371,  -0.288487,  -0.095178,   0.048660,  -0.076838,  -0.137012,  -0.095598,   0.250691,  -0.194594,   0.004815,  -0.030136,  -0.097082,   0.035915,  -0.113552,  -0.199247,   0.029568,  -0.068079,  -0.305769,  -0.183141,   0.041967,   0.195844,  -0.027083,   0.132234,   0.336639,   0.086430,   0.128675,   0.088029,   0.231234,   0.297165,   0.244890,   0.194555,
         0.003783,  -0.319907,  -0.128312,  -0.086890,  -0.179889,  -0.104015,   0.042107,  -0.083311,   0.213589,   0.088607,  -0.232646,   0.020327,  -0.153292,  -0.127072,  -0.138466,  -0.149854,  -0.357037,  -0.114717,   0.140656,   0.143360,  -0.027034,   0.094772,   0.108162,   0.139634,   0.058936,   0.070883,   0.191729,   0.362051,  -0.002283,  -0.111913,   0.037543,   0.194788,
         0.096398,  -0.341387,   0.072010,   0.250312,   0.227962,  -0.076685,  -0.093118,  -0.063817,   0.040611,  -0.092715,   0.106701,   0.006479,  -0.045267,  -0.057081,   0.061989,  -0.115451,   0.043469,  -0.007068,  -0.087072,   0.130737,  -0.147811,   0.156804,   0.093201,   0.015594,   0.142393,   0.260701,   0.301713,   0.083228,   0.082933,   0.471837,   0.243254,  -0.146746,
         0.205242,   0.261108,   0.158353,   0.332182,  -0.380009,  -0.322457,   0.181418,  -0.008537,  -0.085990,  -0.302822,   0.130377,   0.036779,  -0.209267,  -0.116951,  -0.288451,   0.054164,   0.034004,  -0.015946,  -0.034312,   0.046072,  -0.008309,   0.059700,   0.020351,   0.139261,  -0.386459,  -0.053357,   0.018443,   0.090925,  -0.092757,  -0.127078,  -0.171287,  -0.137586,
         0.126654,   0.023456,  -0.084043,  -0.198881,   0.048631,  -0.060452,   0.180542,  -0.003178,  -0.048196,  -0.081151,  -0.210979,   0.150970,  -0.025106,   0.052831,   0.178219,   0.147049,   0.135516,   0.011602,   0.110136,   0.138231,   0.102589,   0.073671,   0.131916,  -0.113699,  -0.016443,   0.092189,   0.038230,  -0.061625,  -0.033882,   0.159165,   0.059826,   0.078599,
        -0.132383,  -0.052921,  -0.018208,   0.011987,   0.013340,  -0.157719,   0.001800,  -0.006117,  -0.013610,  -0.013340,  -0.114353,  -0.028966,   0.118124,   0.076429,   0.254363,  -0.005121,   0.051976,   0.017598,   0.041261,   0.462939,   0.326166,   0.271470,  -0.082775,  -0.301579,  -0.101447,  -0.150973,  -0.096083,   0.105294,  -0.005024,   0.061695,   0.046379,   0.207277,
        -0.022852,   0.014227,  -0.019166,  -0.207968,  -0.062330,  -0.059695,   0.061206,   0.206085,  -0.067256,   0.152647,  -0.034116,  -0.040511,   0.492257,   0.174641,   0.086109,   0.122567,   0.107862,   0.231049,   0.500728,   0.054373,  -0.208243,  -0.163241,  -0.455194,   0.043307,  -0.235832,  -0.376883,  -0.076853,   0.167359,   0.141400,   0.099166,  -0.272748,   0.083780,
         0.109962,  -0.546104,  -0.093217,   0.028009,  -0.084717,  -0.007699,   0.041094,   0.106427,   0.096881,   0.275758,  -0.175696,   0.156645,  -0.007590,  -0.020939,   0.037515,   0.397104,   0.258883,  -0.159019,   0.296139,   0.158282,   0.158645,  -0.334258,   0.349945,   0.263291,   0.122136,   0.336430,  -0.237503,  -0.412189,   0.170276,  -0.157922,  -0.065077,  -0.260602,
         0.048005,   0.099844,  -0.029570,  -0.358339,  -0.089697,   0.116045,   0.068512,   0.176331,   0.024652,   0.071291,   0.051737,   0.042374,   0.157977,   0.321609,   0.197118,  -0.006237,   0.166589,   0.010088,   0.349517,   0.005953,   0.004231,   0.119085,   0.134601,   0.302414,   0.131503,  -0.195880,   0.173875,  -0.128086,  -0.088223,  -0.090136,   0.116124,   0.049993,
         0.050688,  -0.124393,  -0.100136,   0.084958,  -0.042869,   0.106346,   0.298374,   0.177776,   0.283164,   0.081789,   0.253456,   0.252392,   0.196727,   0.323230,  -0.148406,  -0.344917,  -0.119372,  -0.019389,  -0.096664,  -0.094714,  -0.123561,  -0.041084,   0.090996,   0.266940,  -0.326259,   0.031883,   0.047709,  -0.321891,  -0.284967,  -0.298218,  -0.294325,   0.151165,
         0.052649,   0.224782,   0.073700,   0.052907,   0.107295,   0.080969,   0.019779,   0.151020,   0.239799,   0.302437,   0.349691,  -0.034296,  -0.001837,  -0.162062,   0.008305,   0.171197,  -0.186446,  -0.126066,  -0.216585,  -0.279517,   0.062617,   0.115354,  -0.154696,   0.029452,   0.111948,  -0.103649,  -0.229205,   0.056213,  -0.157033,  -0.097247,   0.142754,   0.026630,
         0.112833,   0.124735,   0.058555,   0.171692,   0.043412,   0.130132,   0.225876,   0.428797,   0.036800,  -0.178959,   0.124015,  -0.001255,   0.141865,  -0.125998,   0.056110,   0.143015,   0.228487,  -0.041491,   0.146261,  -0.175978,   0.295252,  -0.114406,  -0.084145,  -0.090062,  -0.013341,  -0.138424,  -0.114448,  -0.082599,  -0.086911,   0.392593,   2.273550,  -0.149059,
        -0.041614,  -0.052991,  -0.110093,   0.000584,   0.062681,  -0.234134,   0.075668,   0.024598,  -0.335670,   0.105606,   0.270632,  -0.058945,   0.048590,   0.289730,   0.289629,   0.077329,   0.182296,   0.393811,   0.027155,  -0.068302,   0.087976,  -0.157719,   0.001319,  -0.210225,  -0.170941,   0.449355,   0.056682,   0.061013,   0.067175,  -0.033294,   0.143437,   0.129751,
        -0.149236,  -0.120437,  -0.034533,   0.001701,  -0.033152,   0.173152,  -0.179290,  -0.538192,   0.214654,  -0.066065,   0.045408,   0.242377,   0.335759,   0.104626,   0.017174,   0.204464,   0.034933,   0.132259,   0.081700,  -0.182790,  -0.133018,  -0.202874,  -0.013417,   0.335521,   0.268185,   0.269252,   0.081379,   0.055971,   0.047058,   0.094487,   0.327255,   0.066084,
         0.005901,  -0.078843,  -0.004353,   0.039119,   0.022052,  -0.589332,  -0.003499,  -0.092429,  -0.246384,   0.070495,   0.148336,  -0.032288,   0.136217,   0.016532,  -0.104056,  -0.239686,  -0.068455,  -0.105927,  -0.069737,   0.103996,   0.144504,   0.110305,   0.130367,   0.012143,  -0.019872,   0.054475,   0.056201,   0.176423,   0.174034,   0.141209,   0.158666,  -0.250934,
         0.082111,   0.134633,   0.515356,  -1.212980,  -0.714540,   0.005670,  -0.173505,   0.031727,   0.278886,  -0.260416,   0.041005,   0.018004,  -0.214803,  -0.372210,   0.120825,   0.068687,  -0.082183,  -0.085850,  -0.056705,   0.422869,  -0.196257,  -0.054464,   0.089398,   0.058807,  -0.025103,  -0.010509,  -0.042372,  -0.048708,   0.003262,  -0.231299,   0.198139,  -0.117572,
         0.584464,  -1.335980,  -1.333480,   0.183627,   0.671379,   0.199099,   0.768003,   0.004313,   0.491839,   0.985037,  -0.047139,   0.477785,   0.090136,  -0.061394,   0.494240,  -0.136017,  -0.058049,   0.178572,   0.254978,   0.078431,   0.196565,   0.237301,  -0.045199,  -0.090689,   0.227266,  -0.048922,  -0.113632,  -0.360575,   0.275849,   0.047178,   0.539669,  -1.478920,
        -1.467300,   0.188522,   0.106667,   0.242537,   0.389785,   0.129924,   0.251303,   0.381522,  -0.546646,   0.159978,   0.044553,   0.215921,   0.100818,   0.321036,   0.070385,  -0.578909,   0.143886,   0.106409,   0.108202,   0.164927,  -0.022546,  -0.148087,   0.109737,  -0.006052,   0.007773,  -0.017827,   0.011027,   0.041956,   0.259704,  -0.754122,  -0.445129,   0.150370,
        -0.099491,   0.044712,  -0.290383,   0.135798,   0.237709,  -0.131496,   0.260991,   0.007444,  -0.018061,   0.207085,   0.274937,   0.315132,   0.236266,  -0.803563,   0.226521,  -0.042690,   0.249792,   0.090400,  -0.122022,  -0.184160,  -0.128215,  -0.107864,  -0.054523,  -0.197416,   0.138568,  -0.027941,   0.212024,  -0.798431,  -0.046982,  -0.193585,   0.173105,  -0.042105,
         0.140716,  -0.104763,   0.289785,   0.118773,   0.268209,  -0.036594,   0.031790,  -0.177315,   0.325750,  -0.058257,   0.026880,  -0.136865,   0.226716,  -0.042242,  -0.099945,  -0.096177,   0.095519,  -0.118889,   0.185784,  -0.022187,  -0.025921,   0.029375,   0.076310,   0.061414,  -0.052688,  -0.433540,   0.047420,   0.004234,   0.284120,   0.221176,   0.322146,  -0.025043,
         0.263885,   0.445997,  -0.120063,   0.433558,  -0.120938,  -0.139942,   0.254999,  -0.225227,  -0.115397,   0.191319,   0.223059,   0.176749,   0.190663,   0.151477,   0.124460,   0.084862,  -0.068719,  -0.033848,   0.275652,   0.142655,   0.024859,   0.019230,  -0.634756,   0.525061,   0.466343,  -0.001931,  -0.283642,  -0.306433,  -0.215212,  -0.036688,  -0.184839,  -0.254283,
         0.114501,  -0.388761,   0.219666,   0.169248,  -0.477758,  -0.108377,  -0.001748,  -0.271323,  -0.053311,  -0.212293,   0.105550,   0.131873,   0.061426,   0.000697,   0.021682,  -0.056971,   0.209847,   0.129576,   0.100918,   0.157230,  -0.291718,   0.564085,   0.119630,   0.135567,  -0.174913,  -0.168396,  -0.617570,   0.134816,  -0.062483,  -0.443410,  -0.050665,  -0.120684,
        -0.003960,  -0.087645,  -0.261674,   0.082427,   0.109511,  -0.313912,  -0.136689,  -0.094418,  -0.111228,  -0.102704,   0.173202,   0.074724,   0.125322,   0.139371,  -0.094920,  -0.191937,   0.033373,   0.083371,  -0.220740,   0.344579,   0.032551,   0.073602,   0.110604,   0.030605,   0.158484,   0.099679,   0.056174,   0.171597,   0.210622,   0.129653,   0.255922,  -0.070647,
         0.195076,   0.250738,   0.093163,   0.184374,   0.208948,   0.208160,   0.023506,  -0.057549,   0.087451,   0.023860,  -0.004007,   0.000420,   0.030723,   0.202635,   0.052689,  -0.014114,   0.073517,   0.333751,   0.008998,   0.126301,  -0.325731,  -0.071546,   0.033832,   0.205722,   0.155071,  -0.178078,   0.299238,  -0.141872,   0.313742,   0.075748,  -0.053328,   0.133118,
         0.046598,  -0.375603,   0.130804,  -0.025335,  -0.004861,  -0.178276,   0.038605,   0.133257,   0.132187,  -0.044149,   0.645202,   0.174227,   0.029956,   0.070578,  -0.041203,   0.235866,  -0.167475,   0.093205,  -0.376024,  -0.601585,  -0.887147,   0.324502,  -0.649940,  -0.550525,  -0.118441,  -0.178962,   0.307195,   0.236646,  -0.771596,  -0.094957,   0.192241,  -0.608476,
         0.042597,  -0.149285,   0.135015,   0.207041,   0.121218,   0.081520,   0.154575,  -0.035909,  -0.027282,   0.068914,   0.115965,   0.029542,  -0.173034,   0.294487,  -0.458537,  -0.088326,  -0.303249,   0.188485,  -0.289251,   0.215846,  -0.197246,  -0.296788,   0.409263,  -0.189214,  -0.207613,  -0.164641,  -0.253234,  -0.059948,  -0.018573,  -0.288258,   0.143080,   0.040729,
        -0.285918,  -0.023624,  -0.004747,   0.079605,   0.257920,   0.042844,   0.036003,  -0.146139,   0.134083,   0.195142,  -0.226966,   0.614988,   0.088292,  -0.181857,  -0.072468,   0.248871,  -0.017899,   0.061897,   0.216382,   0.058750,   0.256994,  -0.156165,   0.147954,  -0.128135,   0.338460,   0.030007,  -0.001688,  -0.388875,   0.334101,   0.261864,  -0.129404,  -0.049621,
         0.227770,  -0.031431,   0.120790,  -0.174599,   0.087042,  -0.064813,   0.090449,   0.109262,  -0.257606,   0.546705,  -0.064968,   0.124123,  -0.542036,  -0.165128,  -0.325362,   0.448514,  -0.308934,  -0.313902,  -0.053140,  -0.118700,   0.144407,  -0.171346,  -0.412354,  -0.053260,   0.080749,  -0.272528,   0.242372,  -0.068160,   0.062752,  -0.025489,   0.083097,   0.114566,
         0.175631,  -0.024959,   0.242242,  -0.005729,  -0.155826,   0.087195,  -0.188510,   0.389507,   0.078198,  -0.076535,  -0.372572,  -0.350636,  -0.499029,   0.273935,  -0.399204,  -0.368027,   0.213941,  -0.291234,   0.079093,   0.276046,  -0.551793,   0.170535,  -0.002604,  -0.119939,   0.014519,  -0.078509,  -0.117414,   0.019833,   0.306239,   0.225243,   1.476640,   0.236461,
        -0.093440,  -0.027020,   0.046736,  -0.144503,   0.074559,  -0.009799,   0.037062,   0.052576,  -0.583396,   0.885578,   0.478875,  -0.059469,  -0.050845,  -0.059681,   0.038836,   0.063901,  -0.074970,   0.017850,   0.036443,  -0.156761,   0.110621,  -0.027003,  -0.317470,  -0.072655,   0.179817,   0.303295,  -0.006630,   0.057814,  -0.183716,   0.164605,   0.084592,   0.078452,
        -0.067398,  -0.096301,   0.002923,   0.083462,  -0.101497,   0.070802,  -0.542151,   0.409462,   0.558408,  -0.164690,  -0.191319,  -0.122284,  -0.261241,   0.099570,  -0.113285,  -0.360920,  -0.053276,  -0.038465,   0.099451,  -0.056939,  -0.166143,  -0.032402,   0.054320,   0.211772,   0.065263,   0.149351,   0.078753,   0.014684,   0.310653,   0.071320,   0.072903,  -0.064992,
         0.038839,   0.003023,   0.003310,   0.066480,  -0.233208,   0.601287,   0.328234,  -0.057350,  -0.046220,   0.030708,   0.057019,   0.039668,   0.146218,  -0.124667,   0.212128,  -0.170067,   0.108857,  -0.055913,   0.147981,   0.162921,   0.086719,   0.393989,   0.269959,   0.138389,   0.280610,  -0.074538,   0.235895,   0.042636,   0.044266,  -0.021185,  -0.161266,  -0.007879,
        -0.029138,  -0.050979,  -0.024945,   0.060369,  -0.021829,   0.058603,  -0.303666,   0.099695,   0.100008,  -0.041201,   0.093754,  -0.050532,   0.428933,  -0.298846,  -0.018452,  -0.079290,  -0.082124,   0.047574,  -0.108917,   0.221464,   0.045939,  -0.054372,  -0.156557,   0.115089,  -0.167413,   0.074976,   0.031161,  -0.060828,   0.575887,  -0.064632,   0.062473,   0.138483,
         0.059848,  -0.020530,  -1.226340,  -0.031129,  -0.240353,  -0.444732,  -0.707020,   0.315116,  -0.428446,  -0.393924,  -0.320915,  -0.336578,   0.369121,   0.287003,  -0.423716,   0.051438,   0.071989,  -0.237389,  -0.082867,  -0.118979,   0.128306,   0.299254,   0.319624,  -0.019720,  -0.003084,  -0.118640,   0.044079,   0.015118,   0.211096,  -0.009703,   0.206842,  -0.318467,
        -1.217580,   0.168167,  -0.117332,  -0.133766,  -0.442172,   0.346165,  -0.073723,  -0.147058,  -0.103152,  -0.089623,   0.144436,   0.213332,  -0.294052,   0.127569,   0.243491,  -0.665970,   0.144566,   0.229364,   0.026339,   0.163915,   0.115064,  -0.128723,   0.149708,  -0.136703,  -0.076872,   0.047307,  -0.014260,  -0.109949,   0.107945,   0.015989,  -0.574360,  -0.157538,
        -0.091171,   0.346260,   0.109733,  -0.017893,   0.185358,  -0.115270,   0.176663,  -0.061381,   0.064994,   0.054950,  -0.096431,   0.330993,  -0.029072,  -0.630358,   0.294049,   0.124700,   0.036522,   0.005952,   0.058790,   0.040319,   0.045258,  -0.129807,   0.082986,   0.051472,   0.127652,   0.092358,  -0.033572,  -0.282005,  -0.154540,  -0.004693,  -0.295681,  -0.121784,
        -0.397648,   0.081765,  -0.145689,  -0.205611,  -0.053220,  -0.012478,   0.008661,  -0.005485,  -0.377301,   0.032621,  -0.027905,  -0.191328,   0.146329,   0.228246,  -0.048735,   0.190063,   0.115919,   0.051705,   0.042605,  -0.166713,   0.157132,  -0.118262,   0.023885,   0.108581,   0.050294,  -0.528731,  -0.203971,  -0.028990,   0.005862,  -0.271266,  -0.137674,   0.195585,
        -0.175954,  -0.015482,  -0.162554,   0.000716,  -0.095037,  -0.097976,  -0.253155,   0.076286,   0.073406,  -0.149900,  -0.043386,   0.010710,  -0.009033,   0.160832,   0.039273,   0.068954,   0.005444,  -0.061613,   0.016105,   0.005735,   0.051381,   0.009146,  -0.456011,  -0.051556,   0.613497,  -0.020852,   0.098272,   0.115740,   0.144220,   0.052367,  -0.013715,   0.174375,
         0.085297,   0.140258,  -0.000766,  -0.124705,  -0.076351,  -0.125337,  -0.113284,   0.275907,   0.192039,   0.084898,   0.070706,   0.133144,   0.145181,   0.253586,  -0.095150,  -0.119277,   0.159636,   0.082526,  -0.003582,   0.096642,  -0.381944,  -0.007076,   0.245791,  -0.050882,  -0.099373,   0.185212,   0.260599,  -0.238016,   0.019135,  -0.056649,   0.100524,   0.027772,
        -0.184284,  -0.238668,  -0.018035,  -0.157960,  -0.105378,   0.501179,   0.052998,   0.114294,   0.216081,  -0.010730,   0.029113,  -0.015505,  -0.027492,   0.074161,   0.041016,  -0.112007,   0.068519,   0.151681,  -0.173108,  -0.522732,   0.172703,   0.157997,   0.005118,  -0.090814,   0.023610,   0.099443,  -0.070948,  -0.113567,   0.082480,   0.050167,   0.039411,  -0.036723,
         0.047216,   0.080993,   0.162067,   0.293475,   0.210501,   0.149438,  -0.015048,   0.017623,  -0.006285,   0.093858,   0.076385,   0.014799,   0.111074,  -0.046913,   0.010309,   0.038584,   0.392600,  -0.738680,  -0.493607,   0.042032,  -0.124175,  -0.078077,  -0.252727,   0.165659,   0.128914,  -0.046671,  -0.011088,  -0.157338,   0.208044,  -0.004481,  -0.173958,   0.177353,
         0.047075,   0.100035,  -0.064059,  -0.173663,   0.126768,   0.085775,   0.131046,   0.032679,   0.167533,   0.047554,   0.095142,   0.152134,   0.117736,   0.033372,   0.431628,  -0.813889,  -1.245010,   0.193601,   0.533276,  -0.004551,   0.091397,  -0.154328,   0.139740,   0.348268,  -0.048815,   0.051010,  -0.025766,  -0.131289,   0.253994,   0.044746,   0.068763,  -0.200757,
         0.389939,   0.202301,   0.025323,   0.246418,  -0.004205,   0.073426,   0.393048,  -0.033567,  -0.026980,  -0.056813,   0.260164,  -0.023959,   0.477826,  -0.684468,  -1.213410,   0.187209,   0.178000,   0.065314,   0.256760,   0.171105,   0.177785,  -0.013769,  -0.406592,  -0.096218,  -0.042728,   0.115780,  -0.054782,   0.041481,   0.038703,  -0.380513,   0.342266,   0.163565,
         0.138099,   0.192333,   0.061012,  -0.246160,   0.043690,  -0.232251,  -0.022977,  -0.092015,   0.123000,  -0.015655,   0.150337,  -0.742820,  -0.274946,   0.041844,  -0.096975,   0.070917,  -0.293130,   0.169664,   0.106225,  -0.249413,   0.255971,   0.001371,   0.048324,   0.041671,   0.150956,   0.177495,   0.192562,  -0.726215,   0.120714,  -0.128634,   0.218598,   0.053540,
        -0.042604,  -0.093180,  -0.015042,  -0.032594,   0.035495,   0.017599,   0.038706,  -0.027070,   0.161608,  -0.510564,  -0.361477,  -0.106901,  -0.038449,   0.150700,  -0.289608,   0.179700,  -0.043576,   0.031904,   0.148529,   0.002158,   0.201594,   0.091974,   0.135958,   0.199156,  -0.049243,  -0.372515,   0.043117,   0.088744,  -0.329719,  -0.039320,   0.073019,   0.092586,
         0.069443,   0.020023,   0.011353,   0.021827,  -0.044435,  -0.011012,   0.011345,  -0.331389,   0.097630,   0.040988,   0.363832,   0.306148,   0.453730,  -0.173414,   0.321545,   0.295518,   0.348283,   0.130734,  -0.366868,  -0.132615,   0.399370,  -0.029597,  -0.303479,   0.111251,   0.229825,   0.257315,  -0.036242,   0.191501,  -0.046419,   0.307016,   1.631610,  -0.294803,
         0.032807,   0.018913,   0.067358,   0.001591,   0.072959,   0.217975,   0.174565,  -0.037326,   0.027449,  -0.035940,  -0.351548,  -0.068370,   0.435194,   0.183069,  -0.256097,  -0.118419,   0.147396,   0.169696,   0.071896,   0.245414,  -0.093816,  -0.020483,   0.288252,  -0.057912,  -0.041490,  -0.093549,  -0.082711,  -0.173139,   0.073955,   0.092411,  -0.059847,  -0.028310,
         0.011722,  -0.059427,  -0.020454,   0.169697,   0.078082,  -0.073527,  -0.036703,   0.461466,  -0.445591,  -0.101735,   0.007745,   0.172740,   0.162790,   0.030675,   0.112635,   0.122549,  -0.037363,   0.197224,  -0.099503,  -0.107394,   0.121213,  -0.126578,   0.190008,   0.014739,   0.004202,  -0.022095,  -0.170425,  -0.256580,   0.047482,  -0.036743,   0.085656,   0.045566,
        -0.137574,   0.053982,  -0.034988,   0.115781,  -0.193365,   1.143190,   0.015845,   0.088183,  -0.061443,   0.000662,  -0.005005,  -0.077127,   0.117334,  -0.147930,   0.023889,  -0.155575,  -0.240849,   0.071589,  -0.013112,   0.023768,  -0.325418,  -0.173140,   0.142841,   0.195111,  -0.543361,  -0.172817,   0.225454,   0.143270,   0.074729,   0.122932,  -0.020361,  -0.169888,
        -0.138533,   0.038773,  -0.139352,   1.199080,   0.476648,  -0.145476,   0.132145,   0.244411,  -0.139509,  -0.017886,   0.129519,  -0.041416,   0.141470,  -0.016869,  -0.173084,  -0.124136,   0.026785,   0.001469,  -0.100856,   0.195131,  -0.020110,   0.008050,  -0.277941,  -0.014020,  -0.019141,   0.045182,   0.222316,   0.039133,   0.270985,   0.195486,   0.122105,  -0.019437,
        -0.411413,   1.750970,   1.854080,  -0.090859,   0.901282,   0.460968,   0.248495,  -0.277916,   0.269093,   0.166863,   0.416189,   0.467645,  -0.274003,  -0.129460,   0.399535,  -0.183540,  -0.067491,   0.410833,  -0.170204,   0.188957,  -0.303664,   0.300960,  -0.193605,   0.031587,   0.099059,   0.078737,   0.148824,  -0.156377,  -0.020886,   0.118870,  -0.347658,   0.957959,
         1.616730,  -0.139559,  -0.062617,   0.120711,   0.465507,  -0.280939,  -0.088403,   0.116302,   0.107613,   0.274508,  -0.079782,  -0.123190,   0.089969,  -0.154886,  -0.170344,   0.539520,   0.069755,   0.066442,  -0.168191,  -0.201545,  -0.102400,   0.078650,   0.129052,   0.072770,  -0.056174,  -0.053720,  -0.112064,   0.097039,  -0.373752,   1.133850,   0.818158,  -0.022035,
         0.025157,   0.029040,   0.046208,  -0.197081,   0.007841,  -0.029677,  -0.111640,  -0.240029,  -0.150803,  -0.174162,   0.027039,  -0.018065,  -0.242194,   0.356193,   0.036343,   0.033413,  -0.310167,  -0.098002,   0.054698,   0.042594,  -0.043958,   0.090347,   0.070619,  -0.191040,  -0.045535,   0.138583,  -0.222872,   1.014710,   0.171265,  -0.318004,   0.216926,   0.160801,
         0.202368,  -0.317478,   0.070167,   0.204193,   0.028957,  -0.101901,  -0.154195,  -0.225781,   0.299636,  -0.312285,   0.047918,   0.373633,   0.022464,  -0.233369,  -0.361497,   0.009089,  -0.064478,  -0.133144,   0.220140,   0.041399,   0.023895,  -0.023730,   0.071392,   0.097675,  -0.082550,   0.628112,  -0.256013,  -0.019399,   0.155896,   0.072946,   0.232776,  -0.045546,
         0.306971,   0.204954,  -0.013016,   0.095442,  -0.093024,  -0.372390,   0.123282,  -0.093121,   0.064651,  -0.020130,   0.304611,   0.190384,  -0.167266,   0.051380,   0.178851,   0.066874,   0.058413,  -0.003569,   0.155682,   0.151543,   0.088911,   0.002204,  -0.603578,   0.186347,   0.647240,  -0.081573,  -0.371444,  -0.152311,  -0.136451,  -0.083223,  -0.190302,  -0.279782,
        -0.024955,  -0.141959,   0.127423,   0.243758,  -0.518865,  -0.017782,   0.047686,  -0.182472,   0.064526,   0.013558,   0.009289,   0.037373,   0.258461,   0.217929,   0.028284,  -0.151538,   0.199393,   0.113788,  -0.027239,   0.120444,  -0.198366,   0.163570,   0.209594,   0.200052,  -0.215755,  -0.210594,  -0.291669,  -0.386395,  -0.290013,  -0.535364,   0.037636,  -0.228466,
         0.338905,   0.109751,  -0.560074,   0.170209,   0.246810,  -0.328021,  -0.072802,   0.013713,   0.199409,  -0.135809,   0.359121,   0.153092,   0.314360,   0.333825,  -0.137605,  -0.043022,  -0.177153,   0.022620,   0.167359,  -0.395370,   0.001667,  -0.090131,   0.089344,   0.116014,   0.349029,  -0.329447,   0.164668,   0.005140,   0.085375,   0.403937,   0.103357,  -0.142237,
         0.123979,   0.152020,   0.362572,   0.171187,   0.009220,   0.146984,  -0.077288,   0.079665,  -0.103097,   0.173456,   0.250539,  -0.002427,   0.010829,  -0.105199,   0.159093,   0.213767,   0.415215,  -1.088050,  -0.509733,   0.187589,  -0.082363,   0.011182,   0.257624,  -0.296734,   0.217292,   0.096838,   0.130157,  -0.076953,   0.334771,  -0.038559,   0.034805,   0.104084,
         0.025557,   0.061811,  -0.089011,  -0.111488,   0.272195,   0.145428,   0.175937,  -0.063727,   0.219946,  -0.138129,   0.453396,   0.073453,   0.195060,   0.009732,   0.379832,  -1.313850,  -0.756266,   0.518079,  -0.492502,  -0.792262,  -0.850328,   0.124782,  -0.588638,  -0.283662,  -0.274975,  -0.143113,   0.348741,   0.562991,  -1.233860,   0.292372,   0.230542,  -0.760935,
        -0.054522,  -0.117907,   0.027780,   0.198729,   0.116033,  -0.044295,   0.098195,  -0.363894,  -0.020492,   0.037097,   0.303052,   0.079825,   0.453999,  -1.307370,  -0.749613,   0.305667,  -0.114102,  -0.081396,  -0.424574,   0.510745,  -0.333651,  -0.610939,  -0.123908,  -0.235357,   0.314507,   0.534694,  -0.371171,   0.131820,   0.190692,  -0.912981,   0.015517,   0.211581,
        -0.193251,   0.302262,   0.332317,   0.033279,   0.422891,   0.123154,   0.085449,  -0.029115,  -0.050172,   0.120003,  -0.014802,  -1.043550,   0.177834,  -0.248630,   0.322357,   0.204858,  -0.168728,  -0.235908,   0.146050,  -0.001766,   0.343513,  -0.142738,   0.258497,  -0.216987,   0.614166,   0.129741,  -0.039447,  -0.857292,   0.300635,   0.010983,  -0.031064,   0.022939,
         0.037878,   0.052275,   0.138710,  -0.133520,   0.254341,   0.040373,   0.225554,   0.103957,   0.005867,  -0.574319,  -0.231484,   0.556027,  -0.333056,  -0.031628,  -0.291653,   0.210829,  -0.146805,  -0.285853,  -0.078010,  -0.045313,   0.284450,   0.060667,  -0.080472,   0.144542,  -0.022544,  -0.235944,   0.039886,  -0.015493,  -0.018408,   0.223143,   0.203606,   0.016252,
        -0.017692,  -0.172028,   0.255943,   0.011849,  -0.019834,   0.131151,   0.000282,  -0.484098,  -0.004299,  -0.008076,  -0.236141,  -0.139594,  -0.160577,   0.190845,  -0.452452,  -0.382999,   0.022509,  -0.099489,   0.115994,   0.324369,  -0.605207,   0.114511,   0.074130,  -0.275228,  -0.020893,  -0.055875,  -0.007791,  -0.017751,   0.110415,   0.185989,   2.029480,   0.300512,
         0.116861,   0.173937,   0.111044,   0.337200,   0.076598,  -0.016640,   0.032220,   0.118944,  -0.018857,  -0.653259,   0.103582,  -0.041217,   0.023186,  -0.004739,   0.271216,  -0.133351,  -0.046889,   0.099052,  -0.040831,   0.180250,  -0.138533,  -0.116346,   0.132305,   0.018994,  -0.001552,  -0.294796,   0.016982,  -0.029795,  -0.051842,  -0.106558,  -0.083178,   0.117964,
         0.012707,   0.014311,   0.292926,  -0.019562,  -0.059979,   0.104215,  -0.090973,   0.344988,   0.258526,   0.225329,  -0.396925,  -0.135417,  -0.243706,   0.170721,  -0.174948,  -0.456179,  -0.107569,  -0.061964,   0.009686,   0.064245,  -0.478268,   0.191335,   0.150289,   0.011316,   0.001357,   0.094757,  -0.111687,   0.060846,  -0.018930,   0.011226,   0.053854,   0.050830,
         0.177111,   0.089702,   0.179250,   0.034551,   0.115421,   0.158595,   0.257395,   0.195188,   0.180387,  -0.376829,   0.016104,  -0.220359,  -0.142982,  -0.096341,   0.008745,  -0.007973,   0.155250,   0.283665,  -0.059587,  -0.023665,  -0.002804,   0.075733,  -0.226364,  -0.158822,   0.016597,   0.078289,  -0.082667,   0.156164,  -0.085142,  -0.056329,  -0.018345,   0.102546,
        -0.003945,   0.007152,  -0.063190,   0.477921,   0.234747,   0.111451,   0.198909,   0.143542,   0.392807,  -0.717515,   0.056751,   0.328687,   0.297930,   0.000252,   0.158106,   0.132048,   0.124340,  -0.058712,   0.013796,  -0.011050,   0.016029,   0.022291,  -0.283352,  -0.118241,   0.148206,   0.019664,   0.215252,   0.281767,  -0.357143,   0.106148,   0.101536,  -0.021761,
        -0.011396,   0.285002,  -0.030498,  -0.129759,   0.172686,   0.356056,   0.316070,  -0.149706,   0.243684,   0.149317,  -0.072632,   0.389283,  -0.263353,  -0.317585,   0.177641,  -0.066218,  -0.079912,  -0.234412,  -0.126933,   0.082778,  -0.255218,  -0.302449,  -0.032595,  -0.018749,  -0.196731,  -0.065406,   0.534120,   0.013414,   0.117397,   0.231744,  -0.007561,   0.382712,
         0.155809,   0.148266,  -0.120515,  -0.184942,  -0.118876,   0.027303,  -0.122294,  -0.349925,   0.044072,  -0.011247,   0.212050,   0.117777,  -0.270185,   0.244530,   0.070955,  -0.113344,  -0.111659,   0.012784,   0.118222,   0.195151,   0.072545,  -0.021161,   0.014163,   0.056646,   0.059219,   0.002329,   0.193664,   0.007538,   0.171051,   0.110154,   0.081343,   0.276350,
         0.240880,  -0.169249,   0.151532,   0.041066,   0.135988,   0.244205,   0.357683,   0.098087,  -0.199366,   0.152042,   0.136442,  -0.083141,  -0.398487,  -0.086844,  -0.122153,  -0.221973,  -0.014231,   0.085114,  -0.045447,   0.136708,   0.081431,   0.130919,  -0.079112,   0.248143,   0.023745,   0.066739,   0.204901,  -0.023401,  -0.291648,  -0.096799,   0.358778,   0.291768,
         0.429887,  -0.443037,   0.139778,   0.275252,   0.231521,  -0.036756,   0.213805,  -0.026568,   0.274901,   0.033489,   0.025584,  -0.332144,  -0.068694,   0.123163,  -0.150375,  -0.100351,  -0.055577,   0.021585,  -0.001856,   0.037858,   0.085316,  -0.023403,   0.167285,   0.076525,   0.232383,  -0.127457,  -0.135713,   0.302363,  -0.041969,   0.110752,  -0.092869,   0.117342,
        -0.068462,  -0.212004,  -0.040399,   0.270692,   0.004884,  -0.059936,  -0.195155,   0.088407,   0.086213,  -0.386439,  -0.058733,   0.053855,  -0.110919,  -0.098323,   0.099368,   0.054360,  -0.122811,   0.172119,  -0.001047,  -0.048123,   0.106072,   0.052775,   0.159481,   0.565576,  -0.994219,  -0.050903,  -0.292571,   0.008316,  -0.109221,   0.032830,  -0.177622,  -0.053197,
         0.096960,  -0.005183,   0.260358,  -0.053814,  -0.064843,   0.078634,   0.076555,  -0.124449,  -0.064468,  -0.100216,  -0.108082,  -0.169655,  -0.039107,   0.052613,   0.070260,   0.071145,   0.369252,  -0.145648,   0.137940,   0.044338,   0.159491,  -0.052141,  -0.282588,   0.054682,  -0.553250,   0.168071,  -0.013644,   0.019914,  -0.140388,  -0.084855,   0.151466,   0.087079,
         0.313514,   0.173057,  -0.493939,   0.202903,   0.321255,  -0.466237,  -0.030003,   0.177324,   0.104704,   0.130065,   0.119302,   0.065068,   0.068455,   0.014131,   0.413014,   0.339184,   0.100276,   0.180027,   0.048418,   0.387360,   0.144253,  -0.032548,   0.011704,  -0.308777,  -0.457341,   0.138986,  -0.177116,  -0.215215,   0.018600,  -0.183584,   0.339765,   0.302083,
        -0.418786,   0.079779,   0.179156,  -0.630156,  -0.075289,  -0.070731,  -0.175679,   0.106195,   0.119806,   0.060748,  -0.353517,  -0.174089,   0.146539,   0.026695,   0.043512,   0.045655,  -0.055651,   0.385245,   0.370806,   0.191021,  -0.349087,   0.155709,   0.165306,   0.188562,  -0.329346,  -0.126927,   0.228866,  -0.276310,  -0.042827,  -0.117589,  -0.408601,  -0.034482,
        -0.089503,  -0.079765,   0.258774,   0.047163,  -0.492985,  -0.406824,   0.040349,   0.154573,  -0.085015,   0.125557,  -0.000306,  -0.119522,   0.055668,   0.149160,   0.045307,   0.490731,  -0.086169,  -0.176241,  -0.601210,   0.301728,   0.243097,   0.230631,   0.064305,   0.108634,  -0.264373,   0.511150,  -0.001022,  -0.499995,  -0.311715,   0.036016,   0.052003,  -0.014443,
         0.209446,   0.119500,   0.044587,  -0.179759,   0.101997,   0.162490,   0.006549,   0.064972,   0.424088,  -0.209257,  -0.015434,   0.196826,   0.012382,   0.161281,   0.046036,   0.056952,  -0.371954,  -0.181681,  -0.054853,  -0.034154,  -0.309174,  -0.130919,   0.060839,  -0.012642,   0.395260,   0.193573,  -0.236977,   0.053437,   0.415983,  -0.117192,  -0.053131,  -0.095976,
         0.323403,   0.263100,   0.182000,   0.170842,  -0.049089,   0.223732,   0.320313,   0.005286,   0.387543,   0.105091,   0.265031,   0.419438,   0.007978,   0.368450,  -0.228083,  -0.221024,  -0.395369,   0.111994,  -0.108204,  -0.023201,   0.140429,  -0.097080,   0.101897,   0.202636,  -0.490712,   0.006238,   0.028280,  -0.638319,  -0.160888,  -0.037298,   0.030701,   0.144647,
         0.195315,  -0.117211,  -0.247473,   0.016824,   0.029595,   0.226999,  -0.043484,   0.092225,   0.207174,   0.017049,   0.258795,  -0.068026,  -0.020748,   0.239467,   0.187064,  -0.019789,   0.026080,   0.070000,   0.362539,  -0.014874,  -0.324617,  -0.160617,   0.147244,  -0.049239,   0.040443,  -0.057435,   0.098784,   0.081031,  -0.260288,  -0.472988,  -0.053992,   0.038793,
         0.023269,   0.114789,   0.072232,   0.174026,   0.011582,   0.124601,   0.204059,   0.286921,  -0.099019,  -0.007524,  -0.137512,  -0.252692,   0.006341,   0.154622,  -0.079212,  -0.145158,  -0.080439,  -0.108983,   0.138316,  -0.328297,   0.172532,  -0.135241,   0.017905,   0.056110,  -0.109626,  -0.073119,  -0.035084,  -0.013359,   0.141330,  -0.313972,  12.761300,  -0.292151,
         0.089405,   0.031290,  -0.032122,   0.128927,  -0.117850,   0.078253,  -0.057571,   0.074376,   0.042148,  -0.312357,   0.183582,   0.048744,   0.058351,  -0.189492,   0.223355,  -0.142385,   0.043312,   0.031676,   0.208084,   0.241129,  -0.125650,  -0.365448,   0.114485,  -0.041374,  -0.007074,  -0.098773,   0.111482,  -0.116154,  -0.146088,  -0.178479,  -0.069737,   0.047553,
        -0.236590,   0.094957,   0.101954,  -0.022120,  -0.005702,   0.173809,  -0.031178,   0.374158,   0.248412,   0.164729,  -0.259221,  -0.022815,   0.221466,  -0.005833,  -0.094237,  -0.046117,  -0.128020,   0.087827,  -0.153753,  -0.086682,  -0.179990,   0.058590,  -0.037630,  -0.076819,   0.017239,   0.117510,  -0.055502,  -0.059077,  -0.031716,   0.120119,  -0.219047,   0.051758,
         0.334973,  -0.048986,   0.122157,   0.189929,   0.038075,   0.155329,   0.458863,   0.272508,  -0.154017,  -0.393529,  -0.211280,   0.046470,  -0.272215,  -0.023233,  -0.110863,  -0.131935,   0.212176,   0.231101,  -0.266499,  -0.090951,  -0.142980,  -0.075355,  -0.279027,   0.061028,  -0.100465,  -0.012123,   0.249768,   0.033932,  -0.058571,  -0.038471,   0.167998,   0.143350,
         0.083944,   0.105175,  -0.009631,   0.261484,   0.696063,   0.034905,  -0.214132,  -0.125428,  -0.184438,   0.117858,  -0.244028,  -0.211956,   0.136039,  -0.065318,   0.032176,   0.077896,  -0.249021,   0.018314,  -0.023216,  -0.356528,   0.062530,  -0.024029,  -0.204402,  -0.160229,   0.085700,   0.021805,   0.300474,   0.321684,  -0.433234,   0.293832,   0.018789,  -0.052280,
         0.080961,   0.182751,   0.521289,  -0.167441,   0.558791,   0.665329,   0.821046,  -0.579023,   0.466464,   0.526477,   0.119497,   1.029860,  -0.373600,  -0.478296,   0.448799,  -0.101083,   0.000885,  -0.315308,  -0.049549,  -0.131903,   0.065740,  -0.080177,  -0.172701,  -0.029186,  -0.079454,   0.101013,   0.328929,  -0.164057,   0.109148,   0.054953,   0.083050,   0.223366,
         0.462968,   0.228296,  -0.065500,  -0.130967,   0.146257,  -0.213124,  -0.010116,  -0.140572,   0.157287,  -0.006139,   0.051396,   0.028087,  -0.163590,   0.049105,  -0.066907,   0.017680,  -0.160739,  -0.137708,  -0.332022,  -0.119335,  -0.029541,   0.091504,  -0.001352,   0.202222,   0.359272,   0.247943,   0.237583,   0.095221,   0.252024,   0.256312,   0.043012,   0.218136,
         0.069503,  -0.421780,  -0.271005,   0.192870,  -0.159199,  -0.031415,   0.058788,  -0.008956,   0.176274,   0.135051,  -0.195945,  -0.019644,   0.045548,  -0.376611,  -0.198818,  -0.184494,  -0.298103,   0.103266,   0.063377,   0.013125,  -0.150769,  -0.017453,   0.121159,   0.087055,   0.105806,   0.066394,   0.178769,   0.289701,   0.063275,   0.095438,  -0.082187,   0.177832,
         0.170008,  -0.368673,  -0.120643,   0.014077,   0.186644,   0.030952,   0.078499,  -0.080909,  -0.104990,  -0.051333,   0.122142,  -0.059738,   0.126194,   0.190689,  -0.281614,  -0.066073,   0.066363,   0.008405,   0.143228,   0.082249,  -0.106407,  -0.027613,   0.173894,   0.024195,  -0.009869,   0.290659,   0.205885,   0.172849,   0.085924,   0.394003,   0.102058,  -0.141058,
         0.096637,   0.070103,   0.278608,   0.289910,  -0.117283,  -0.203383,   0.109411,   0.021140,  -0.121858,  -0.281600,  -0.065015,  -0.001953,  -0.126286,  -0.130985,   0.026152,   0.044723,   0.063040,   0.226747,  -0.125909,  -0.002306,   0.030601,   0.018530,   0.048114,   0.031152,  -0.014886,  -0.161616,  -0.012648,  -0.078928,   0.077061,  -0.148575,  -0.117602,  -0.024389,
         0.277585,   0.194742,  -0.054415,  -0.200214,   0.038537,   0.028225,   0.068707,   0.037026,   0.052730,  -0.088805,  -0.222825,  -0.105774,  -0.067697,   0.200925,  -0.126375,   0.096630,   0.120438,  -0.151351,  -0.030533,   0.136997,  -0.009933,  -0.334377,   0.311234,  -0.103368,  -0.290904,   0.382416,   0.287526,  -0.083274,   0.011288,   0.022545,   0.084782,   0.184512,
         0.116108,   0.104091,  -0.378349,  -0.117705,  -0.129982,   0.028255,  -0.029116,   0.254313,  -0.234247,  -0.080620,   0.047553,   0.021688,  -0.042463,   0.016162,   0.255293,   0.085733,   0.053483,   0.190647,   0.008896,  -0.043063,   0.287951,   0.057213,  -0.007028,  -0.249107,  -0.189771,   0.102655,  -0.057292,   0.007451,  -0.084728,  -0.098151,   0.123153,   0.357568,
        -0.119500,   0.002058,  -0.106970,  -0.459359,  -0.034183,  -0.150930,  -0.076462,   0.183385,  -0.080062,   0.044993,   0.031603,  -0.025500,   0.182419,   0.147449,   0.140643,   0.187124,  -0.043922,  -0.048792,   0.576830,  -0.054059,  -0.044587,   0.133692,   0.030101,  -0.181973,  -0.244319,  -0.049659,   0.137391,  -0.192610,  -0.123326,   0.013633,  -0.252159,   0.006090,
        -0.091456,  -0.478069,   0.152593,   0.027669,  -0.395405,  -0.238794,  -0.099926,  -0.020739,   0.134997,   0.340495,  -0.307065,   0.226492,   0.071193,   0.103487,   0.120058,   0.286674,   0.044603,   0.023616,   0.365962,   0.632061,   0.409986,  -0.302601,   0.513676,   0.439051,   0.035513,   0.867612,  -0.317106,  -0.558950,   0.353128,  -0.059124,   0.027763,   0.027485,
         0.377679,   0.153183,   0.039482,  -0.050061,  -0.188671,   0.192670,  -0.104594,   0.052196,   0.269196,  -0.060247,   0.012818,   0.142847,   0.090457,   0.006917,   0.314916,   0.180393,  -0.087692,  -0.069551,   0.193292,   0.026763,  -0.361309,  -0.085803,   0.062117,  -0.167580,   0.456056,  -0.008626,  -0.025298,  -0.005941,  -0.079080,  -0.070788,   0.059860,  -0.092377,
        -0.254796,  -0.104374,   0.049256,  -0.031816,  -0.001846,   0.047585,   0.182044,   0.120952,   0.330063,   0.122252,   0.151926,   0.245299,   0.275602,   0.441855,  -0.033795,  -0.364510,  -0.213327,   0.023007,  -0.034446,  -0.019776,   0.113816,  -0.109989,   0.078548,   0.260608,  -0.278475,   0.012998,  -0.204266,  -0.336972,  -0.126367,  -0.240079,   0.074659,   0.184799,
         0.119529,  -0.062964,  -0.130301,   0.177060,   0.043532,   0.026199,   0.125303,  -0.034172,   0.020529,   0.196227,   0.690963,  -0.082522,   0.004518,   0.158393,   0.352700,  -0.263851,  -0.014363,   0.114406,   0.057153,   0.043723,  -0.238923,   0.011000,   0.135845,  -0.018389,  -0.113861,  -0.196731,  -0.111006,   0.056451,  -0.244360,  -0.391255,  -0.121225,  -0.101179,
        -0.188002,   0.094908,   0.011125,   0.084347,   0.065021,  -0.048315,   0.004471,   0.468679,   0.361291,   0.014465,   0.121881,   0.137172,   0.316911,   0.009979,   0.200043,   0.112957,   0.331317,   0.247418,  -0.078665,  -0.304049,   0.218037,  -0.210674,  -0.182341,   0.014665,   0.057785,   0.016380,  -0.103671,  -0.138109,  -0.149511,   0.326462,   1.417680,   0.537558
    };

    //! weights for layer 2
    const double ANN_WEIGHTS_CONTACT_HELIX_HELIX_LAYER_1[ 17] =
    {
        -0.267875,  -0.747771,  -0.894413,   1.028650,   0.707112,   0.693675,   0.478090,  -0.824753,   0.939537,  -0.630738,   0.344558,  -0.740079,   0.659889,   0.557211,   1.061230,  -0.743564,  -0.765808
    };
    //! ANN CONTACT_HELIX_HELIX definition
    double ANN_CONTACT_HELIX_HELIX( const linal::Vector< double> &INP)
    {

      // declare variables
      int nora( 0), norb( 0), wei( 0);
      double *hid;

      // test net size
      BCL_Assert( INP.GetSize() == 543, "wrong input size!");

      // allocate memory
      linal::Vector< double> hidden[ 3];
      hidden[0] = linal::Vector< double>( 543);
      hidden[1] = linal::Vector< double>( 16);
      hidden[2] = linal::Vector< double>( 1);

      // normalize data
      hid = hidden[ 0].Begin();
      for( const double *inp = INP.Begin(); inp != INP.End(); inp++)
        ( *( hid++)) = ( *inp) * ANN_NORMALIZE_CONTACT_HELIX_HELIX_A[ nora++] + ANN_NORMALIZE_CONTACT_HELIX_HELIX_B[ norb++];

      // calculate network
      // calculate layer 1
      wei = 0;
      hid = hidden[1].Begin();
      for( size_t i = 0; i < 16; i++)
      {
        *hid = ANN_WEIGHTS_CONTACT_HELIX_HELIX_LAYER_0[ wei++];
        for( const double *inp = hidden[ 0].Begin(); inp != hidden[ 0].End(); inp++)
          ( *hid) += ANN_WEIGHTS_CONTACT_HELIX_HELIX_LAYER_0[ wei++] * ( *inp);
        *hid = double( 1.0) / ( double( 1.0) + exp( -( ( *hid)))); hid++;
      }

      // calculate layer 2
      wei = 0;
      hid = hidden[2].Begin();
      for( size_t i = 0; i < 1; i++)
      {
        *hid = ANN_WEIGHTS_CONTACT_HELIX_HELIX_LAYER_1[ wei++];
        for( const double *inp = hidden[ 1].Begin(); inp != hidden[ 1].End(); inp++)
          ( *hid) += ANN_WEIGHTS_CONTACT_HELIX_HELIX_LAYER_1[ wei++] * ( *inp);
        *hid = double( 1.0) / ( double( 1.0) + exp( -( ( *hid)))); hid++;
      }

      // denormalize data
      // end
      return hidden[ 2]( 0) * ANN_NORMALIZE_CONTACT_HELIX_HELIX_A[ nora] + ANN_NORMALIZE_CONTACT_HELIX_HELIX_B[ norb];

    }

    //! normalization values A*x+b
     const double ANN_NORMALIZE_CONTACT_HELIX_SHEET_A[ 424] =
    {
       0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.000833333, 0.000833333, 0.000833333, 1.2
    };

    //! normalization values a*x+B
     const double ANN_NORMALIZE_CONTACT_HELIX_SHEET_B[ 424] =
    {
       0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, -0.1
    };

    //! weights for layer 1
    const double ANN_WEIGHTS_CONTACT_HELIX_SHEET_LAYER_0[ 6784] =
    {
         0.105774,   0.075385,   0.077020,  -0.024217,   0.489225,  -0.126658,   0.083324,   0.193281,  -0.317094,   0.816847,   0.441671,   0.107409,  -0.562153,   0.204011,  -0.338469,   0.349911,  -0.320586,  -0.527628,   0.069711,  -0.005750,   0.052190,   0.058665,  -0.154673,  -0.140112,   0.154337,  -0.208273,  -0.052732,  -0.066942,   0.463920,  -0.313384,   0.015087,   0.318622,
         0.292043,   0.040481,   0.259695,   0.176816,   0.061508,   0.055597,   0.048330,  -0.602834,   0.263671,   0.078197,  -0.250392,  -0.216353,  -0.348023,   0.107135,  -0.315775,  -0.474739,   0.401750,   0.204037,  -0.221455,   0.169597,  -0.118857,   0.161480,  -0.160716,   0.032223,   0.115949,  -0.164152,  -0.321349,  -0.174225,   0.425074,   0.008604,  -0.043814,   0.004333,
        -0.093261,   0.021390,   0.135667,  -0.015992,   0.408284,  -0.772105,  -0.196123,  -0.357301,   0.356935,  -0.197092,   0.639488,  -0.643541,  -0.027572,   0.477410,   0.193988,  -0.332055,  -0.285269,   0.224295,   0.607935,  -0.239924,   0.046983,   0.513809,   0.499591,  -0.036471,   0.310165,   0.340415,   0.169783,  -0.163065,   0.248678,   0.146351,   0.284184,   0.419539,
         0.473150,  -0.202983,   0.462663,  -1.624260,  -0.390967,   0.722593,  -0.245905,  -0.124340,  -0.359682,   0.077133,  -0.103395,   0.087462,  -0.197662,  -0.042030,  -0.382632,   0.090607,  -0.546808,   0.129390,   0.204105,   0.243037,  -0.488977,   0.121811,  -0.005514,   0.094813,   0.002596,  -0.098672,   0.157136,  -0.300761,   0.780256,  -0.071081,   0.008534,   0.162590,
         0.336068,  -0.712351,  -0.865031,   0.036230,  -0.573138,  -0.301246,  -1.806790,   0.253232,  -0.305365,  -1.290930,  -0.123481,  -0.454063,  -0.186037,   0.226354,  -0.885507,   0.031423,  -0.174694,  -1.654460,   0.196644,   0.204400,   0.027294,  -0.383076,  -0.391606,  -0.186106,   0.321353,  -0.155739,   0.171688,  -0.298694,   0.169155,   0.057413,   0.273769,  -0.606529,
        -0.868213,   0.474641,  -0.089519,  -0.077064,  -0.419440,   0.126664,   0.199258,  -0.430875,   0.418742,  -0.378751,   0.112104,  -0.076873,   0.095635,  -0.106962,  -0.513625,  -1.012240,   0.112497,  -0.049659,  -0.523274,  -0.100952,   0.220877,   0.201547,   0.211279,  -0.001375,  -0.062703,   0.379687,   0.381095,   0.369836,  -0.106119,  -0.353370,   0.044797,   0.353978,
         0.010933,   0.236838,   0.281550,  -0.459115,   0.375769,   0.281847,   0.164974,  -0.663708,   0.548944,   0.133970,   0.716524,   0.602055,  -0.192905,  -0.271917,   0.148509,   0.382714,   0.139264,   0.333168,   0.557814,  -0.289017,   0.231529,  -0.148644,  -0.162509,  -0.422488,   0.332809,  -0.157820,   0.072987,  -0.945917,   0.338195,   0.353038,  -0.281725,   0.107052,
        -0.534015,  -0.121456,  -0.470291,  -0.393473,  -0.497368,   0.199214,  -0.373002,  -0.094338,  -0.487115,  -0.031423,   0.269697,  -0.746323,   0.337745,   0.646181,  -0.377899,  -0.021583,  -0.304526,   0.122761,   0.038750,  -0.167459,   0.326841,  -0.047556,  -0.090330,   0.228814,  -0.165314,   0.068240,   0.117955,  -0.230519,  -0.463823,   0.136280,  -0.982440,   0.891558,
        -0.276977,  -1.387350,   0.039242,   0.137632,   0.324553,  -0.066546,  -0.265961,  -0.031165,  -0.266137,  -0.179983,  -0.041626,   0.167637,  -0.412900,   0.322432,   0.127493,   0.093886,   0.041314,  -0.020238,   0.338659,  -0.114915,   0.520299,  -0.081493,  -0.316493,   0.524732,  -0.439335,   0.080561,   0.286989,   0.128763,   0.339332,  -0.638414,   0.524211,   0.580734,
         0.239756,   0.131932,  -0.106626,  -0.434520,   0.289918,  -0.318770,  -0.408928,   0.359359,   0.245940,   0.150425,  -0.022579,   0.292135,   0.103146,  -0.124814,   0.078222,   0.291387,   0.094768,  -0.186037,   0.015797,   0.204951,  -0.422685,   0.305383,  -0.147856,   0.073209,   0.076636,   0.309917,   0.693239,   0.048785,  -0.138729,   0.122225,  -0.238618,  -0.219430,
         0.028912,  -0.203055,   0.477101,  -0.636760,  -0.062401,   0.250622,   0.148956,   0.028782,  -0.451394,   0.229759,  -0.354705,   0.065045,   0.191418,   0.089510,   0.184022,   0.144401,   0.011745,   0.425707,   0.216897,  -0.041194,  -0.063522,  -0.201575,   1.084150,   1.517370,   1.873780,   0.129207,   0.595251,   0.185459,   0.147338,   0.713274,  -0.495769,  -0.524842,
         0.933609,   0.018797,   0.169507,   0.193028,   0.419346,   0.551048,   0.329775,   0.593159,  -0.410326,  -0.040387,   0.128503,   0.248112,  -0.057721,   0.252558,  -0.017356,   0.171321,   0.220844,  -0.577184,   0.246643,  -0.114295,   0.001879,  -0.007842,  -0.065047,  -0.162008,  -0.270076,   0.095792,  -0.208156,  -0.487995,  -0.146807,   0.312525,   0.547936,  -0.185512,
         0.004320,   0.875103,  -0.124557,  -0.092623,   0.027812,   0.086921,   0.151403,   0.005989,  -0.272100,  -0.217960,   0.414064,  -0.164063,   0.153511,   0.158091,   0.444513,  -0.855389,   0.020998,  -0.368749,   0.890417,   0.318091,   0.763376,  -0.463187,  -0.026217,   0.101538,  -0.101364,   0.343416,  -0.292080,  -0.044602,   0.393005,  -0.121164,  -0.364122,   0.134642,
         0.191252,   0.222944,   0.349218,   0.407697,   0.098786,   2.837850,  -1.361040,   0.209890,
         0.059661,   0.275002,  -0.093821,  -0.009275,   0.209699,  -0.275447,  -0.003508,   0.215774,  -0.688444,   1.470300,   0.631229,  -0.127055,  -0.696808,  -0.089075,  -0.418259,   0.094684,  -0.736574,  -0.803780,  -0.598808,  -0.475073,   0.309827,   0.225451,  -0.666661,   0.130778,   0.194921,  -0.230910,  -0.147177,  -0.023005,   0.165198,  -0.391719,   0.191793,   0.197284,
         0.127716,   0.009379,   0.211563,   0.043722,   0.051841,   0.131617,  -0.319306,   0.214940,   0.668992,   0.127413,   0.399314,  -0.054675,  -0.508985,   0.141048,   0.007560,  -0.318050,   0.323330,   0.202661,  -0.133404,   0.045804,   0.099086,   0.040368,  -0.372174,  -0.034726,   0.479663,  -0.037988,   0.698575,  -0.269380,   0.139196,   0.299791,   0.172514,   0.099599,
        -0.102903,  -0.120694,   0.062912,   0.326282,  -0.173526,   0.286687,   0.664199,  -0.639210,   0.142693,  -0.022663,   0.468227,  -0.178073,   0.124361,   0.281392,   0.147144,   0.164298,  -0.059851,   0.012410,   0.519844,  -0.225658,   0.222940,  -0.174067,   0.383370,   0.196268,   0.160305,   0.288983,   0.252897,  -0.056664,   0.049630,  -0.026925,   0.035632,   0.168207,
         0.087646,   0.017050,  -0.057368,  -1.364910,   0.850493,   0.105597,  -0.650413,  -0.181992,   0.212625,  -0.182056,   0.053630,   0.003311,   0.171109,  -0.149586,  -0.587463,  -0.032642,  -0.558396,   0.295654,  -0.350915,   0.712825,   0.175238,   0.239390,   0.275703,   0.227128,  -0.192485,   0.107204,   0.472901,  -0.037854,   0.534936,   0.010103,   0.026751,   0.253789,
         0.113966,  -1.388150,  -0.245892,  -0.096208,  -0.281855,  -0.552632,  -0.828771,   0.175901,  -0.498477,  -1.138230,  -0.810156,  -0.270049,  -0.143547,   0.276675,  -0.542200,   0.137668,  -0.236082,  -0.343527,   0.254192,  -0.326238,  -0.094427,  -0.563628,  -0.490823,   0.050118,   0.253126,   0.068035,   0.096477,  -0.140985,   0.148775,   0.169724,   0.476894,  -1.626610,
        -0.470277,   0.169297,   0.217077,   0.434693,   0.458553,   0.005389,   0.236823,  -0.052694,   0.146776,   0.161541,   0.086478,  -0.434802,   0.331304,  -0.159819,  -0.503916,   0.234204,   0.457658,   0.218088,  -0.064668,   0.042009,   0.223857,  -0.276659,   0.334589,   0.012079,  -0.263637,  -0.232796,   0.174548,   0.147783,   0.176560,  -1.444570,  -0.026292,   0.017130,
         0.160201,   0.475837,   0.825393,  -0.268531,   0.412854,   0.619051,   0.026153,   0.167289,   0.160224,   0.099086,   0.440792,   0.075872,   0.248073,  -0.290697,   0.299312,   0.545950,   0.869887,   0.535295,   0.470943,  -0.017253,   0.363994,  -0.119991,   0.125973,  -0.160952,   0.110071,  -0.042212,   0.215058,  -1.780790,   0.088348,   0.094704,  -0.806773,  -0.387082,
        -0.511874,   0.228341,  -0.700335,  -0.564748,  -0.640880,  -0.048320,  -0.093959,   0.360449,  -0.936694,  -0.073979,   0.110339,  -0.834982,   0.166611,   0.403120,  -0.118528,   0.392773,   0.071114,   0.013831,  -0.092544,  -0.084002,   0.046568,  -0.149509,  -0.340691,  -0.028083,   0.148421,  -0.884424,  -0.115626,  -0.527439,  -0.603599,  -0.141756,  -0.780834,  -0.263155,
        -0.103762,  -0.760367,   0.198993,  -0.361309,   0.112016,  -0.183856,  -0.346789,  -0.070187,  -0.256163,   0.489984,  -0.214927,  -0.104250,  -0.519698,  -0.156329,  -0.138343,   0.226683,   0.420045,   0.253991,  -0.026710,   0.007990,   0.265837,  -0.103048,  -0.299388,   0.309825,  -0.111067,   0.233478,   0.632414,   0.177047,   0.277633,   0.155266,   0.480816,   0.421570,
         0.279124,  -0.028562,   0.032945,  -0.400245,   0.785041,  -0.177403,  -0.214025,   0.394334,   0.372455,  -0.009064,  -0.325412,   0.226556,  -0.070067,   0.071266,   0.001093,   0.176383,   0.155819,   0.078636,   0.030574,   0.090504,  -0.294345,   0.127876,  -0.295010,  -0.011027,   0.086484,   0.281640,   0.175083,  -0.535270,   0.026461,  -0.202069,  -0.321883,  -0.006876,
         0.382941,   0.088891,   0.269794,  -0.346871,  -0.139488,   0.796095,  -0.028804,   0.120456,  -0.228211,  -0.130830,  -0.315338,   0.151281,  -0.055518,  -0.220435,   0.189866,   0.182407,   0.282073,   0.215571,   0.660765,  -0.226847,  -0.159432,   0.023666,   0.256911,   0.066178,  -0.003883,   0.062784,  -0.239829,   0.026939,  -0.506195,  -0.433640,  -0.082635,  -0.099057,
         0.228840,  -0.095093,   0.173947,   0.636317,   0.400308,   0.584865,  -0.143819,   0.171085,   0.042192,   0.198155,   0.148066,   0.234196,   0.215141,  -0.015955,  -0.038817,   0.225060,  -0.298595,  -0.359506,   0.214507,   0.009587,  -0.158265,   0.470326,   0.553674,  -0.282664,  -0.026208,   0.195680,  -0.297851,  -0.367813,  -0.004961,  -0.082388,   0.269367,  -0.008424,
        -0.145621,   0.664871,   0.125130,  -0.011374,  -0.024646,   0.272800,   0.171271,   0.099636,  -0.137231,  -0.072809,   0.473915,   0.055445,  -0.173451,   0.237699,  -0.066412,  -0.439487,  -0.186818,  -0.075842,   0.658756,   0.142411,   0.044014,  -0.235866,  -0.130468,   0.087921,   0.311696,   0.041154,  -0.157571,  -0.026174,   0.096946,  -0.179418,  -0.243158,  -0.193472,
         0.171159,  -0.188116,  -0.099588,   0.169558,   0.161248,   1.085740,  25.497000,  -0.019668,
         0.026470,  -0.004839,   0.154295,   0.273320,  -0.123372,  -0.016559,   0.075240,   0.060770,  -0.113277,   0.107462,   0.194480,   0.157977,  -0.110732,  -0.345547,   0.212365,  -0.340091,   0.027418,   0.438613,   0.342956,  -0.267198,   0.068980,   0.047181,   0.160611,   0.159069,   0.433022,  -0.021436,  -0.081800,  -0.139312,   1.040310,   0.115911,   0.016756,   0.139605,
        -0.494456,  -0.148273,  -0.195652,  -0.049671,  -0.092259,   0.103876,  -0.354870,   0.775876,   0.390118,   0.171263,   0.022411,   0.380765,   0.513657,  -0.348122,  -0.099665,   0.376611,   0.268858,   0.030128,   0.271206,   0.011486,   0.236367,  -0.053337,  -0.044692,   0.251771,   0.090463,   0.211341,  -0.025768,   0.193095,   0.163023,  -0.399313,  -0.080065,  -0.036479,
        -0.153874,  -0.023945,  -0.095731,  -0.216770,  -0.181693,   0.159455,   0.853715,   0.246292,  -0.221387,   0.178694,   0.214270,   0.262364,  -0.257020,   0.041611,   0.214407,   0.265582,  -0.370135,  -0.229757,  -0.057963,  -0.322098,  -0.096509,  -0.028319,   0.153086,   0.180830,   0.092806,  -0.299254,  -0.198957,   0.128597,   0.078065,   0.082371,  -0.066076,   0.131281,
         0.266460,   0.170990,  -0.071003,   0.838602,   0.070670,   0.681427,   0.231602,   0.119180,  -0.324409,   0.586279,  -0.097730,   0.042299,  -0.165862,  -0.256416,  -0.005347,  -0.342779,   0.119670,  -0.090793,  -0.196675,  -0.468127,   0.100924,   0.276573,   0.050100,  -0.373804,   0.058025,   0.228454,   0.007156,   0.143826,   0.199907,   0.155246,   0.134988,   0.100185,
         0.153085,   0.318899,  -0.454774,   0.389437,  -0.108450,  -0.052767,  -0.044890,   0.093352,  -0.431503,  -0.118828,   0.748681,  -0.298826,  -0.051061,  -0.040454,   0.009374,  -0.084749,   0.353272,   0.071574,  -0.127165,   0.063939,  -0.106633,   0.198206,   0.069794,   0.141416,  -0.203880,   0.042507,   0.275709,  -0.023334,  -0.056942,  -0.060721,   0.019898,  -0.281708,
        -0.066968,  -0.202080,   0.422020,   0.224796,   0.327142,  -0.015952,   0.411564,   0.472575,   0.174638,   0.233400,   0.184303,   0.094087,   0.745663,  -0.135654,   0.124720,   0.003867,   0.134526,   0.243476,  -0.435888,  -0.146896,   0.051464,   0.097864,  -0.271541,   0.084638,   0.046684,  -0.045848,  -0.217860,   0.008789,   0.050711,  -0.913149,   0.144628,   0.198439,
        -0.386497,  -0.111103,   0.082569,  -0.392472,  -0.234607,  -0.517162,   0.518811,   0.334150,  -0.242190,  -0.499403,  -0.282443,  -0.461406,  -0.163261,   0.300786,   0.214016,   0.063250,  -0.631381,   0.321024,  -0.119476,   0.166393,  -0.382818,  -0.105693,   0.048900,   0.038122,   0.279298,   0.078618,   0.239382,   2.916330,  -0.283901,   0.152021,   0.055787,  -0.498644,
         0.190311,   0.223045,  -0.246449,  -0.100299,   0.092173,  -0.130983,   0.102260,   0.293848,   0.233477,  -0.310315,  -0.205186,  -0.134273,  -0.476545,   0.074071,   0.044491,  -0.181810,   0.244335,   0.151643,  -0.142633,   0.037128,   0.152023,   0.088527,   0.376318,   0.066833,   0.428979,  -0.658012,  -0.389421,   0.288221,  -0.011318,  -0.313885,  -0.183003,   0.440833,
        -0.250218,   0.466759,  -0.512193,  -0.041660,   0.273330,   0.397929,   0.288846,   0.038163,  -0.251306,  -0.132060,  -0.271599,  -0.133111,  -0.212610,  -0.068783,   0.174386,   0.332033,  -0.272463,   0.003433,   0.477726,   0.168288,  -0.011943,   0.359400,  -1.887210,   0.601717,  -4.450630,   0.206672,  -0.361701,  -0.295456,  -0.166116,  -0.117654,  -0.399731,   0.311585,
        -0.469380,  -0.388942,   0.170070,  -0.152498,  -0.078982,   0.464674,   0.320790,   0.343150,   0.444427,   0.620662,  -0.439813,  -0.005079,   0.126022,  -0.265285,   0.090487,  -0.133525,  -0.200083,  -0.017860,   0.082140,  -0.165798,  -0.277138,  -0.011531,   0.208717,  -0.070585,   0.517057,   0.635950,   1.227760,  -0.189910,   0.030296,   0.732231,   0.105090,   0.496779,
        -0.586091,  -0.140153,   0.189852,  -0.024356,  -0.472181,  -0.152733,   0.613052,   0.337570,   0.085961,  -0.062353,   0.093035,   0.485001,   0.343947,   0.365229,   0.221521,  -0.375334,  -0.207801,   0.211670,  -0.539068,   0.334225,  -0.388322,  -0.248059,  -0.824399,   0.323576,  -0.078357,  -0.180604,  -0.256994,  -0.400719,  -0.622179,  -0.695035,   0.329872,   0.243423,
        -0.433791,   0.364311,   0.024289,  -0.056891,  -0.013864,   0.520560,   0.714555,   0.111236,  -0.008223,  -0.131773,   0.420362,   0.244587,  -0.166271,   0.277255,  -0.171041,  -0.134136,  -0.042884,  -0.095767,   0.006209,  -0.330259,   0.858589,   0.120294,   0.594143,  -0.094134,   0.293865,   0.614742,  -0.020219,   0.545975,   0.178186,  -0.058920,   0.410492,   0.275467,
         0.477552,  -0.996493,   0.106895,  -0.092359,   0.100600,   0.362439,  -0.068972,   0.145787,  -0.225883,   0.049571,  -0.015947,  -0.181290,   0.179646,  -0.054897,  -0.900602,   0.154591,  -0.023731,   0.008334,  -0.457874,   0.043043,  -0.266922,   0.038006,   0.002983,   0.459806,  -0.206105,  -0.129754,   0.250991,   0.021771,   0.303537,   0.383335,   0.196241,   0.908959,
         0.033536,   0.071304,   0.470551,   0.297092,  -0.118236,  -0.791766,  -0.194075,  -0.313479,
         0.207900,   0.025226,   0.129581,   0.102419,  -0.052944,   0.094750,   0.026500,   0.007845,   0.372917,   0.217714,  -0.011324,  -0.056849,  -0.455293,   0.214742,   0.039784,  -0.171742,  -0.160305,   0.067872,   0.434808,   0.097919,  -0.003615,  -0.099491,  -0.307423,  -0.213974,  -0.009017,   0.045978,   0.265716,   0.157244,   0.463464,   0.161687,  -0.023815,   0.165661,
        -0.142273,   0.117728,   0.275608,   0.206668,   0.073743,   0.005265,   0.319689,   0.144641,   0.025543,   0.139764,  -0.035897,  -0.100682,   0.064707,  -0.167539,  -0.371273,  -0.315666,   0.281372,  -0.200853,   0.142939,  -0.111205,  -0.595896,  -0.218432,  -0.087901,   0.181530,   0.172889,   0.149469,  -0.040450,  -0.260329,  -0.310797,   0.063633,   0.013343,   0.084608,
         0.209151,   0.322528,   0.123101,   0.274240,   0.380309,  -1.027210,   0.168161,   0.248022,  -0.055533,  -0.010814,   0.021681,  -0.233437,  -0.000228,  -0.114572,  -0.248260,  -0.102217,   0.179925,  -0.177655,   0.098758,  -0.012738,  -0.189238,  -0.251013,  -0.259903,  -0.353308,  -0.368318,  -0.467657,  -0.035619,  -0.022628,   0.066215,   0.237719,   0.027572,   0.121282,
         0.047246,   0.073813,   0.185689,  -0.357976,   0.453588,  -0.141525,   0.197416,   0.230372,   0.494214,  -0.694982,   0.178379,   0.272638,   0.375259,   0.418474,   0.034232,   0.055447,   0.059313,  -0.016331,  -0.067756,  -0.118347,   0.531497,   0.203101,   0.684721,   0.173791,  -0.016718,   0.155520,  -0.030863,   0.269833,  -0.092525,   0.143520,   0.159813,   0.211194,
        -0.025906,   0.011888,   0.182809,   0.229193,  -0.154629,  -0.191306,  -0.057846,  -0.066838,  -0.039700,   0.016053,  -0.271120,  -0.115250,   0.145118,   0.145968,  -0.346322,   0.134904,   0.277342,  -0.110521,   0.145492,   0.148293,  -0.132274,  -0.111817,  -0.013101,   0.105087,   0.067516,   0.263567,   0.307940,   0.282083,   0.224673,   0.224093,   0.011155,  -0.363761,
         0.558323,   0.196349,  -0.270046,  -0.209559,  -0.222824,  -0.199894,   0.008796,  -0.179589,  -0.075159,   0.034243,   0.171130,   0.568980,  -0.362966,   0.125227,   0.278724,  -0.115331,  -0.167632,   0.038974,  -0.044600,   0.453955,   0.269290,  -0.121342,  -0.150003,   0.045012,  -0.050415,   0.222506,   0.062994,  -0.057617,   0.127191,  -0.203909,   0.274848,   0.226698,
         0.247435,  -0.004663,   0.216617,  -0.082723,   0.058544,   0.245986,   0.053128,  -0.029514,   0.068617,   0.330195,   0.416320,  -0.067637,  -0.498470,   0.130421,   0.338681,   0.074490,  -0.463484,  -0.284827,  -0.066276,   0.014557,   0.194612,   0.290277,   0.065890,   0.331926,   0.133110,  -0.020321,   0.119517,  -0.230195,   0.228138,  -0.089869,   0.170564,   0.559147,
         0.291717,  -0.325762,   0.286077,   0.054865,   0.364485,   0.019052,   0.103273,   0.095445,   0.160666,  -0.011999,  -0.334812,  -0.056975,   0.119002,  -0.058735,   0.079157,  -0.214741,  -0.483228,  -0.092443,  -0.178857,  -0.015614,   0.101616,   0.221219,   0.259068,   0.157865,  -0.301617,  -0.033133,   0.785520,   0.112000,   0.164571,   0.139434,  -0.021691,  -0.626702,
         0.073332,  -0.220283,   0.275321,  -0.035370,  -0.042553,   0.014486,   0.102945,   0.057389,   0.225839,   0.069813,   0.217603,  -0.048838,   0.067182,   0.187782,  -0.000892,  -0.076341,   0.232688,   0.219011,  -0.170965,   0.284338,   0.244191,  -0.168467,   0.623843,   0.273425,  -0.033354,  -0.321134,   0.689478,   0.603482,  -0.032238,  -0.409826,   0.042556,   0.330375,
         0.260586,   0.439579,   0.108997,  -0.182281,   0.420471,   0.058445,   0.151890,  -1.068990,  -0.235370,  -0.249167,  -0.223046,  -0.090958,  -0.352541,   0.021365,  -0.115507,  -0.068078,   0.183629,   0.262394,   0.228887,  -0.062762,   0.407208,   0.180341,   0.299419,  -0.213957,  -0.810570,  -0.014248,  -0.395737,   0.193058,  -0.450526,  -1.134170,   0.200522,  -0.423142,
         0.018540,  -0.021628,  -0.967606,   0.199239,  -0.034415,  -0.433124,  -0.159821,   0.045250,  -0.096829,  -0.370765,  -0.194722,  -0.068788,   0.035598,   0.342145,  -0.239339,   0.351541,   0.019382,  -0.010039,   0.231265,   0.339043,   0.130281,  -0.234254,   0.099356,   0.260994,   0.307163,  -0.276142,   0.042687,  -0.109530,   0.037840,  -0.138418,  -0.336762,  -0.232474,
         0.241426,  -0.076746,   0.006890,  -0.441820,  -0.114578,  -0.021001,  -0.038093,  -0.098529,  -0.695339,   0.006896,  -0.191903,   0.073648,   0.318664,   0.097122,   0.231154,  -0.002144,   0.403392,   0.117913,   1.078250,   0.385948,  -0.093760,   0.157058,  -0.508456,  -0.111271,  -0.556358,  -0.787274,   0.078286,  -0.179034,  -0.251324,  -0.189177,  -0.435835,  -0.256547,
        -0.155247,   1.290210,  -0.408605,  -0.405602,  -0.032853,  -0.086489,   0.006618,   0.316925,  -0.021525,   0.189234,   0.048248,   0.168054,   0.059230,   0.310735,   0.608679,  -0.253017,   0.788420,  -0.142549,   0.116330,   0.145838,  -0.348428,  -0.104076,  -0.210174,  -0.334751,   0.033624,  -0.210198,   0.294827,   0.147183,  -0.326909,  -0.072910,   0.156703,  -0.655281,
        -0.280164,   0.248559,  -0.446664,  -0.188677,   0.254104,   0.564591,   6.413870,   0.974492,
        -0.048195,  -0.059712,   0.019361,  -0.220543,   0.234169,  -0.506036,  -0.076346,  -0.179696,  -0.703501,   0.596753,   0.706517,  -0.009681,  -0.710933,   0.157001,   0.065914,  -0.117011,  -0.805352,  -0.008073,   0.193634,   0.206526,   0.214194,  -0.362391,  -0.142195,  -0.352116,  -0.050210,   0.071675,  -0.024678,   0.203116,   0.561579,  -0.539266,   0.300648,   0.161600,
         0.011598,  -0.010105,  -0.137652,  -0.087627,  -0.019854,   0.133600,  -0.501307,  -0.538665,   0.828511,  -0.039399,   0.170107,  -0.225369,  -0.256834,   0.523421,  -0.052418,  -0.133508,   0.512939,  -0.056874,  -0.413889,  -0.445461,   0.437316,   0.092877,   0.046700,   0.316704,   0.453201,   0.385250,   0.284561,  -0.158620,   0.425943,   0.079686,   0.032523,  -0.114064,
        -0.145914,  -0.071919,  -0.012033,   0.132144,  -0.055794,  -0.534451,   0.108823,  -0.060310,  -0.023395,   0.069668,   0.911746,  -0.528028,   0.056020,   0.168737,   0.787157,   0.118363,  -0.013016,  -0.148875,  -0.084938,  -0.022426,  -0.027324,   0.965075,   0.527464,   0.572514,   0.285388,   0.618793,   0.387780,   0.222050,   0.307912,   0.125749,   0.104918,   0.032901,
         0.235884,   0.038659,   0.765968,  -3.051530,  -0.503448,   0.086969,  -0.381095,  -0.136239,   0.839522,   0.296503,  -0.064440,  -0.452872,  -0.206010,  -0.872169,   0.029251,  -0.011465,  -0.797792,  -0.158375,   0.319783,   0.207570,  -0.274487,  -0.072025,   0.083981,  -0.028382,  -0.320911,   0.180938,   0.513660,  -0.015105,   0.204907,   0.421360,   0.477750,   0.143924,
         1.001920,  -1.883140,  -1.943950,   0.555881,   0.451650,   0.264137,   0.205662,   0.366123,   0.060601,  -0.478807,  -0.584356,  -0.235677,   0.070536,   0.406743,  -0.075087,   0.332797,   0.082146,   0.089162,   0.425354,   0.350252,  -0.010274,  -0.032060,  -0.051700,  -0.322235,   0.611840,  -0.055785,  -0.337508,  -0.536892,   1.108970,  -0.179095,   0.915712,  -2.009970,
        -2.027550,   0.699182,   0.172200,   0.046583,  -0.588994,   0.140461,   0.153090,   0.221282,  -0.318318,   0.156303,   0.218663,   0.038690,  -0.260815,   0.492041,   0.498316,  -0.345405,   0.633389,   0.242489,  -0.056016,   0.437754,  -0.004405,  -0.202976,   0.407841,  -0.032324,  -0.458028,   0.044845,   0.506662,   0.040196,   1.218750,  -1.605920,  -1.561500,   0.166390,
         0.344360,   0.390499,   0.114694,   0.373506,   0.307469,   0.061145,  -0.032239,   0.103054,  -0.464145,  -0.060024,   0.245619,   0.337971,  -0.491730,  -0.474997,   0.271557,  -0.118713,  -0.066565,   0.037991,   0.198997,  -0.141507,   0.596156,  -0.056267,  -0.090964,  -0.059670,  -0.285375,  -0.084949,   0.894874,  -2.276030,  -0.486360,  -0.464996,   0.248388,   0.564198,
        -0.187342,   0.429108,   0.847348,  -0.079898,   0.364662,   0.139649,  -0.550792,   0.000599,  -0.341263,   0.174288,   0.075928,  -0.557177,   0.051488,   0.480442,  -0.454445,   0.731250,  -0.473451,  -0.182664,   0.038386,  -0.126667,  -0.100129,  -0.061005,  -0.019079,   0.070495,   0.497450,  -1.710900,  -0.153026,  -0.369115,  -0.399135,   0.366049,   0.208323,  -0.593800,
         0.170558,  -0.095921,   0.333274,   0.535308,  -0.138310,  -0.678689,  -0.307195,   0.185798,  -0.147978,  -0.854580,   0.265834,   0.385626,  -0.132890,   0.340855,   0.050415,   0.047357,   0.389509,   0.085609,   0.060300,  -0.158460,   0.253857,  -0.048348,  -3.042390,   1.143950,  -0.318925,  -0.242996,   0.030344,   0.071879,  -0.024807,   0.179252,   0.524082,  -0.205354,
        -0.450068,  -0.162161,   0.214557,  -0.083702,   0.140628,   0.242847,   0.426780,   0.326326,   0.333112,   0.351136,   0.182663,  -0.054340,  -0.172925,  -0.387771,  -0.482525,  -0.167049,   0.236096,  -0.161815,   0.026022,  -0.161769,  -3.196280,   0.510700,   0.220656,   0.076607,  -0.273288,   0.671132,  -0.228958,   0.262830,   0.261182,  -0.329580,   0.116345,  -0.124028,
         0.268978,   0.003804,  -0.188514,   0.340141,  -0.018869,   0.295259,  -0.041389,   0.208578,   0.232054,   0.134153,  -0.166921,   0.099210,   0.027406,  -0.116721,  -0.015516,   0.062973,   0.140720,   0.128586,  -1.877400,   0.000025,   0.590009,  -0.299077,   0.805608,   0.707499,   0.802057,   0.178611,   0.792647,   0.507833,  -0.271333,  -0.049941,  -0.227776,  -0.400155,
         0.557957,   0.083962,  -0.562986,  -0.090226,   0.448489,  -0.194959,  -0.251450,  -0.231051,  -0.555758,   0.017754,   0.439401,   0.273136,  -0.285449,  -0.066702,  -0.104406,   0.073850,  -2.667780,  -0.123785,   1.140080,  -0.348812,  -0.299958,   0.270224,   0.373271,  -0.304509,  -0.173610,   0.404720,  -0.230592,   0.008800,  -0.220994,  -0.693739,   0.217413,   0.020788,
         0.069050,   0.865370,   0.090628,  -0.159998,  -0.038709,   0.169038,   0.129400,   0.218676,   0.017312,  -0.179498,   0.467720,   0.174401,   0.138366,   0.037830,  -1.631560,  -0.508529,   0.539743,  -0.142884,   0.224376,   0.388804,   0.156472,   0.038326,   0.160013,   0.040607,  -0.011324,  -0.630925,  -0.301661,   0.253969,   0.187224,   0.214932,  -0.022037,  -0.011213,
         0.192186,   0.076709,   0.302770,   0.039999,  -0.053609,   0.674826,   0.152645,   1.595380,
         0.115658,   0.287969,   0.199822,   0.310379,  -0.025473,   0.084066,   0.119631,   0.201900,   0.334867,  -1.183970,  -0.048821,   0.150475,   0.049711,  -0.018198,  -0.319360,  -0.072157,  -0.355294,   0.249692,  -0.141280,  -0.126239,   0.310392,   0.353783,   0.068605,   0.236483,  -0.024516,   0.261722,   0.133092,   0.159673,  -0.353558,  -0.095612,   0.236032,   0.182706,
        -0.000941,   0.176869,   0.197956,  -0.240829,   0.260317,   0.029250,   0.146590,  -0.719901,  -0.077715,   0.056073,   0.167372,  -0.583213,  -0.243943,  -0.297854,  -0.345107,  -0.085258,   0.142376,  -0.069112,   0.196721,   0.179197,  -0.054055,   0.248319,   0.072878,   0.094756,  -0.114570,   0.008751,   0.354753,  -0.300661,  -0.057196,   0.232451,   0.188889,   0.097578,
         0.230611,  -0.009017,   0.179551,   0.324585,   0.282631,  -0.362594,  -0.358914,   0.129675,   0.094714,  -0.285706,  -0.203597,  -0.334566,  -0.345511,  -0.377762,  -0.521043,  -0.137313,   0.463979,  -0.136079,   0.205435,   0.060713,   0.301636,  -0.401902,  -0.351556,  -0.423857,  -0.117676,  -0.036474,   0.038077,   0.001871,   0.070312,   0.200427,   0.169404,  -0.099668,
        -0.144000,   0.125880,   0.017201,   0.433077,   0.240802,  -0.223658,  -0.134011,  -0.340812,   0.381052,  -0.355265,   0.025991,   0.061814,   0.190737,  -0.204474,   0.078863,   0.211551,  -0.431732,  -0.002019,   0.158524,   0.040382,   0.262752,  -0.042966,  -0.198053,   0.140647,   0.223341,   0.236270,   0.326676,   0.475758,  -0.275422,   0.288766,   0.043462,   0.009136,
         0.128119,  -0.727492,   0.081717,  -0.466076,   0.415886,   0.195295,  -0.185382,  -0.269429,   0.124871,   0.088509,  -0.664464,   0.272808,   0.271983,   0.237466,   0.529422,   0.093827,   0.215501,  -0.314959,  -0.249482,  -0.474634,  -0.193057,   0.266470,  -0.148574,   0.099330,   0.054835,   0.199824,  -0.000518,   0.354936,   0.139504,  -0.140482,  -0.011368,   0.037861,
         0.100049,   0.217446,   0.076179,  -0.264998,   0.020053,  -0.522935,  -0.064351,  -0.132807,   0.086255,   0.034495,   0.086000,   0.091051,  -0.108017,   0.162174,  -0.087613,  -0.455638,  -0.015473,  -0.004512,  -0.310236,  -0.426688,   0.070206,   0.072748,  -0.012851,  -0.012097,   0.168233,  -0.026728,   0.109894,   0.157951,  -0.030045,   0.079611,   0.258781,   0.133559,
        -0.232880,  -0.530388,  -0.079310,   0.061858,  -0.214704,  -0.202362,  -0.346906,  -0.142603,   0.188545,  -0.076936,  -0.217403,  -0.241522,  -0.014805,  -0.573480,   0.127953,  -0.115074,   0.538159,  -0.068372,  -0.013442,  -0.145012,  -0.391635,  -0.085614,  -0.025465,  -0.000635,  -0.146942,   0.038073,  -0.012734,   0.022894,   0.380403,   0.180369,   0.233626,   0.329908,
         0.198914,   0.126169,   0.224275,  -0.015914,   0.238928,   0.223581,   0.215834,   0.384635,  -0.081221,   0.316221,  -0.077692,  -0.137750,   0.159669,   0.022523,  -0.019791,  -0.185801,  -0.021162,  -0.043387,   0.017797,   0.093724,  -0.034965,  -0.043638,  -0.230485,   0.047616,  -0.393375,   0.935425,   0.619895,  -0.611210,   0.356524,   0.230701,   0.705092,  -0.682297,
         0.168363,   0.197321,  -0.153834,   0.589469,   0.053609,  -0.288649,   0.427568,  -0.302964,   0.043490,   0.025870,  -0.122597,  -0.190786,   0.157641,   0.061183,  -0.125051,  -0.095072,  -0.121241,   0.137295,  -0.294181,   0.263541,   0.096176,  -0.132503,   0.367372,   0.039858,   0.128549,  -0.600142,   0.560691,   0.538459,  -0.821441,  -0.469597,   0.027035,  -0.144408,
         0.427702,   0.401042,  -0.182299,  -0.483635,   0.258841,   0.035389,  -0.035720,  -0.060421,   0.143158,   0.096807,  -0.083044,  -0.043825,  -0.721223,  -0.169942,  -0.345733,  -0.093077,  -0.267674,   0.251215,   0.167855,  -0.211585,  -0.169769,   0.171226,   0.673545,  -0.140588,   0.162007,  -0.067265,   0.520382,  -0.303418,   0.515619,   0.143307,  -0.341731,   0.215548,
        -0.666110,  -0.360215,   0.108545,  -0.035379,  -0.295752,   0.361691,   0.011063,   0.302307,  -0.061398,  -0.224136,  -0.160399,   0.334680,   0.577715,   0.428530,  -0.017423,   0.278307,   0.203321,   0.218810,  -0.144933,   0.315658,  -0.263757,  -0.698845,   1.067060,   0.215279,  -0.206807,  -0.116450,   0.458485,   0.355023,  -1.030810,   0.200711,   0.146335,  -0.141363,
         0.834375,   0.398750,  -0.106631,  -0.639633,  -0.303961,   0.315744,   0.565562,   0.372649,  -0.194678,   0.025088,   0.209078,   0.240930,  -0.057345,   0.225313,  -0.156220,   0.080225,   0.671417,  -0.008697,   0.423078,   0.224676,   0.627730,   0.472742,   0.473437,  -0.439929,  -0.010160,   0.021524,  -0.169683,   0.035368,  -0.101159,   0.247405,   0.133925,  -0.280530,
         0.064118,   0.602051,  -0.003625,   0.164682,  -0.582999,   0.132853,   0.084123,   0.217672,   0.333600,   0.222079,   0.207621,   0.133614,   0.112370,   0.203374,  -0.099847,   0.001860,   0.555103,  -0.223765,   0.662927,  -0.171239,  -0.303856,  -0.397019,  -0.171664,   0.169049,   0.362948,   0.091474,   0.178228,  -0.026737,   0.297468,   0.204969,   0.118715,  -0.479869,
        -0.305904,   0.124747,   0.001172,   0.074192,   0.034007,   0.265217,   3.163840,   0.366257,
         0.028035,  -0.049547,   0.432049,   0.103468,   0.062125,   0.231623,   0.040399,  -0.138936,   0.149661,  -0.914472,  -0.196664,  -0.025695,   0.471337,   0.190706,   0.252027,  -0.200782,   0.237638,   0.344462,  -0.197928,   0.383327,  -0.159045,  -0.067096,   0.475984,  -0.142912,  -0.022858,   0.109362,   0.078793,   0.284838,   0.073557,   0.078547,  -0.128331,   0.087047,
         0.166239,   0.129891,  -0.041577,  -0.065772,   0.039217,   0.123761,   0.195449,  -1.051950,  -0.441980,  -0.195719,  -0.039215,   0.068021,   0.632198,   0.073681,   0.041455,   0.281922,  -0.177252,   0.198431,  -0.370005,  -0.384613,   0.013425,  -0.235707,   0.077817,   0.174342,   0.290987,   0.464463,   0.667759,  -0.164240,  -0.304373,   0.176567,  -0.109467,  -0.168205,
         0.051916,   0.074831,   0.101554,   0.069581,   0.085748,  -0.661222,   0.031993,   0.161743,  -0.382138,   0.226172,  -0.290284,  -0.114993,   0.199275,  -0.316856,   0.048044,  -0.222429,  -0.342426,  -0.058796,  -0.098462,   0.129612,  -0.353546,   0.103747,   0.207806,   0.422783,   0.355455,   0.169821,  -0.089728,  -0.059942,   0.258091,   0.021707,   0.042517,   0.424375,
         0.180478,  -0.009782,   0.103421,  -1.581640,  -0.079122,  -0.370237,  -0.026757,  -0.058047,  -0.201396,  -0.136355,   0.147754,   0.200681,   0.171945,  -0.048034,  -0.177205,   0.201439,   0.463906,   0.152595,  -0.209026,  -0.465587,   0.094895,  -0.198369,   0.090766,   0.094753,   0.080354,  -0.169055,  -0.017508,  -0.069925,   0.037565,  -0.007352,   0.176008,  -0.096274,
        -0.032437,  -0.537192,   0.599499,   0.014002,   0.109878,   0.863031,   2.366100,  -0.250221,   0.943209,   0.801221,  -0.272522,   0.937698,  -0.862272,  -0.285371,   0.534714,   0.423537,  -0.160734,   0.233799,   0.271704,   0.370529,  -0.343923,   0.315323,  -0.299564,   0.079114,  -0.176895,  -0.209856,   0.151119,   0.055499,   0.133181,   0.144386,   0.025640,  -0.925233,
         1.110140,   0.110747,  -0.027771,   0.392468,  -0.083026,  -0.229785,   0.026984,   0.095588,  -0.111417,   0.024633,  -0.443086,   0.104653,  -0.111768,   0.159790,   0.194541,   0.103133,   0.227876,   0.392372,   0.350006,   0.139819,  -0.169552,  -0.331538,   0.142183,  -0.101153,   0.208596,   0.087514,   0.219776,  -0.142612,  -0.040387,  -0.206827,   0.648058,  -0.270518,
        -0.138582,  -0.143686,  -0.115476,  -0.442046,  -0.021119,  -0.257066,   0.408853,  -0.089775,  -0.192262,   0.008256,  -0.234770,  -0.232632,  -0.173426,  -0.555636,  -0.065667,  -0.244061,  -0.256672,   0.061123,   0.143525,   0.014681,   0.383298,   0.164543,   0.203850,   0.010526,   0.022468,   0.077133,  -0.147002,   0.066686,  -0.156384,  -0.400375,  -0.039881,   0.061627,
         0.108053,  -0.033477,   0.283689,   0.109036,   0.372264,  -0.118920,  -0.202512,  -0.188068,   0.209334,  -0.268978,  -0.256192,   0.056019,   0.161847,  -0.012259,  -0.332454,  -0.401877,  -0.012185,   0.068696,   0.029613,   0.052038,  -0.197410,   0.070809,   0.030151,  -0.006698,  -0.414822,   0.948395,  -0.898571,   0.102729,   0.114417,   0.152812,   0.100034,  -0.048976,
         0.210816,   0.497382,  -0.087153,   0.017471,  -0.023550,  -0.173724,  -0.370007,   0.146873,   0.026607,   0.024786,   0.434633,   0.311396,  -0.237344,   0.001645,   0.073525,   0.192053,  -0.283846,  -0.238066,   0.500232,   0.200277,  -0.308226,   0.384070,   0.283592,  -1.604040,   0.977954,   0.340763,  -0.902826,  -0.336823,  -0.286545,  -0.136275,  -0.326587,  -0.794914,
         0.329607,  -0.186087,   0.233238,   0.396135,  -0.307622,  -0.100833,   0.151659,  -0.238116,  -0.257843,  -0.061332,  -0.002592,  -0.356251,   0.375721,   0.268128,   0.037626,   0.117236,  -0.023507,   0.131102,   0.384469,   0.219431,   0.114165,   0.056970,  -0.246332,   0.587472,   1.039790,   0.083811,  -0.629779,   0.413328,   0.299944,   0.730644,   0.198692,   0.165539,
         0.721749,   0.195684,   0.524449,   0.184543,  -0.218278,  -0.529616,   0.676972,   0.248503,   0.114559,   0.303431,   0.106515,  -0.067535,  -0.291009,  -0.513325,   0.516024,  -0.041409,   0.423974,   0.110359,  -0.504258,  -0.168071,  -0.327040,   0.693628,  -1.500030,  -0.647108,  -1.344870,   0.375084,  -0.436312,  -0.728525,  -0.087095,  -0.600687,   0.379031,   0.223088,
        -1.289030,  -0.050192,   0.334355,  -0.414361,  -0.554678,  -1.006730,  -0.790950,  -0.022464,   0.505343,   0.115094,   0.098286,  -0.017028,   0.011013,   0.039823,   0.042040,   0.300792,  -1.182750,  -0.093824,  -0.480185,  -0.037198,  -0.047199,  -0.030354,  -0.427658,   0.316487,   0.659720,  -0.235084,  -0.305047,   0.636651,   0.029513,  -0.114072,  -0.209658,   0.106504,
         0.025087,  -0.807472,   0.548884,   0.396560,   0.179565,   0.045350,  -0.076940,  -0.045038,   0.289905,   0.298082,   0.054278,   0.083987,   0.261212,  -0.084120,  -0.589574,  -0.089445,  -0.480780,   0.367704,  -0.214883,   0.010133,  -0.433689,   0.253317,  -0.098012,  -0.492000,   0.055945,  -0.768449,  -0.180757,   0.165694,  -0.545820,   0.259543,   0.291109,   0.178099,
         0.365567,  -0.329376,  -0.334896,  -0.186661,  -0.133363,   0.460910,   1.721700,  -0.340697,
         0.191162,   0.086394,   0.062904,   0.163729,   0.083245,   0.344888,   0.191022,   0.173440,   0.337780,   0.756612,   0.213484,   0.298083,   0.062904,  -0.164975,  -0.354983,  -0.202565,  -0.263482,   0.114889,   0.332008,   0.099940,  -0.095354,  -0.009766,   0.169388,   0.206749,  -0.111808,  -0.288887,   0.129370,  -0.117981,  -0.106499,  -0.266949,  -0.046878,  -0.150885,
        -0.118587,  -0.009896,  -0.086578,   0.130258,   0.084734,  -0.039420,   0.234204,   0.406766,   0.476533,   0.509599,  -0.169659,   0.216232,  -0.065689,   0.549118,  -0.070943,   0.087047,   0.490984,  -0.422113,  -0.180315,  -0.173134,  -0.228394,   0.036646,  -0.112598,  -0.779985,   0.185928,   0.273930,   0.058692,  -0.144147,  -0.029376,   0.243415,   0.288349,   0.313784,
         0.654912,   0.050271,  -0.052386,   0.293530,   0.085631,   0.016827,   0.641948,  -0.120567,  -0.531691,  -0.321783,  -0.454751,  -0.155128,  -0.344870,  -0.442173,  -0.122703,  -0.011277,   0.478001,  -0.105196,  -0.621367,   0.051719,   0.326735,  -0.977598,  -0.235963,  -0.197529,   0.337496,   0.143767,  -0.122240,   0.130382,  -0.151092,   0.122121,   0.308612,   0.120933,
        -0.038011,   0.197207,   0.019844,   0.518090,   0.575465,   0.205046,  -0.024372,  -0.333282,  -0.273242,   0.141886,  -0.228005,  -0.265232,   0.199985,   0.219500,   0.346435,   0.061049,  -0.289425,   0.159937,  -0.008618,  -0.895822,   0.245320,   0.080692,  -0.402834,  -0.220966,   0.150254,  -0.027955,  -0.137935,   0.174240,  -0.342728,   0.261925,  -0.094680,   0.094293,
         0.132275,   0.313540,   0.453221,  -0.153065,   0.424720,   0.128915,  -0.432206,  -0.094454,   0.144793,   0.512230,   0.003422,  -0.014260,  -0.329406,  -0.271373,   0.116721,  -0.186755,  -0.270528,  -1.119480,  -0.023132,  -0.067423,  -0.275236,  -0.002308,  -0.269234,  -0.076069,  -0.022410,   0.103847,   0.302314,   0.033703,  -0.009011,   0.112916,   0.231766,   0.249918,
         0.218129,   0.302854,  -0.222046,  -0.104577,  -0.086532,   0.068512,  -0.042543,  -0.128117,   0.464187,   0.298185,  -0.404418,  -0.131768,  -0.148274,  -0.198955,   0.159237,  -0.087858,  -0.049420,  -0.171090,   0.182913,   0.305219,  -0.237672,   0.696507,   0.547869,   0.477914,   0.661873,   0.077758,  -0.033811,   0.388680,   0.345189,   0.182266,  -0.142175,   0.057572,
        -0.275037,  -0.311863,  -0.212190,   0.091814,  -0.272670,  -0.491110,  -0.043897,   0.098203,   0.421509,   0.392666,  -0.197738,  -0.041770,   0.229539,  -0.559840,  -0.242821,  -0.359473,   0.443005,   0.153426,   0.333681,   0.224086,  -0.167565,   0.172010,   0.272319,  -0.012205,   0.196460,   0.177267,   0.530670,   0.059850,  -0.365005,   0.663154,   0.092292,  -0.488004,
        -0.064614,   0.106274,  -0.161335,   0.077606,  -0.167381,  -0.136192,   0.234717,   0.215071,  -0.213360,   0.139752,   0.086535,  -0.661387,  -0.146277,  -0.289719,   0.061066,   0.312204,   0.072416,   0.253692,  -0.099007,   0.158297,   0.369995,   0.246487,   0.068100,   0.411967,   0.594613,   0.260168,  -0.161558,   0.291889,  -0.257346,  -0.210128,  -0.449290,   0.601151,
        -0.427452,  -0.099851,   0.214816,   0.166730,   0.106820,  -0.009960,  -0.266142,  -0.097163,  -0.211419,  -0.446301,  -0.060226,   0.039788,  -0.010464,   0.103810,   0.015814,   0.111186,  -0.428451,   0.034062,   0.553720,  -0.088122,   0.019406,   0.200670,   0.006075,   0.896145,  -0.927244,  -0.008801,  -0.456926,  -0.191052,  -0.890313,   0.013570,  -0.246357,   0.139814,
         0.105891,  -0.275191,   0.323740,  -0.390059,  -0.691497,  -0.183471,   0.138042,   0.430859,   0.393569,   0.295931,  -0.325566,   0.070891,   0.075351,   0.092949,   0.050606,   0.115343,   0.087848,   0.252470,   0.043274,   0.068497,   0.385337,   0.048883,   0.403882,  -0.599694,  -0.241273,   0.261991,   0.825384,  -0.314302,  -0.014897,  -0.293020,  -0.089213,   0.241494,
        -0.053889,  -0.053750,  -0.028313,  -0.122579,   0.159876,  -0.366701,   0.092789,   0.021421,  -0.873173,  -0.313339,   0.201776,   0.087935,   0.229674,   0.425983,  -0.041151,   0.184080,   0.042065,  -0.038219,   0.271644,   0.427878,  -0.548668,  -0.182358,   0.109458,   0.561349,   0.119476,   0.037959,   0.213601,   0.381325,   0.415813,  -0.120288,  -0.245932,  -0.277853,
        -0.244822,  -0.019546,  -0.286228,  -0.496843,   0.450057,   0.275202,  -0.324631,   0.002534,  -0.305073,  -0.093240,  -0.039653,   0.161204,   0.054852,   0.266418,   0.134433,   0.021564,   0.456958,   0.224697,   0.383803,  -0.144931,   0.103270,   0.010302,   0.591615,   0.044660,  -0.029003,   0.205843,  -0.303449,   0.363014,   0.006074,   0.052915,   0.317358,  -0.021953,
         0.293342,   0.295345,   0.178897,   0.252466,   0.173851,   0.316124,   0.137781,  -0.161349,  -0.081071,   0.040322,  -0.036951,   0.211623,  -0.094141,   0.018904,  -0.462953,   0.208459,   0.543926,  -0.132407,  -0.045069,  -0.115847,  -0.159609,  -0.122297,   0.037607,   0.067524,   0.203769,   0.807174,  -0.299178,  -0.646677,  -0.053664,  -0.244517,  -0.371751,   0.112897,
         0.183937,   0.226448,   0.087817,  -0.069162,  -0.327513,   0.422165,   1.397430,   0.739720,
        -0.000381,  -0.164204,  -0.089905,  -0.105782,   0.173075,  -0.038118,   0.048986,  -0.053297,  -0.290274,   0.167400,  -0.120456,  -0.032094,   0.285521,   0.116362,   0.282855,  -0.124066,   0.101562,   0.240403,  -0.015441,   0.138594,  -0.121582,   0.006295,   0.294600,   0.175539,   0.049183,   0.491164,   0.133751,  -0.033613,   0.109758,  -0.125967,  -0.249235,   0.104569,
        -0.076922,  -0.164371,  -0.061135,  -0.122121,   0.034139,  -0.046695,   0.007199,   0.379694,  -0.640479,  -0.112728,   0.067915,  -0.027140,   0.290560,   0.254052,  -0.109593,  -0.022160,   0.365461,  -0.078972,  -0.168866,  -0.127767,  -0.393554,  -0.232205,   0.159081,  -0.110926,   0.101780,   0.048478,  -0.294355,  -0.174919,  -0.030709,   0.066727,  -0.060868,  -0.113082,
        -0.158654,  -0.031202,   0.206062,  -0.015024,  -0.023005,   0.866857,  -0.401216,  -0.052412,   0.143313,   0.121795,   0.389986,   0.247104,   0.111594,   0.078906,   0.054239,  -0.092644,  -0.081536,   0.189591,   0.269054,  -0.022233,  -0.437326,  -0.149619,   0.112781,  -0.031356,  -0.253664,  -0.002359,   0.173796,  -0.008932,   0.091481,  -0.025608,  -0.051415,   0.067809,
         0.183176,  -0.038427,  -0.035930,   0.931956,  -0.583209,   0.290984,   0.333517,  -0.019131,  -0.053307,   0.128142,   0.173376,   0.542155,   0.332883,  -0.079758,  -0.046263,   0.073708,   0.510177,   0.050223,  -0.208769,   0.236313,   0.111680,   0.024475,  -0.143501,   0.094665,   0.107217,  -0.034524,   0.079348,  -0.184514,   0.192310,  -0.058615,   0.037489,   0.112879,
        -0.191726,   1.913630,  -0.018767,  -0.201560,  -0.592462,  -0.094799,  -0.108776,   0.288958,  -0.088615,  -0.020641,   0.010920,  -0.285384,  -0.143553,  -0.302364,  -0.347208,   0.140692,   0.140607,  -0.337522,   0.423999,   0.268900,   0.253964,   0.103836,  -0.106599,  -0.140114,   0.010647,  -0.114152,  -0.054368,  -0.055126,   0.003382,   0.074650,  -0.165022,   1.620370,
         0.434562,   0.066265,  -0.200716,   0.044811,  -0.018733,   0.005732,  -0.121628,  -0.242561,  -0.103823,  -0.143950,  -0.311844,   0.144423,  -0.333278,   0.327984,   0.333738,  -0.428609,   0.146584,  -0.182082,  -0.354733,  -0.018493,   0.154084,  -0.171076,   0.149995,   0.047567,  -0.036513,   0.398484,   0.061814,  -0.084902,  -0.016826,   1.602540,   0.109152,  -0.153716,
         0.160265,   0.003621,   0.040257,  -0.238755,   0.046169,   0.035060,   0.243092,  -0.128248,   0.241337,  -0.016253,   0.159947,  -0.094954,  -0.063812,  -0.214970,   0.031768,   0.045240,   0.050667,   0.060181,   0.361475,  -0.093697,   0.083036,  -0.042455,   0.140951,   0.187533,  -0.067996,   0.036475,  -0.443709,   0.884813,  -0.194712,   0.035280,  -0.267773,  -0.293734,
        -0.061260,  -0.189087,   0.014291,  -0.184511,   0.179427,   0.014996,  -0.163478,  -0.196202,  -0.353539,  -0.063966,  -0.182418,   0.299374,   0.277649,  -0.106101,  -0.175668,  -0.108482,  -0.026245,   0.065691,  -0.180033,  -0.025142,  -0.124702,  -0.091333,   0.008113,  -0.129706,  -0.520505,   0.243489,  -0.280544,   0.087251,  -0.061138,  -0.050268,  -0.328289,  -0.028414,
         0.288590,  -0.025758,  -0.125209,  -0.068660,  -0.102098,  -0.009414,  -0.163201,   0.332149,   0.074428,   0.203553,   0.302688,   0.271651,   0.037459,   0.065543,   0.018506,   0.047856,  -0.091553,  -0.311783,   0.300539,  -0.140647,   0.093166,   0.032010,  -0.195233,  -0.693818,   0.560504,   0.354350,  -0.573794,  -0.505023,  -0.003139,   0.497563,  -0.332816,  -0.641308,
        -0.095407,  -0.064094,   0.317788,   0.544994,  -0.374234,   0.086103,   0.428155,  -0.545667,  -0.261267,  -0.602443,   0.252457,  -0.072111,   0.595256,   0.274996,   0.201584,   0.114470,   0.016723,   0.046587,   0.089930,   0.239134,  -0.178765,  -0.027370,  -0.317189,   0.462611,   0.865341,   0.382303,   0.171773,   0.164311,   0.437724,   0.772434,   0.386280,   0.567111,
         0.067195,  -0.140457,   0.417296,   0.275876,  -0.089904,  -0.744346,   0.347255,   0.232557,   0.155900,   0.601565,  -0.060870,  -0.047427,  -0.384077,  -0.535409,   0.349938,   0.014594,   0.247602,   0.149703,  -0.428949,  -0.116878,  -0.724223,   0.533796,  -0.882710,  -0.706724,  -0.558439,  -0.125539,  -0.214674,  -0.766759,   0.204331,  -0.243632,   0.025052,   0.094273,
        -0.658517,  -0.137013,   0.527325,  -0.149790,  -0.310475,  -0.525862,   0.203811,   0.433405,   0.584112,   0.031375,  -0.152243,  -0.211180,  -0.135069,  -0.093050,  -0.124900,   0.134764,  -1.021580,  -0.084261,  -0.736002,  -0.031906,   0.273904,   0.007921,   0.198779,  -0.003805,   0.380597,   0.107689,   0.188133,   0.326278,   0.035535,  -0.089913,   0.062933,  -0.107620,
         0.096523,  -1.105440,   0.372608,   0.405218,   0.219408,  -0.004834,  -0.073636,  -0.104714,  -0.127537,  -0.007610,   0.382688,  -0.035823,   0.157018,  -0.004424,  -0.586816,   0.135641,  -0.325071,   0.184180,  -0.088044,   0.177043,  -0.388280,   0.188880,  -0.260830,  -0.386581,  -0.239407,  -0.033243,   0.146722,   0.122341,  -0.474379,  -0.106343,   0.269829,  -0.072060,
         0.125895,  -0.383226,  -0.021600,  -0.067959,   0.152492,   1.772830,   1.679990,  -0.519395,
         0.152629,   0.036927,   0.137427,   0.137501,   0.036477,   0.151359,   0.142229,   0.069194,   0.196460,   0.052007,  -0.522439,   0.029237,  -0.098588,  -0.011130,   0.106375,   0.141399,   0.135806,   0.195474,  -0.005007,  -0.057772,  -0.014361,   0.031995,   0.147005,   0.165382,  -0.036837,   0.368810,  -0.081300,  -0.011288,   0.065510,   0.060205,   0.051571,   0.028090,
         0.016645,   0.134948,   0.186990,   0.018798,   0.120934,   0.164662,   0.252755,   0.403908,  -0.598388,   0.088635,  -0.197851,  -0.159410,  -0.015429,   0.343094,  -0.542746,   0.002000,   0.203529,   0.302920,  -0.178060,  -0.097238,  -0.307031,  -0.060989,   0.105470,  -0.302169,   0.093188,   0.066350,   0.177575,   0.324739,  -0.098971,  -0.091855,   0.092013,   0.022123,
         0.046189,   0.164221,   0.243213,   0.010620,   0.078301,   0.255649,  -0.127072,   0.136288,  -0.033799,  -0.129120,  -0.140645,  -0.351005,   0.007709,  -0.075074,  -0.219346,   0.080883,  -0.160315,   0.023896,  -0.237958,  -0.231374,  -0.096348,  -0.392129,  -0.133729,   0.234103,   0.270979,  -0.135380,  -0.026604,   0.030318,   0.093347,   0.271886,  -0.141013,   0.226496,
         0.092539,   0.093394,   0.035629,   0.196605,   0.093227,   0.004935,   0.435007,   0.025434,   0.104196,  -0.196123,   0.299664,   0.085159,  -0.075945,   0.115021,   0.032215,   0.011714,   0.369250,   0.055118,  -0.196353,  -0.056282,   0.291955,   0.111197,   0.111685,   0.024715,   0.000307,   0.160046,   0.176641,   0.163567,  -0.005952,  -0.143378,   0.064389,   0.053468,
         0.021371,   0.351162,   0.253262,   0.098632,  -0.045535,  -0.018415,   0.304326,  -0.033670,  -0.049624,  -0.019444,  -0.036449,  -0.117934,   0.149473,  -0.057720,  -0.173588,  -0.041317,   0.210942,   0.018652,  -0.028417,   0.179348,   0.163119,  -0.002440,  -0.070525,   0.284569,   0.111159,   0.148458,   0.397449,   0.024909,   0.071090,   0.230898,  -0.024451,   0.405276,
         0.472343,   0.196337,  -0.579574,  -0.357157,  -0.482658,   0.210674,  -0.283922,  -0.071028,  -0.132387,  -0.111770,   0.001506,   0.142040,  -0.381605,  -0.211093,   0.099651,  -0.391796,  -0.097650,   0.284450,  -0.558299,   0.174843,   0.100208,  -0.418539,  -0.153807,   0.099935,   0.146339,   0.215671,  -0.079495,  -0.253363,   0.005972,   0.322652,   0.455436,  -0.028702,
        -0.006218,  -0.139259,  -0.091339,   0.142483,  -0.083455,  -0.108356,   0.653226,  -0.239701,  -0.075781,   0.295818,   0.018071,   0.001171,  -0.096217,  -0.523571,  -0.130567,  -0.180663,   0.128168,   0.020674,  -0.388509,   0.014609,   0.041121,   0.336267,  -0.214051,  -0.004024,   0.102762,   0.051066,   0.111377,   0.959621,  -0.008334,  -0.238279,   0.095142,   0.243572,
         0.470966,  -0.138674,   0.438942,   0.113637,  -0.307263,  -0.061861,  -0.115789,  -0.055672,   0.445266,  -0.133703,  -0.422328,   0.316230,  -0.068159,  -0.161010,  -0.580857,   0.052160,  -0.172906,   0.100012,  -0.172620,   0.046512,   0.316801,  -0.117426,   0.295815,   0.146254,  -0.124052,  -0.009868,   0.206326,   0.358269,  -0.019649,  -0.218258,   0.352258,   0.058672,
         0.055003,   0.219430,  -0.459961,   0.038099,  -0.114302,   0.222496,   0.001247,   0.220003,  -0.047259,   0.128998,   0.306006,   0.212391,   0.101334,  -0.066744,  -0.006276,   0.002940,  -0.319858,   0.086126,   0.000799,   0.147792,   0.125876,   0.110643,  -0.044654,   0.725858,  -1.436800,  -0.136508,   0.262469,   0.404027,  -0.208813,  -0.046656,   0.272763,   0.685021,
        -0.104772,   0.335366,   0.025778,   0.091971,  -0.190832,   0.240167,  -0.350499,   0.701086,   0.268625,   0.287586,   0.131524,   0.206202,  -0.351530,  -0.163094,   0.256434,  -0.079746,   0.195931,  -0.259026,   0.206020,   0.035158,   0.190256,   0.412698,  -0.216456,  -0.298454,  -0.341950,  -0.029860,  -0.205594,   0.165205,  -0.143651,  -0.622508,   0.172665,  -0.295112,
        -0.272393,  -0.161409,  -1.116720,   0.239556,  -0.027817,  -1.184190,   0.204620,   0.241873,  -0.175656,  -0.413424,   0.326260,  -0.020385,   0.238515,   0.311295,   0.071267,  -0.189947,   0.165533,  -0.066567,   0.244947,   0.075032,   0.280472,   0.096427,  -0.534689,   0.150972,   0.522540,  -0.158211,   0.113559,  -0.475522,   0.152882,  -0.172208,  -0.404648,  -0.045664,
        -0.797840,  -0.080610,   0.003464,  -0.612462,  -0.076086,  -0.535436,   0.082960,  -0.103166,  -0.633356,   0.159679,  -0.143963,   0.145267,   0.188661,  -0.029209,   0.123172,   0.227112,   0.014816,   0.089940,   0.403689,   0.010574,  -0.012247,   0.582769,   0.410123,   0.047156,  -0.001053,   0.022514,   0.091894,  -0.040019,  -0.562602,  -0.084888,  -0.185963,  -0.065309,
        -0.215718,   0.576995,   0.115137,   0.069112,   0.266055,  -0.185266,  -0.307024,   0.128019,   0.214994,   0.349678,   0.056172,   0.020212,  -0.053345,   0.074639,   0.695510,  -0.396794,   0.632395,  -0.343376,   0.038841,   0.093066,  -0.592773,  -0.135693,  -0.035374,  -0.418306,  -0.190781,  -0.205924,  -0.081997,   0.007944,   0.023477,   0.125024,   0.091741,  -0.242319,
        -0.144365,  -0.002658,   0.028978,   0.018729,  -0.403698,   0.400428, 114.878000,  -1.613180,
         0.266416,   0.212652,  -0.012796,   0.189770,   0.116980,  -0.026259,  -0.308267,   0.208427,   0.072611,  -0.287970,   0.837999,   0.035710,  -0.564521,  -0.030475,  -0.429168,  -0.336319,  -0.237257,   0.251613,   0.519548,   0.071354,  -0.194962,   0.037264,   0.015121,  -0.065225,  -0.093548,   0.305507,   0.312643,   0.002778,   0.577498,  -0.174050,   0.033200,   0.037581,
         0.015833,   0.055686,   0.015687,  -0.103591,  -0.000512,   0.166838,   0.027352,  -0.182774,   0.750173,  -0.074810,   0.363286,   0.198815,   0.267908,  -0.131570,  -0.023652,   0.173986,  -0.226085,   0.082572,  -0.252578,  -0.295803,  -0.027600,   0.102213,   0.089475,  -0.560008,  -0.052019,   0.374262,   0.558306,   0.050144,  -0.204512,   0.208550,  -0.030461,   0.237356,
         0.164802,  -0.016550,  -0.263514,   0.339158,   0.082275,  -0.455410,   0.931880,   0.017531,  -0.039598,  -0.008397,  -0.436165,   0.062071,   0.002077,  -0.349702,  -0.047790,  -0.186938,   0.178503,  -0.503145,   0.108238,   0.024933,  -0.402523,  -0.351257,  -0.052869,  -0.115202,  -0.429298,  -0.413630,  -0.243642,   0.002465,   0.096276,   0.193952,   0.262752,   0.018485,
         0.033181,   0.069777,  -0.128857,   0.876490,   0.874444,  -0.004455,  -0.120813,  -0.441996,  -0.438713,   0.114530,  -0.270016,  -0.198672,   0.041618,  -0.321220,   0.212312,   0.161505,  -0.520305,  -0.082320,  -0.296248,  -0.277536,   0.070561,  -0.134184,  -0.182938,  -0.324253,   0.238590,  -0.046993,   0.101106,   0.347148,  -0.376920,   0.305683,  -0.047184,   0.060688,
        -0.012409,  -0.112544,   0.750148,  -0.156681,   0.749714,   0.517639,  -0.409622,  -0.036632,   0.283471,   0.561906,   0.157449,   0.094936,  -0.151123,  -0.256358,   0.453355,  -0.089490,  -0.217356,  -0.709777,   0.076677,   0.163982,  -0.284356,   0.131629,  -0.325211,   0.204611,   0.005032,   0.109858,   0.110978,   0.126316,  -0.031940,   0.173065,  -0.108530,   0.025832,
         0.751187,   0.149156,  -0.046286,   0.028087,   0.560935,  -0.194412,   0.190882,   0.455066,   0.336029,   0.238950,  -0.364017,  -0.124450,  -0.126248,  -0.314200,  -0.178908,  -0.432115,   0.116067,   0.230727,   0.196432,  -0.081026,  -0.236485,  -0.164851,  -0.196554,   0.185509,   0.180998,   0.184249,  -0.212046,  -0.006255,   0.102579,   0.667964,   0.267487,   0.227582,
        -0.339694,  -0.179260,  -0.114864,  -0.268775,  -0.213052,  -0.374178,   0.258024,   0.183481,  -0.271974,  -0.456721,  -0.240713,  -0.298395,  -0.708535,  -0.221873,   0.029392,  -0.227761,  -0.268668,  -0.191619,  -0.254559,   0.147532,  -0.003833,   0.130109,   0.013252,   0.382519,  -0.007863,   0.303483,   0.048004,   0.522259,  -0.018463,   0.335117,   0.681586,  -0.466539,
        -0.365403,   0.116672,   0.050946,   0.219666,  -0.099778,  -0.383243,   0.165148,   0.191257,   0.245985,   0.144120,  -0.446913,  -0.209664,  -0.303447,  -0.257872,   0.266233,   0.137466,  -0.108897,  -0.033113,  -0.045429,   0.126990,  -0.077771,   0.242056,   0.147325,   0.030937,  -0.318461,   0.746900,   0.420176,  -0.175824,   0.385960,   0.061909,   0.106889,  -0.097515,
         0.024357,   0.454617,   0.480428,   0.502462,  -0.223823,   0.004946,   0.614598,  -0.169851,   0.255213,  -0.068206,  -0.110487,  -0.388680,   0.567608,   0.477443,  -0.318400,  -0.090991,  -0.523666,  -0.177967,  -0.042925,  -0.059840,  -0.061469,   0.188484,   1.118070,   0.432596,  -0.996589,  -0.083989,  -0.054177,   0.137828,  -1.010750,  -0.075919,  -0.003582,   0.180268,
         0.278503,  -0.585438,   0.095654,  -0.430799,  -0.192686,  -0.239293,   0.273179,   0.414322,   0.281819,   0.313475,  -0.516393,   0.041231,  -0.160998,  -0.042568,  -0.125764,  -0.184144,   0.030698,   0.134810,   0.131905,   0.022682,   0.349830,  -0.080436,   0.669923,  -0.126838,   0.165887,   0.182120,   0.659456,   0.196729,   0.139499,  -0.446974,  -0.194984,   0.483793,
         0.006258,   0.257037,   0.178041,   0.474594,   0.160900,  -0.589155,   0.090943,   0.066548,  -0.145099,  -0.420744,   0.354737,  -0.042239,   0.119500,   0.208510,  -0.127012,   0.220277,   0.199724,  -0.098645,   0.517118,   0.246366,  -0.423271,  -0.370466,   0.669811,   1.730340,   0.836857,  -0.232120,   1.050370,   0.436081,   0.761116,   0.364197,  -0.445627,  -0.492744,
         0.716270,   0.282752,  -0.259881,  -0.176467,   0.590073,   0.348432,   0.579890,   0.145842,  -0.645794,  -0.226339,  -0.277146,  -0.102573,  -0.134448,   0.056653,   0.284989,  -0.057574,   0.624893,   0.058960,   0.402454,   0.448059,   0.881456,  -0.104166,   0.346085,  -0.235949,  -0.045227,   0.098030,   0.094171,   0.322550,   0.027480,   0.165335,   0.156771,  -0.117370,
         0.671834,   0.018801,   0.183229,   0.387063,  -0.079463,  -0.097695,   0.661005,   0.353370,   0.382970,   0.322248,  -0.233528,   0.397854,   0.090327,   0.189051,  -0.691684,   0.152653,   0.547675,  -0.405107,   0.815007,  -0.186171,  -0.159611,  -0.214979,  -0.009266,   0.321979,  -0.842669,   0.640562,  -0.032743,  -0.403376,   0.079489,  -0.139825,  -0.205092,   0.404782,
         0.220068,   0.273519,  -0.318513,  -0.181951,  -0.331097,   0.522530,   0.105415,   1.178960,
         0.064164,   0.005437,  -0.125440,   0.009121,  -0.125244,   0.178260,  -0.339396,   0.006083,  -0.108163,  -1.367540,   0.375562,  -0.135934,   0.039056,   0.398993,  -0.174075,   0.341041,   0.287839,  -0.231490,   0.186377,  -0.459233,   0.063998,   0.353459,   0.540589,   0.302794,   0.157976,  -0.503637,   0.151568,   0.043138,   0.421945,   0.428429,  -0.163402,   0.118966,
         0.177565,   0.127207,  -0.245956,  -0.064326,   0.039407,  -0.043612,  -0.124501,  -0.736727,   0.290811,  -0.344576,  -0.105373,  -0.048599,   0.061885,   0.256219,   0.148868,  -0.078697,  -0.486060,  -0.203230,   0.053987,  -0.191011,  -0.284146,  -0.337912,   0.152084,   0.046948,   0.188780,   0.189505,  -0.451767,   0.161634,  -0.250420,   0.103899,  -0.115856,  -0.130493,
        -0.149195,   0.020530,   0.306779,   0.158287,  -0.070888,   0.114006,   0.043391,   0.018725,   0.127707,  -0.282350,  -0.078277,  -0.825586,   0.002214,  -0.118702,   0.065337,  -0.038358,  -0.168166,   0.082285,   0.131526,  -0.307190,   0.229278,   0.153284,  -0.152929,  -0.218295,  -0.115894,   0.319507,  -0.234763,   0.012872,  -0.131263,  -0.216550,  -0.223593,   0.078518,
         0.160119,   0.126878,  -0.172294,   0.718696,   0.119922,  -0.076110,   0.193610,   0.216926,  -0.023233,   0.134712,   0.074960,  -0.036158,   0.298713,   0.154633,   0.121914,   0.429702,   0.932632,  -0.033996,   0.251258,   0.629348,  -0.103254,  -0.146924,   0.474190,   0.151027,   0.057343,  -0.164995,   0.099337,  -0.021586,   0.085974,   0.484720,  -0.068570,  -0.096368,
        -0.120097,   1.996300,   1.578650,  -0.466340,   0.079270,   0.178910,   0.121429,   0.028983,   0.225877,  -0.315540,   0.463603,  -0.213901,  -0.307843,  -0.006713,   0.040573,   0.261734,   0.205896,   0.242889,   0.098465,  -0.256667,   0.194851,   0.064160,  -0.163203,   0.016486,  -0.290525,  -0.197095,   0.098259,   0.314736,  -0.161495,   0.028939,  -0.057779,   1.906650,
         0.972200,  -0.210965,  -0.077161,   0.101717,   0.151347,  -0.221166,  -0.158916,  -0.284874,  -0.377728,   0.080612,  -0.518639,  -0.218112,  -0.451874,  -0.085774,   0.108043,   0.184499,   0.077752,   0.043831,  -0.352202,  -0.079336,  -0.112575,   0.121682,  -0.032157,  -0.116844,   0.161958,   0.142790,   0.201290,   0.119221,   0.068534,   1.787820,   0.322876,  -0.095433,
         0.288937,  -0.164372,  -0.127432,  -0.332803,  -0.066711,  -0.122966,   0.018929,   0.020087,  -0.166997,   0.103158,   0.009540,  -0.160149,   0.163218,   0.069287,  -0.074907,  -0.518856,  -0.095810,  -0.084519,   0.085976,  -0.001493,   0.424051,   0.179711,   0.347479,   0.100102,  -0.131212,   0.082369,  -0.239542,   1.320290,  -0.199160,  -0.433152,   0.120278,  -0.179198,
         0.162375,  -0.087764,   0.002193,   0.584151,   0.494380,   0.150784,  -0.087031,  -0.061336,   0.175262,   0.175777,   0.015910,   0.426291,   0.028623,  -0.075563,   0.298014,  -0.177303,   0.261974,  -0.127296,  -0.057264,   0.012580,  -0.323725,  -0.179563,  -0.102771,  -0.209025,  -0.848546,  -0.238689,  -0.998652,  -0.264157,   0.113371,   0.000464,  -0.228030,  -0.162310,
         0.191299,   0.275390,  -0.092166,  -0.051362,  -0.278055,  -0.178146,  -0.544566,   0.104063,  -0.077826,  -0.004260,   0.100146,   0.205460,  -0.230614,   0.237945,  -0.014840,   0.239946,   0.196176,  -0.084009,   0.216709,   0.339977,   0.023366,   0.115795,  -0.630596,  -0.688967,   1.361530,   0.282380,   0.035199,  -0.083672,  -0.362315,   0.103393,  -0.038592,   0.126410,
        -0.399948,   0.274169,   0.401466,   0.193586,   0.441271,  -0.584485,   0.497815,  -0.087605,  -0.352049,  -0.186492,   0.192455,   0.149776,   0.477872,   0.292387,  -0.102144,  -0.005565,   0.360140,  -0.185961,   0.420008,   0.365194,  -0.944350,   1.003440,  -0.478746,   0.426539,   0.325188,  -0.698370,  -1.063480,   0.384925,   0.319523,   0.520858,  -0.017531,  -0.188084,
         0.823000,   0.215296,   0.426312,   0.268630,   0.294754,  -0.621486,   0.007757,  -0.109141,  -0.022508,   0.192073,   0.560781,   0.092360,   0.601438,  -0.042698,   0.257100,   0.317750,   0.643879,   0.391815,  -1.255460,   0.793920,  -0.799996,   0.752827,   0.361200,  -0.622648,  -0.457237,   0.585981,  -0.037095,   0.301268,  -0.651656,  -0.120271,   0.878953,  -0.263884,
        -0.242242,  -0.162245,   0.252319,  -0.709098,  -0.025307,  -0.520853,  -0.109552,   0.340187,   1.256460,   0.410364,  -0.049001,  -0.266057,   0.147863,   0.078202,   0.052722,   0.496724,  -2.190130,   1.954530,  -1.401800,   0.273714,  -0.414941,  -0.495575,  -0.589330,   0.480849,  -0.018398,  -0.219245,  -0.360752,   0.060149,   0.361337,  -0.457669,  -0.138730,   0.030430,
        -0.069285,  -0.994234,   0.536136,   0.015351,   0.485321,   0.208104,   0.419165,  -0.021207,   0.299262,   0.282392,   0.136648,   0.123940,   0.153134,   0.023879,  -1.898160,   2.115450,  -0.752008,   0.087388,  -0.079341,   0.324628,  -0.780623,   0.426908,  -0.086168,  -0.484683,  -0.198408,  -0.637200,   0.213690,   0.325301,  -0.745916,   0.117574,   0.335877,  -0.404622,
         0.536026,  -0.355016,  -0.398117,  -0.243850,   0.213250,   1.427300,   0.179184,  -0.057012,
        -0.139401,  -0.084112,   0.079519,  -0.106381,   0.167346,  -0.233101,  -0.062392,  -0.026706,  -0.328454,   0.141427,   0.040466,  -0.057448,  -0.223185,   0.366555,   0.467295,  -0.533772,  -0.385565,   0.199288,   0.422601,  -0.120649,   0.193862,  -0.094606,  -0.183804,  -0.151005,   0.038665,   0.503944,   0.073553,   0.258780,   0.783520,  -0.109441,   0.307229,  -0.202638,
         0.043043,  -0.115917,  -0.290277,  -0.086674,   0.046675,  -0.114823,  -0.386664,   0.323241,   0.138021,  -0.245840,   0.444034,   0.287799,   0.062761,   0.174053,  -0.076999,   0.009109,   0.252220,  -0.006791,  -0.096055,  -0.426343,   0.094489,  -0.226982,  -0.147709,  -0.332676,   0.203885,   0.232120,   0.564153,   0.317791,   0.216284,   0.143219,  -0.216232,  -0.415904,
        -0.215402,  -0.020035,  -0.222838,   0.240539,  -0.241872,   0.368606,   0.045885,  -0.296968,  -0.074332,   0.079011,  -0.179868,  -0.318230,   0.126413,  -0.578992,   0.167002,   0.111565,  -0.287740,   0.137149,   0.089085,   0.072841,   0.056931,   0.060473,   0.236582,   0.347212,   0.637127,   0.128469,   0.307258,   0.051004,   0.113003,  -0.104292,  -0.086741,  -0.066252,
         0.147906,  -0.019182,   0.400936,  -2.005490,  -0.584256,  -0.312725,   0.044417,   0.245236,   0.199999,   0.446656,  -0.164535,  -0.250302,  -0.180434,  -0.383286,   0.081780,  -0.068026,  -0.050838,  -0.348415,  -0.074052,  -0.386458,  -0.009196,   0.013523,  -0.167444,   0.211447,  -0.090322,  -0.078304,   0.388024,  -0.056835,   0.386378,   0.133511,   0.194017,   0.013560,
         0.668522,  -1.154150,  -1.393670,   0.206012,   0.589881,   0.165047,   0.546128,  -0.120493,   0.279397,   0.235698,   0.229350,   0.015579,  -0.272015,   0.043738,   0.413420,   0.202556,  -0.135640,   0.822235,   0.403960,  -0.087147,   0.010165,   0.177588,  -0.281203,  -0.234170,   0.506578,  -0.181785,  -0.288745,  -0.229502,   0.493111,  -0.081905,   0.626576,  -1.516220,
        -1.416570,   0.233295,   0.419983,   0.372671,  -0.118554,  -0.099383,   0.069365,   0.492579,   0.030772,   0.647594,   0.044102,  -0.298573,   0.175900,   0.335062,   0.059524,  -0.329769,   0.080890,   0.267042,   0.653324,   0.088551,  -0.091668,  -0.683507,   0.145457,  -0.283009,  -0.216610,   0.035734,   0.066336,  -0.268469,   0.953603,  -1.261070,  -1.326170,   0.004845,
         0.015068,   0.097262,  -0.332516,   0.578528,  -0.156592,  -0.293919,   0.148340,  -0.248937,  -0.052293,   0.386812,  -0.054865,  -0.099659,  -0.102574,  -0.733764,   0.272605,  -0.028693,   0.409235,  -0.246226,   0.435308,  -0.482872,  -0.081151,  -0.388856,  -0.166157,  -0.053727,  -0.312873,   0.033788,   0.458166,  -1.843730,  -0.371377,  -0.421603,   0.290232,   0.032049,
         0.064679,   0.504900,   0.156526,  -0.207556,   0.504784,  -0.026562,   0.083222,   0.238752,   0.339754,   0.192627,  -0.427158,  -0.398293,  -0.002614,  -0.068541,  -0.132768,  -0.075729,   0.312639,  -0.052968,   0.063545,   0.116425,  -0.106982,  -0.236158,  -0.060208,  -0.091674,  -0.045466,  -1.429040,   0.139544,  -0.046281,  -0.015582,   0.509463,   0.500572,   0.014150,
         0.881471,   0.665614,   0.573470,   0.771576,   0.035114,  -0.031908,   0.331497,   0.516084,   0.070358,  -0.131424,   0.227215,   0.101556,   0.250054,   0.400890,  -0.219357,  -0.047404,   0.459231,  -0.163512,   0.005244,   0.083430,   0.208325,   0.106018,  -0.116505,  -0.294565,   0.393573,  -0.083992,   0.332166,  -0.144988,  -1.133160,   0.598220,   0.104439,  -0.635710,
         0.056020,  -0.763491,   0.673633,   0.137904,   0.059949,   0.322724,   0.320355,  -0.675526,   0.083567,  -0.289843,  -0.583596,  -0.298137,  -0.047796,  -0.120651,  -0.292485,  -0.102421,   0.018106,  -0.260536,   0.365589,  -0.168072,  -0.200040,   0.441711,  -0.795647,   0.496913,  -0.058410,  -0.098327,  -1.624360,   0.859509,   0.120574,  -0.324994,  -0.253700,  -0.017348,
         0.761300,   0.268919,  -0.343914,   0.161521,   0.145240,  -0.501350,  -0.100540,  -0.000774,   0.874905,  -0.169179,   0.002336,   0.504802,   0.231648,  -0.020395,   0.362322,   0.077634,   0.290560,   0.405150,  -0.553898,   0.346822,  -0.658231,   0.464063,  -0.233566,  -0.346879,  -1.252220,   0.802210,   0.101515,  -0.405681,  -1.251500,  -0.034708,   0.565491,  -0.093634,
        -0.194372,   0.385593,   0.134862,  -0.328999,  -0.286860,   0.169577,  -0.133964,  -0.096680,   0.404803,   0.171778,  -0.122884,  -0.071970,  -0.155699,  -0.300937,   0.100053,   0.220519,  -1.682420,  -0.068249,  -0.394885,   0.043377,  -0.030591,   0.016705,  -0.052691,   0.818546,   0.313902,  -0.003314,  -0.456166,   0.164703,   0.209667,  -0.274121,  -0.155084,   0.175348,
         0.337748,  -0.517040,   0.261000,   0.475299,   0.661113,   0.377933,   0.086993,   0.323542,   0.335659,   0.162593,   0.398835,   0.076497,   0.275108,   0.212842,  -1.512440,  -0.314095,  -0.014716,   0.037259,  -0.057054,   0.460759,  -0.019717,   0.412257,  -0.342463,  -0.235585,  -0.015624,  -0.288630,   0.016393,   0.356302,  -0.470394,   0.044007,   0.184055,  -0.317133,
         0.193141,  -0.302627,  -0.182792,   0.053438,  -0.131155,  -0.036779,  -0.518096,   1.696030,
         0.058062,  -0.085029,   0.163279,   0.057985,   0.026320,   0.197644,  -0.035693,  -0.069734,  -0.454952,   0.832666,   0.682292,   0.006884,  -0.118751,   0.030256,  -0.129968,  -0.142457,   0.040002,  -0.024950,   0.354567,  -0.129048,   0.134348,   0.272110,  -0.015248,   0.001691,  -0.162431,   0.170882,   0.156749,   0.072889,   0.196621,  -0.106022,   0.076104,  -0.133288,
         0.285478,   0.047116,  -0.101628,   0.394607,  -0.052761,   0.009971,  -0.500494,   0.921747,   0.696751,  -0.110210,   0.252620,   0.132169,  -0.168401,   0.231323,  -0.022908,  -0.253703,   0.198448,  -0.018786,  -0.057299,  -0.255229,   0.017138,  -0.258949,   0.143710,  -0.493428,   0.190555,   0.384142,  -0.141358,   0.315926,   0.019910,   0.355004,  -0.220540,  -0.203565,
        -0.082109,  -0.127328,   0.168205,   0.247550,  -0.324814,   1.629680,   0.369325,  -0.171703,   0.208788,  -0.091797,  -0.047964,   0.157869,  -0.098351,  -0.218288,  -0.223139,  -0.057218,  -0.075200,  -0.017952,  -0.064796,  -0.212818,  -0.128024,  -0.122621,  -0.019181,  -0.100902,  -0.281771,  -0.479688,   0.157852,   0.059302,  -0.090972,  -0.223678,   0.012322,  -0.000951,
         0.008750,   0.150242,  -0.135585,   0.634525,   0.167815,  -0.197640,   0.023406,   0.018103,   0.177180,  -0.032859,  -0.264929,  -0.384663,   0.587790,  -0.350569,  -0.093021,   0.025817,   0.299057,  -0.217460,   0.074760,   0.012488,   0.343110,   0.065257,   0.076929,   0.042422,   0.034323,  -0.164341,   0.541082,   0.257817,   0.089829,   0.162824,   0.028584,  -0.111561,
        -0.099454,   1.076890,  -0.181534,  -0.533165,  -0.018113,   0.122879,   0.138617,  -0.155516,   0.023704,  -0.246769,  -0.031489,  -0.126045,  -0.215113,   0.004320,  -0.010787,   0.089444,  -0.132203,   0.465024,   0.094232,   0.079341,  -0.245544,   0.252976,  -0.080429,  -0.165737,   0.009398,  -0.228241,  -0.066373,   0.144527,   0.009114,   0.051074,   0.115234,  -0.155034,
        -0.026945,  -0.375425,  -0.059687,   0.339844,  -0.033511,  -0.159642,  -0.057921,   0.025293,  -0.254685,   0.129031,  -0.266858,  -0.057136,  -0.221073,   0.198389,   0.142281,  -0.198423,   0.237300,   0.187573,   0.310975,   0.011666,  -0.194349,  -0.081293,   0.121787,  -0.107339,   0.127688,   0.077552,   0.101855,   0.043140,   0.189393,  -0.962291,  -0.133096,  -0.095863,
        -0.076791,   0.078374,   0.142082,   0.045114,  -0.179604,  -0.135632,   0.027547,  -0.192207,  -0.151502,   0.057866,  -0.076617,  -0.069667,   0.114320,  -0.015511,   0.224529,  -0.196382,   0.037507,  -0.267625,   0.276111,  -0.201667,   0.288607,   0.133404,   0.063304,   0.225232,  -0.170061,  -0.108800,   0.241958,  -1.561200,  -0.175846,  -0.528845,   0.111697,  -0.237915,
        -0.112592,  -0.273847,  -0.136117,  -0.160090,   0.230394,  -0.097230,  -0.185515,   0.053239,  -0.133743,   0.041142,  -0.190197,  -0.235469,   0.026583,   0.170810,   0.145461,  -0.288951,   0.303083,   0.018284,   0.056051,   0.155153,  -0.019614,  -0.221620,   0.048203,   0.032747,   0.500589,  -1.560420,  -0.726719,  -0.166419,   0.044839,  -0.050807,  -0.500381,  -0.181927,
         0.371119,   0.703782,  -0.030402,  -0.087425,  -0.286968,  -0.038405,  -0.187109,   0.281920,  -0.106258,   0.078886,   0.080787,   0.113733,  -0.043941,   0.073165,   0.103721,  -0.010573,   0.514746,   0.104509,  -0.025920,   0.430871,   0.192028,   0.056251,  -0.811753,  -0.923083,   1.622600,   0.211726,   1.006350,   0.700885,   0.062161,   0.093571,   0.380606,   0.228190,
         0.296362,  -0.144539,   0.080174,  -0.129244,   1.208340,  -0.273911,  -0.203650,  -0.186613,   0.150057,  -0.180151,   0.166046,  -0.095994,   0.044720,  -0.014807,  -0.173592,  -0.189316,   0.233504,  -0.420160,   0.009255,   0.057281,  -0.766401,   0.057274,   0.061465,  -0.003182,   0.123143,  -0.275942,  -0.535044,  -0.462532,   0.449245,   0.510442,   0.562477,  -0.234061,
         0.857438,   0.108307,   0.443937,  -0.004487,  -0.358463,  -0.386747,  -0.051910,  -0.175340,   0.092208,   0.159565,  -0.037212,   0.016901,   0.136158,  -0.149658,   0.180072,   0.278698,   0.725796,   0.200240,  -0.799615,   0.157723,  -0.498839,   0.314136,   0.723894,  -0.557567,  -0.615384,  -0.130322,   0.154423,  -0.077014,  -0.782020,   0.096819,   0.490099,  -0.124075,
         0.447443,   0.106176,  -0.071249,  -0.239707,  -0.313431,  -0.274659,   0.349307,   0.064292,   0.733561,   0.245265,  -0.045869,  -0.212960,  -0.118266,  -0.241540,   0.289732,   0.225751,  -1.770120,   0.956270,  -0.992141,   0.149076,   0.392608,  -0.226885,  -0.270074,  -0.150252,   0.486126,   0.364688,  -0.640429,  -0.084566,   0.015244,   0.010220,   0.467446,   0.226604,
        -0.089307,  -0.358783,   0.200081,  -0.296977,   0.557870,   0.219334,   0.311295,   0.115905,   0.452580,   0.381511,   0.151176,   0.315754,   0.158107,  -0.025822,  -1.893120,   1.272050,  -0.657141,   0.148046,   0.446896,   0.342516,  -0.873723,  -0.159191,   0.102070,  -0.332931,  -0.340487,  -0.598312,   0.504582,   0.491103,  -0.285496,   0.386596,   0.236562,  -0.041295,
         0.065254,  -0.178008,  -0.169761,   0.037738,   0.340516,   0.436528,   1.719130,  -0.357672,
         0.071082,  -0.094131,  -0.018413,  -0.042657,   0.309024,  -0.082426,   0.040547,   0.005810,  -0.312465,   1.571460,   0.056757,  -0.028056,  -0.270766,   0.429837,   0.297606,   0.365366,  -0.148905,  -0.385791,   0.216383,   0.376601,  -0.091000,  -0.098403,  -0.312938,  -0.207956,  -0.336810,  -0.127848,  -0.223417,   0.175561,   0.078273,   0.006701,   0.094404,   0.122798,
         0.204193,  -0.031830,   0.234090,  -0.385973,   0.102114,   0.142330,  -0.059760,  -0.042458,   0.239599,   0.045521,  -0.417305,  -0.300123,  -0.526059,  -0.039818,  -0.383115,  -0.505567,   0.269558,   0.180793,  -0.052352,   0.013097,  -0.462719,   0.234963,  -0.005749,  -0.315692,   0.075635,  -0.199441,  -0.092034,  -0.026110,   0.418442,   0.359099,   0.197144,  -0.014441,
         0.024356,   0.042483,  -0.035859,   0.269580,   0.047096,   0.068123,  -0.274522,  -0.502852,   0.141409,  -0.116747,   0.157833,   0.037244,   0.107444,   0.030669,   0.299389,   0.139407,  -0.084174,   0.001435,   0.237610,  -0.102870,   0.065155,  -0.379911,   0.407295,   0.057788,   0.146297,   0.174996,   0.218060,  -0.054693,   0.285588,  -0.045212,   0.117540,   0.098249,
         0.237676,  -0.039836,   0.234134,  -1.266080,  -0.366198,   0.109626,  -0.283891,  -0.062476,   0.084403,  -0.320377,   0.055452,   0.127286,   0.253356,  -0.076309,  -0.276815,  -0.033329,  -0.481327,   0.205840,  -0.140101,  -0.138617,  -0.413142,   0.137210,   0.500544,   0.069509,   0.092176,  -0.088930,  -0.108993,  -0.349899,   0.632099,  -0.144683,  -0.107962,   0.139487,
         0.174117,  -0.716137,  -0.876933,   0.154671,  -0.723365,  -0.419274,  -0.925450,  -0.056670,  -0.750221,  -0.902874,  -0.332829,  -0.405127,   0.013346,  -0.259723,  -1.118940,  -0.109596,   0.078453,  -1.362650,   0.110355,   0.080152,  -0.453010,  -0.467922,  -0.123000,   0.313298,   0.258451,   0.009563,   0.178080,  -0.176361,   0.334885,   0.240506,   0.315603,  -0.803654,
        -0.876154,   0.315791,  -0.227722,   0.069513,  -0.542578,   0.667704,  -0.071516,  -0.351273,  -0.180141,   0.114361,   0.426607,   0.201191,  -0.048225,   0.479209,   0.035762,  -0.518576,   0.219512,   0.143775,  -0.151021,   0.175222,   0.454757,  -0.060235,   0.337453,  -0.111229,  -0.171380,   0.177081,   0.249140,   0.186215,   0.095760,  -0.660457,  -0.174604,   0.105017,
         0.084689,   0.408792,   0.157730,  -0.431933,   0.495744,   0.170772,   0.101061,  -0.175565,   0.245335,  -0.145737,   0.455676,   0.187588,  -0.289597,  -0.377955,   0.061478,   0.193878,   0.348536,   0.020937,   0.251923,  -0.118340,   0.231240,  -0.108238,   0.111731,  -0.137395,   0.211499,  -0.072275,   0.265101,  -1.019250,  -0.147345,   0.573988,  -0.415586,  -0.443125,
        -0.752508,   0.069484,  -0.459176,  -0.379723,  -0.417041,  -0.197489,   0.027670,   0.313046,  -0.451169,   0.191097,   0.186054,  -0.778348,   0.286966,   0.387059,  -0.317748,  -0.071704,   0.147789,   0.239603,   0.020050,   0.009233,   0.245086,   0.076308,  -0.164554,   0.231417,   0.176534,  -0.395508,  -0.090453,  -0.030412,  -0.608227,  -0.205459,  -0.922514,   0.480938,
        -0.215054,  -0.947179,  -0.335795,  -0.069079,   0.306004,   0.101661,  -0.179192,   0.115280,   0.086055,   0.297180,   0.045097,  -0.009152,  -0.113242,   0.034338,   0.123680,   0.061708,  -0.019404,  -0.237382,   0.205603,  -0.206405,   0.145607,   0.073938,   0.682690,  -0.202302,   0.007456,   0.438620,   0.274394,   0.023405,   0.524659,  -0.153897,   0.457529,  -0.036543,
         0.476501,  -0.250512,  -0.172017,  -0.449793,   0.212931,  -0.191883,  -0.167622,  -0.069394,   0.340115,   0.277841,   0.077749,   0.149716,  -0.121558,  -0.086855,   0.021163,   0.243082,   0.023870,   0.014566,   0.000974,   0.064262,   0.143661,   0.188873,  -0.523363,   0.167845,   0.190725,   0.209434,   0.232180,  -0.063443,   0.008761,  -0.199412,   0.133658,  -0.347907,
        -0.043736,   0.007832,   0.190973,  -0.268358,  -0.019654,  -0.147200,   0.060267,   0.133095,  -0.034160,   0.000603,  -0.395774,   0.221268,   0.215130,   0.019603,   0.270765,  -0.185883,   0.103128,   0.252828,   0.998028,  -0.400389,  -0.202605,   0.004488,   0.147873,   0.462774,   0.417083,  -0.544002,   0.701370,   0.413325,   0.426682,   0.268865,  -0.396966,  -0.129753,
         0.333884,  -0.057481,   0.168605,   0.229299,   0.327540,   0.436412,   0.147215,   0.270873,  -0.098769,  -0.041432,   0.064748,  -0.166859,  -0.122617,  -0.025865,   0.132181,   0.121463,   0.911587,  -0.904903,  -0.036365,   0.122780,   0.334685,   0.136694,   0.288243,   0.014492,   0.045496,   0.109742,   0.212752,  -0.369227,  -0.075057,   0.286609,   0.387050,  -0.113236,
         0.289943,  -0.745584,  -0.200113,  -0.017101,  -0.035449,   0.328370,   0.131360,  -0.076328,  -0.137913,  -0.028999,   0.342814,  -0.031343,   0.058464,   0.089657,   0.360088,  -0.683066,  -0.325517,  -0.009462,   0.983245,   0.241929,   0.349555,   0.004622,   0.106833,   0.722748,   0.038044,   0.764955,  -0.084353,   0.217482,   0.602807,   0.091557,  -0.237273,   0.102444,
        -0.126598,   0.061479,   0.040119,   0.201557,   0.260871,   1.255490,   1.833440,  -0.143680,
         0.296121,  -0.086516,   0.152490,   0.098008,  -0.097969,   0.204056,  -0.061968,   0.050028,   0.349913,   0.797225,   0.120388,  -0.087140,  -0.104691,   0.365014,  -0.015003,  -0.200829,   0.090603,   0.421711,   0.152728,   0.465923,  -0.485559,  -0.437349,   0.241684,  -0.233435,  -0.139906,  -0.460528,   0.371717,   0.173510,   0.154048,   0.226183,  -0.233183,   0.041919,
        -0.093173,   0.163822,   0.138984,   0.319658,   0.134825,   0.206500,   0.306328,   0.298392,   0.182974,   0.470091,  -0.294533,  -0.129175,   0.099847,   0.199930,  -0.318896,  -0.336302,   0.638883,  -0.047480,  -0.133924,  -0.163901,  -0.473618,  -0.165921,  -0.077230,   0.025551,   0.237438,   0.364726,  -0.218477,  -0.308038,  -0.115646,   0.133739,  -0.061908,   0.142136,
         0.514210,  -0.052995,  -0.092292,   0.486088,   0.205906,  -0.346435,   0.408307,  -0.072261,  -0.361819,   0.082433,  -0.106927,   0.032985,  -0.194514,  -0.330355,  -0.405760,  -0.034970,   0.275897,  -0.138805,  -0.329716,   0.144207,   0.152958,  -0.717249,  -0.235265,  -0.086302,  -0.384450,  -0.202029,   0.365228,  -0.054973,   0.026889,   0.058634,   0.121688,   0.082374,
         0.103772,   0.231446,   0.167613,   0.003918,   0.772435,   0.181271,  -0.059647,  -0.121328,   0.149605,  -0.399365,   0.076881,   0.089585,   0.369462,   0.127658,   0.161860,   0.202699,  -0.250117,   0.290730,  -0.215377,  -0.215633,   0.413606,   0.236457,  -0.227867,  -0.099937,   0.165795,  -0.130041,  -0.152912,   0.134802,  -0.299421,   0.206818,   0.028657,  -0.043508,
         0.083407,  -0.041716,   0.588150,   0.169960,   0.081037,   0.339502,   0.074678,  -0.112661,   0.288033,   0.470099,   0.177987,   0.191180,  -0.354532,  -0.286623,   0.050907,  -0.212626,  -0.193745,  -0.696940,   0.395541,   0.196852,  -0.228596,   0.049225,  -0.269989,   0.059763,  -0.160791,   0.124300,   0.348406,   0.133237,  -0.027543,   0.242484,   0.141401,  -0.237310,
         0.700290,   0.398789,  -0.200157,  -0.200314,   0.159186,  -0.297042,  -0.111040,  -0.208042,   0.227512,  -0.165141,  -0.006209,   0.372979,  -0.422930,  -0.066742,   0.014672,   0.146835,  -0.006086,  -0.079018,  -0.198332,   0.308436,   0.157035,   0.231301,  -0.189754,  -0.003357,   0.482515,   0.081859,   0.108217,   0.243980,   0.184228,  -0.112051,   0.098417,   0.359675,
        -0.156622,  -0.152633,  -0.328338,   0.020430,  -0.233067,  -0.113881,   0.039177,  -0.004694,   0.325730,   0.453809,   0.076837,  -0.100529,  -0.093562,  -0.398207,   0.120459,  -0.079519,  -0.647018,   0.136311,  -0.007767,   0.133756,  -0.077834,   0.256537,   0.161164,   0.179243,   0.090573,   0.132303,   0.176042,   0.117429,   0.154958,   0.262565,  -0.290083,   0.084955,
         0.352231,  -0.063189,  -0.031291,   0.260051,   0.304914,  -0.216257,   0.028769,  -0.073868,  -0.275303,  -0.001392,  -0.323743,  -0.431655,   0.195797,  -0.031642,  -0.365604,  -0.167406,  -0.165229,   0.014745,  -0.211630,   0.048731,   0.097135,   0.450315,   0.050341,   0.100214,  -0.037680,   0.456821,   0.379531,   0.429783,   0.285227,   0.318769,   0.161563,   0.089084,
         0.043359,   0.101261,   0.413930,  -0.001150,   0.064120,  -0.176851,   0.258334,   0.010813,  -0.223249,   0.057724,   0.144297,   0.229788,  -0.538433,   0.099702,   0.120100,  -0.104623,  -0.350349,  -0.142156,  -0.121108,   0.148416,   0.113438,  -0.091052,   0.491898,   0.558418,  -0.326292,   0.121904,   0.052184,   0.210449,  -0.049468,   0.084607,  -0.132845,   0.124289,
         0.562946,   0.185542,   0.060867,  -0.591576,  -0.138886,  -0.304577,  -0.223695,  -0.051587,  -0.024122,  -0.070270,  -0.702030,  -0.450074,  -0.154039,   0.074441,   0.027218,   0.133237,   0.447217,   0.303059,   0.053919,   0.133870,   0.313949,   0.216157,   0.565761,  -0.335038,  -0.594915,   0.093428,  -0.306768,  -0.166754,  -0.455477,  -1.262450,   0.080909,  -0.130796,
         0.278425,  -0.075923,  -0.889099,  -0.038587,   0.152098,  -0.359951,  -0.382504,  -0.066918,  -0.040332,  -0.404883,  -0.037181,  -0.206615,  -0.078806,   0.118887,  -0.192619,   0.278344,  -0.034363,  -0.239933,   0.180108,   0.461607,  -0.343717,   0.127031,   0.333801,   0.692888,   0.648268,  -0.192057,   0.076744,   0.202381,   0.717392,   0.156835,  -0.577457,  -0.405925,
         0.401810,  -0.175438,  -0.056406,  -0.010592,   0.511043,   0.189260,  -0.662142,  -0.138765,  -0.534538,   0.129882,  -0.264471,   0.022623,   0.282682,   0.241423,   0.263723,   0.018465,   0.134043,   0.295872,   0.757601,   0.254845,  -0.065590,   0.028487,  -0.052110,  -0.283315,  -0.600359,  -0.573713,  -0.623609,  -0.056518,   0.177626,  -0.058822,  -0.185219,  -0.138106,
         0.195694,   0.474566,  -0.428591,  -0.299572,  -0.037690,  -0.110349,   0.307536,   0.178823,  -0.293618,   0.079199,  -0.100353,   0.276187,  -0.167743,   0.091829,   0.165209,  -0.095741,   0.861508,  -0.275066,  -0.102205,   0.311967,   0.038418,  -0.139449,  -0.129171,  -0.192477,   0.180451,   0.229808,  -0.039821,  -0.310701,  -0.380454,  -0.238441,  -0.300188,  -0.037487,
         0.270476,   0.164527,  -0.166565,  -0.276642,   0.002969,   1.063550,   3.988010,   1.276220
    };

    //! weights for layer 2
    const double ANN_WEIGHTS_CONTACT_HELIX_SHEET_LAYER_1[ 17] =
    {
        -0.107348,   1.002390,   0.937804,  -1.245940,  -0.779245,   0.977986,  -0.782728,   1.040400,  -0.974420,   0.807166,  -1.637150,  -0.765903,   0.945661,   0.716065,   0.683338,   0.876319,  -0.869628
    };
    //! ANN CONTACT_HELIX_SHEET definition
    double ANN_CONTACT_HELIX_SHEET( const linal::Vector< double> &INP)
    {

      // declare variables
      int nora( 0), norb( 0), wei( 0);
      double *hid;

      // test net size
      BCL_Assert( INP.GetSize() == 423, "wrong input size!");

      // allocate memory
      linal::Vector< double> hidden[ 3];
      hidden[0] = linal::Vector< double>( 423);
      hidden[1] = linal::Vector< double>( 16);
      hidden[2] = linal::Vector< double>( 1);

      // normalize data
      hid = hidden[ 0].Begin();
      for( const double *inp = INP.Begin(); inp != INP.End(); inp++)
        ( *( hid++)) = ( *inp) * ANN_NORMALIZE_CONTACT_HELIX_SHEET_A[ nora++] + ANN_NORMALIZE_CONTACT_HELIX_SHEET_B[ norb++];

      // calculate network
      // calculate layer 1
      wei = 0;
      hid = hidden[1].Begin();
      for( size_t i = 0; i < 16; i++)
      {
        *hid = ANN_WEIGHTS_CONTACT_HELIX_SHEET_LAYER_0[ wei++];
        for( const double *inp = hidden[ 0].Begin(); inp != hidden[ 0].End(); inp++)
          ( *hid) += ANN_WEIGHTS_CONTACT_HELIX_SHEET_LAYER_0[ wei++] * ( *inp);
        *hid = double( 1.0) / ( double( 1.0) + exp( -( ( *hid)))); hid++;
      }

      // calculate layer 2
      wei = 0;
      hid = hidden[2].Begin();
      for( size_t i = 0; i < 1; i++)
      {
        *hid = ANN_WEIGHTS_CONTACT_HELIX_SHEET_LAYER_1[ wei++];
        for( const double *inp = hidden[ 1].Begin(); inp != hidden[ 1].End(); inp++)
          ( *hid) += ANN_WEIGHTS_CONTACT_HELIX_SHEET_LAYER_1[ wei++] * ( *inp);
        *hid = double( 1.0) / ( double( 1.0) + exp( -( ( *hid)))); hid++;
      }

      // denormalize data
      // end
      return hidden[ 2]( 0) * ANN_NORMALIZE_CONTACT_HELIX_SHEET_A[ nora] + ANN_NORMALIZE_CONTACT_HELIX_SHEET_B[ norb];

    }

    //! normalization values A*x+b
     const double ANN_NORMALIZE_CONTACT_STRAND_STRAND_A[ 304] =
    {
       0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.000833333, 0.000833333, 0.000833333, 1.2
    };

    //! normalization values a*x+B
     const double ANN_NORMALIZE_CONTACT_STRAND_STRAND_B[ 304] =
    {
       0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, -0.1
    };

    //! weights for layer 1
    const double ANN_WEIGHTS_CONTACT_STRAND_STRAND_LAYER_0[ 4864] =
    {
         0.227475,  -0.019315,  -0.001585,   0.060740,  -0.006008,   0.149528,   0.100242,   0.091311,   0.725419,   0.212082,  -0.055834,   0.036140,   0.105981,   0.194112,  -0.272369,  -0.185740,   0.118993,   0.130545,   0.636119,   0.158725,  -0.080637,   0.027451,  -0.158590,  -0.054645,   0.069224,  -0.304906,   0.099841,  -0.136135,   0.253685,   0.175862,  -0.032109,   0.010065,
         0.213654,   0.134268,   0.137008,   0.030672,   0.207170,   0.032845,   0.352762,   0.092457,   0.458491,  -0.154006,  -0.171828,   0.110282,  -0.058827,  -0.029607,   0.108340,  -0.092330,   0.081781,   0.044079,   0.083452,   0.106448,  -0.265686,   0.124162,   0.291225,  -0.264831,  -0.230818,  -0.060644,   0.069549,   0.104020,   0.067403,  -0.197128,   0.005042,   0.102109,
         0.052437,   0.214470,  -0.072140,  -0.112698,   0.287506,   0.337481,   0.559451,  -0.208096,   0.327101,   0.317652,   0.506418,  -0.187324,   0.202596,   0.128897,  -0.045349,   0.356146,  -0.392572,   0.079445,  -0.004418,  -0.073766,   0.121132,  -0.136852,   0.049622,  -0.017244,   0.128269,  -0.134618,  -0.228344,   0.038560,   0.001932,   0.011047,   0.122307,   0.099480,
         0.154267,   0.067742,   0.426089,   0.022570,   0.957642,  -0.013910,   0.083311,   0.117029,   0.165815,   0.028456,   0.056643,   0.167322,   0.085888,   0.355890,  -0.201358,  -0.123864,   0.013545,  -0.097109,  -0.071358,   0.588316,  -0.086796,   0.038200,   0.105366,  -0.048334,  -0.142114,  -0.240988,  -0.260669,   0.103989,  -0.015435,   0.155525,   0.206256,   0.054970,
         0.723509,   0.003710,   0.417667,   0.135618,   0.229474,   0.064635,  -0.430999,  -0.033463,  -0.000568,   0.007917,   0.443512,   0.308413,  -0.162044,  -0.278114,  -0.051380,  -0.241566,   0.060020,  -0.095754,  -0.003410,   0.024374,  -0.011620,  -0.010689,  -0.256636,   0.072247,  -0.108635,   0.041780,   0.258636,   0.089599,  -0.079790,   0.182836,   0.596335,   0.485030,
        -0.611893,   0.036989,  -0.411912,  -0.226399,   0.097977,   0.242081,  -0.442539,  -0.217257,  -0.499288,  -0.174191,  -0.052183,  -0.154300,  -0.386590,  -0.097904,   0.109283,  -0.169716,   0.290145,  -0.030840,  -0.226609,  -0.027376,   0.185637,  -0.168658,  -0.111852,   0.135098,   0.197899,   0.077236,  -0.121795,   0.026083,   0.703810,   0.034546,   0.268273,   0.005713,
        -0.212523,   0.174366,   0.171776,  -0.027186,  -0.226184,  -0.426616,   0.670459,   0.011502,  -0.244008,  -0.316894,  -0.188435,   0.155158,  -0.012671,  -0.176234,   0.110456,  -0.253363,  -0.001322,  -0.189197,  -0.346502,  -0.093086,   0.061307,   0.045044,   0.121184,   0.169157,  -0.007914,  -0.130862,   0.407425,   0.333230,  -0.282907,   0.010963,  -0.252475,   0.276090,
         0.339110,  -0.038243,  -0.205915,  -0.132701,   0.515020,   0.101484,  -0.209912,  -0.250548,  -0.114247,   0.115989,  -0.205483,  -0.251972,  -0.034624,   0.026324,   0.183633,  -0.187686,  -0.162461,  -0.067840,   0.000523,   0.025421,   0.145657,   0.159118,  -0.067262,   0.072962,   0.621142,   0.177034,   0.313957,  -0.153974,  -0.042623,   0.075132,   0.099986,  -0.022523,
        -0.116367,   0.082790,  -0.077923,  -0.078833,  -0.268858,  -0.209063,   0.102594,  -0.025192,  -0.067619,   0.017794,   0.015551,  -0.023656,   0.147629,  -0.099228,  -0.042508,   0.059195,   0.113490,   0.183829,   0.112018,   0.012178,   0.018045,   0.203269,   0.230549,  -0.073488,   0.493718,   0.037099,   0.132392,   0.018723,   0.191506,  -0.163632,  -0.053873,  -0.040424,
        -0.073892,  -0.247995,   0.045471,   0.037972,   0.044505,   0.151801,   0.085336,   0.002033,   0.028319,  -0.127917,   0.178480,   0.056386,  -0.036560,   0.366480,  56.644300,  -0.930477,
         0.241221,   0.095995,  -0.049208,   0.089781,   0.195438,   0.039622,   0.185352,   0.190825,   0.339519,   0.236346,  -0.222874,   0.054032,  -0.129523,   0.219693,   0.003094,   0.072610,   0.006536,   0.018280,   0.166360,   0.037586,   0.149771,  -0.143222,  -0.023630,   0.004500,   0.109388,  -0.073178,   0.168586,  -0.055294,   0.184223,  -0.129315,   0.190105,   0.029714,
         0.173805,   0.140113,   0.025562,   0.055196,   0.153466,   0.079024,   0.161106,   0.203395,   0.417737,  -0.041027,  -0.121121,   0.150232,   0.159168,  -0.066413,   0.048693,  -0.067524,   0.378508,  -0.031980,   0.072969,   0.024428,  -0.169371,  -0.089349,   0.066556,  -0.253316,  -0.105403,   0.000936,  -0.086175,   0.007040,  -0.006720,  -0.157194,  -0.205893,  -0.027458,
         0.220581,  -0.015807,  -0.056494,  -0.089302,   0.312285,   0.371503,   0.609927,   0.129867,  -0.589062,   0.130389,   0.260098,   0.143136,  -0.257413,  -0.494054,   0.282510,  -0.084005,  -0.081323,   0.100609,  -0.397001,   0.027939,   0.284338,  -0.287312,  -0.085838,  -0.153585,  -0.040318,  -0.032248,  -0.047813,  -0.132039,  -0.240971,   0.068801,   0.177812,   0.063222,
         0.106684,   0.029719,   0.382016,   0.213420,   0.974396,   0.176832,  -0.385271,  -0.062767,  -0.061007,   0.345667,  -0.124337,   0.008098,   0.451373,   0.014300,   0.055628,   0.160095,  -0.168167,  -0.065126,   0.053741,   0.384127,  -0.104584,  -0.052708,   0.060804,  -0.077438,  -0.047149,  -0.254842,  -0.337674,  -0.125485,   0.225468,   0.258467,   0.128347,  -0.033056,
         0.576130,   0.301860,   0.065761,   0.218827,  -0.077345,  -0.192205,  -0.437744,   0.192926,  -0.136286,  -0.052582,   0.577731,  -0.073814,  -0.155661,  -0.045048,  -0.020445,   0.025691,   0.145270,  -0.076863,  -0.047453,  -0.042481,   0.219079,  -0.125499,  -0.052640,   0.424692,   0.123589,   0.252580,   0.394967,   0.018342,   0.124959,   0.518648,   0.404688,   1.217700,
        -0.935250,   0.138197,  -0.293932,  -0.170744,   0.143781,   0.269732,  -0.179310,  -0.040621,  -0.813699,   0.083219,   0.134516,  -0.091258,  -0.447519,  -0.062152,   0.013954,   0.041429,   0.201126,   0.085595,  -0.085698,   0.312010,   0.197597,  -0.031838,  -0.318602,  -0.075201,   0.176146,   0.064123,  -0.197381,  -0.126407,   0.344808,   0.180975,   0.326842,   0.090758,
        -0.345672,  -0.001776,   0.196894,   0.245881,  -0.432441,  -0.277580,   0.301371,   0.037795,   0.014333,  -0.255831,  -0.238643,   0.019145,  -0.064581,   0.110767,  -0.010034,  -0.132976,  -0.173784,  -0.078720,  -0.207375,  -0.121081,  -0.079883,   0.065653,   0.200843,   0.037652,  -0.015634,   0.077193,   0.282591,   0.525643,   0.135587,   0.172987,  -0.504190,   0.014505,
         0.140039,   0.039089,  -0.561024,  -0.416338,   0.302853,  -0.045109,  -0.146958,  -0.057570,  -0.326216,   0.195797,  -0.020240,  -0.318707,  -0.160267,  -0.254818,  -0.188966,  -0.226403,   0.093822,   0.036982,  -0.058206,   0.019183,   0.105674,  -0.006757,  -0.011621,   0.109373,   0.476119,   0.233927,   0.373048,   0.065959,  -0.034496,   0.042566,   0.105287,   0.008730,
        -0.201448,  -0.090822,  -0.038001,   0.074618,  -0.066853,  -0.006813,  -0.095873,  -0.094332,  -0.044677,   0.262147,   0.045419,  -0.041855,  -0.098340,  -0.240922,   0.081373,   0.174575,   0.024914,   0.118872,   0.176737,   0.181279,   0.193875,   0.134170,   0.256736,   0.107598,   0.302057,   0.183265,  -0.137384,  -0.145842,  -0.118829,   0.069638,  -0.167802,  -0.141607,
         0.062150,  -0.006249,   0.054205,   0.012396,  -0.146203,   0.238996,   0.124314,  -0.078309,   0.001638,  -0.134490,  -0.122945,  -0.130892,   0.085436,   0.018896,  67.538700,  -1.107760,
         0.021036,   0.000255,   0.405102,   0.086917,   0.037153,  -0.247177,   0.391074,  -0.158392,  -1.170520,   3.047300,  -0.955792,  -0.102406,   0.019574,  -0.501777,  -0.750429,  -0.081869,  -0.187870,  -0.314800,   0.366934,  -0.344631,   0.132503,   0.151019,  -0.486933,  -0.096842,   0.311534,  -0.904902,  -0.300771,  -0.215989,   0.584313,   0.214436,   0.746724,   0.175051,
         0.069928,   0.037957,  -0.166402,  -0.060705,   0.613666,   0.139692,  -0.595771,   3.107720,  -1.696080,   0.826855,   0.167343,  -0.959819,  -0.450103,   0.555061,  -0.052768,   0.045453,  -1.010220,  -0.124979,   0.548320,  -0.074516,  -0.329209,   0.113646,   0.191257,  -0.326793,   0.395903,  -0.053653,   0.296984,   0.287941,   0.612809,   0.245111,   0.298215,   0.024490,
        -0.061597,  -0.364782,   0.555450,   0.240633,  -1.083980,   0.187415,  -0.532234,  -0.115176,   0.305034,  -0.335256,  -0.050625,   0.188512,   0.193936,   0.579787,  -0.416261,   0.014672,   0.257704,  -0.207542,  -0.350667,  -0.173812,   0.070842,   0.093718,   0.464509,   0.461563,   0.171957,   0.367784,   0.423971,  -0.065623,   0.080394,  -0.160006,  -0.180841,  -0.241766,
         0.198596,  -0.075145,  -1.579810,   0.122453,  -0.082008,  -0.021031,  -0.401512,  -0.257345,   0.033023,   0.243028,  -0.258343,  -0.312351,   0.274586,  -0.325319,   0.041583,   0.158034,  -0.516275,   0.144323,   0.364993,  -0.966395,   0.303664,  -0.216760,  -0.320743,   0.090992,  -0.099636,  -0.031687,   0.040739,  -0.151090,   0.461003,  -0.325484,   0.086443,   0.109023,
        -1.740820,  -0.953988,   1.374400,  -0.092292,  -0.077677,   0.064073,   0.247381,   0.195554,  -0.118167,  -0.097617,   0.555621,  -0.168040,   0.273112,   0.349428,  -0.345398,   0.528460,   0.539324,  -0.239364,   0.365149,  -0.010053,  -0.280071,  -0.000722,  -0.059090,   0.189375,  -0.046257,  -0.280069,  -0.174662,  -0.012740,  -0.276251,   0.044581,   0.491666,   0.604453,
        -1.148510,  -0.098127,  -0.276990,   0.069326,   0.113425,   0.458877,  -0.160474,  -0.176626,   0.481666,  -0.454635,  -0.272855,  -0.238837,  -0.310749,   0.118168,   0.034706,  -0.423389,  -0.173585,  -0.111312,   0.022905,  -0.002309,  -0.268764,   0.509582,   0.236682,   0.124594,   0.055512,   0.025788,   0.092702,   0.299942,   0.178737,   0.341575,  -0.893781,  -0.065370,
         0.143975,   0.158475,  -0.047825,   0.251192,   0.057187,  -0.152698,   0.100808,  -0.037916,  -0.391034,  -0.312486,  -0.076631,   0.156236,   0.047345,  -0.295862,  -0.011572,   0.134430,   0.047531,   0.082052,   0.030125,   0.509168,   0.267894,   0.125582,  -0.134519,   0.243761,   0.357212,   0.238402,  -0.231219,  -0.304253,   0.353505,   0.043979,   0.329986,   0.355881,
        -0.210432,   0.128829,  -0.248674,  -0.086031,  -0.301067,  -0.259578,  -0.092088,  -0.100751,   0.054007,  -0.139113,  -0.018222,   0.941198,  -0.176730,   0.117747,  -0.178932,   0.025581,   0.081506,   0.208680,   0.227141,   0.158048,   0.284632,   0.121868,   0.122293,   0.142202,  -0.076246,  -0.746512,   1.979310,   0.002053,   0.290559,  -0.016541,   0.005271,  -0.046675,
         0.344670,  -0.026768,   0.031844,  -0.364999,   0.119342,   0.024827,   0.390858,   0.070159,   0.018844,  -0.132306,  -0.101624,  -0.064817,  -0.046482,   0.259066,   0.217840,   0.474538,   0.109441,   0.029628,   0.268251,  -0.001113,   0.108420,   0.220571,  -0.192935,  -0.903140,   1.327470,   0.010460,  -0.043142,   0.372901,   0.400907,   0.185148,   0.126557,  -0.035714,
        -0.127811,  -0.022128,   0.092933,  -0.074574,   0.505948,   0.166478,  -0.088595,   0.773472,  -0.130340,   0.076081,   0.002696,   0.461269,   0.005706,   2.083430,  -9.536830,  -1.272500,
        -0.103909,   0.192056,   0.273183,  -0.034152,   0.208109,  -0.047791,   0.117016,   0.079500,  -0.695969,  -0.956508,   0.536659,   0.236688,  -0.189650,  -0.010561,  -0.421496,   0.183467,  -0.164430,  -0.328277,   0.251315,   0.146395,   0.243687,   0.390078,  -0.231519,  -0.023062,  -0.036001,  -0.249056,   0.108831,  -0.185390,   0.009698,  -0.128642,   0.020609,   0.127966,
        -0.071447,  -0.112339,   0.284076,   0.089223,   0.494524,   0.139768,  -0.162339,   0.119736,  -0.457155,   0.645504,  -0.466308,  -0.798861,  -1.008770,   0.539965,  -0.368901,  -0.074189,  -0.305961,   0.073975,   0.192858,   0.225057,  -0.370731,   0.148434,  -0.147824,  -0.028097,   0.289699,  -0.240896,  -0.119656,  -0.219449,   0.154990,   0.231797,   0.277227,  -0.224834,
        -0.014575,  -0.243697,   0.751943,   0.255008,  -0.227420,  -0.284033,   0.070189,   0.508883,  -0.575161,  -0.565445,  -0.804239,   0.670350,  -0.072230,  -0.137677,  -0.283408,  -0.142917,   0.420714,   0.093430,  -0.655099,  -0.080313,   0.036062,  -1.543740,   0.299287,  -0.042161,  -0.215028,  -0.105854,   0.266783,   0.112994,   0.052456,  -0.088691,   0.010527,  -0.113228,
         0.256408,   0.093628,  -0.538689,   0.250211,  -0.783733,  -0.016396,  -0.360308,  -0.746932,  -0.466173,   0.209364,  -0.292946,  -0.282839,  -0.293474,  -0.321098,   0.107485,   0.087195,  -0.479920,  -0.147437,   0.020137,  -0.633853,   0.219486,  -0.250323,   0.286910,   0.160889,  -0.074045,  -0.003160,  -0.107288,  -0.223745,   0.210917,  -0.115048,   0.102410,   0.152411,
        -0.284835,   0.446619,  -0.511477,  -0.003929,  -0.126972,  -0.401708,  -0.577695,   0.137882,  -0.169418,  -0.336424,  -0.364291,  -0.180264,  -0.073795,   0.106761,  -0.736513,   0.334187,   0.097319,  -0.734919,   0.146789,  -0.334670,   0.159371,   0.075816,  -0.284127,   0.346934,   0.165346,  -0.014341,  -0.021148,   0.158234,   0.098947,   0.047120,  -0.229783,  -1.361540,
         0.381220,   0.224902,   0.278624,   0.532918,   0.191881,   0.146280,   0.118913,   0.094684,   0.217125,  -0.258288,  -0.352748,  -0.183304,  -0.181085,   0.257831,   0.105314,   0.566011,   0.243695,   0.059182,  -0.481439,   0.203273,  -0.034947,   0.052224,   0.025612,  -0.024749,  -0.105220,   0.342073,   0.026890,  -0.043833,   0.231284,  -0.493611,   0.937483,  -0.195825,
         0.276236,   1.096660,   1.026770,  -0.268962,   0.184710,   0.345947,   0.729338,   0.107673,  -0.129463,   0.142888,   0.102172,   0.115343,  -0.093117,  -0.011048,   0.362256,   0.126510,  -0.009227,   0.186655,  -0.310331,   0.269912,   0.149957,   0.102245,  -0.203630,   0.062504,   0.122863,   0.264904,  -0.006129,  -0.252503,   0.876646,  -0.144747,   0.459153,   0.967090,
         1.070010,  -0.040632,   0.284457,   0.581452,   0.666422,   0.131371,  -0.163335,  -0.136708,   0.475410,   0.013464,  -0.115524,   0.610058,   0.207276,   0.260525,  -0.040383,   0.205819,  -0.398436,   0.033304,   0.139994,   0.079267,   0.036152,   0.061560,   0.145255,   0.167257,   0.118842,  -0.273491,  -0.244852,   0.045062,   0.496577,   0.274379,   0.516686,   0.028467,
         0.232011,   0.356607,  -0.115094,  -0.156595,  -0.112921,   0.136262,   0.489637,   0.071095,  -0.221052,   0.508042,   0.110014,   0.499616,   0.010402,   0.122451,   0.232603,   0.143355,   0.144243,   0.086773,   0.085955,   0.219909,  -0.006490,   0.049933,   0.285666,   0.101279,  -0.823092,   0.083713,   0.336227,   0.096109,   0.284001,  -0.068284,   0.221802,   0.228968,
         0.047088,  -0.068646,   0.029485,  -0.051551,   0.496317,  -0.030747,  -0.171297,   0.356577,  -0.002589,   0.120568,   0.130228,   0.076647,   0.022277,   0.861011,   1.108040,  -0.504523,
         0.259877,  -0.058906,  -0.208777,  -0.150781,   0.096733,   0.037297,  -0.037531,   0.004524,  -0.296369,   0.242084,   0.354184,  -0.205283,  -0.183892,   0.236580,  -0.203174,  -0.040125,  -0.141741,   0.275087,   0.331471,  -0.136140,  -0.236219,  -0.055443,  -0.250962,  -0.302766,  -0.059198,  -0.730026,  -0.058609,  -0.038486,  -0.331431,  -0.457171,  -0.051195,   0.103336,
         0.168088,   0.111314,  -0.016384,  -0.085648,   0.053959,   0.118704,  -0.382981,   0.220738,   0.206073,  -0.498529,   0.913798,   0.467760,   0.059027,  -0.022610,   0.419451,   0.348648,  -0.104120,   0.309102,  -0.291281,  -0.106058,   0.551299,   0.090008,  -0.096878,   0.650251,  -0.130081,   0.119085,  -0.046638,  -0.527135,  -0.206259,  -0.032197,  -0.075867,   0.096729,
        -0.078427,   0.288062,  -0.101985,  -0.163865,   0.461631,   0.037903,   0.476567,  -0.325817,   0.500601,  -0.032432,   0.972366,  -0.190747,   0.059072,   0.367334,  -0.010546,   0.471866,  -0.353589,  -0.262656,   0.264807,  -0.177916,  -0.279321,   0.842553,  -0.343274,   0.274815,  -0.134382,   0.024227,  -0.173926,   0.076019,  -0.199699,   0.075605,  -0.048532,   0.069410,
        -0.013251,   0.016953,   0.548736,  -0.122389,   1.763280,  -0.172316,   0.132690,   0.340349,   0.445007,  -0.073822,  -0.178154,   0.571613,   0.265597,   0.206407,  -0.050238,  -0.367274,   0.258521,  -0.338904,  -0.242911,   0.019130,   0.347189,   0.652299,  -0.048968,  -0.347228,   0.239958,   0.035606,   0.131808,   0.128964,  -0.136332,   0.056699,   0.057881,   0.105556,
         0.705148,  -0.181220,   1.168170,   0.239972,   0.410043,   0.012742,   0.092864,  -0.031863,   0.155559,   0.181328,   0.024060,   0.146540,   0.199306,   0.024076,   0.496051,  -0.029291,  -0.154598,  -0.124312,  -0.074648,   0.364634,  -0.191118,   0.073869,   0.217002,  -0.208326,   0.046868,   0.138578,   0.063444,   0.120712,   0.102876,   0.095090,   0.026433,  -0.000532,
         0.693664,  -0.112290,  -0.189676,   0.136179,  -0.008151,  -0.013574,  -0.141226,  -0.272818,  -0.031888,  -0.275478,  -0.260717,   0.106602,  -0.166651,  -0.442019,   0.590656,  -1.306440,  -0.030748,  -0.232402,  -0.350732,  -0.171253,   0.474765,   0.119060,   0.326751,   0.291139,   0.000959,   0.180503,   0.037291,   0.080765,   0.048156,   0.007500,   1.175040,  -0.114968,
         1.038480,   0.632839,   1.424790,   0.347430,   0.540822,   0.511814,  -1.159430,   0.558452,  -0.039472,  -0.550966,   0.526009,   0.124054,  -0.189578,  -0.027349,   0.083921,  -0.080904,   0.128420,   0.154424,  -0.413049,  -0.381517,  -0.312133,  -0.009743,   0.026131,   0.255579,  -0.187208,  -0.112678,   0.167572,   0.149687,   0.006664,  -0.309045,   0.274320,   0.715520,
         1.070020,   0.152312,   0.040548,   0.272638,  -0.070929,   0.127284,  -0.271513,   0.078399,   0.497390,  -0.014670,  -0.367383,   1.963870,   0.125512,   0.675593,  -0.531976,  -0.725430,  -0.199162,   0.182606,   0.092598,   0.202704,  -0.036380,   0.196178,  -0.050150,   0.155651,   0.772038,  -0.430337,   1.507000,  -0.218045,   0.348874,   0.150130,   0.355574,   0.000047,
         0.054036,   0.163111,  -0.707512,   0.352510,  -0.176205,  -0.282418,   0.662684,  -0.134293,  -0.482237,  -0.288594,  -0.037957,   0.148563,   0.216298,  -0.075164,   0.248802,   0.061884,  -0.188333,  -0.090612,   0.027577,   0.323169,  -0.066878,   0.125084,   0.515071,  -0.459770,   1.170770,   0.055175,   0.236844,   0.281274,   0.380716,   0.253905,  -0.165967,  -0.088304,
         0.195511,  -0.171448,  -0.019361,  -0.017996,   0.290028,   0.036612,  -0.099903,   0.032298,   0.210750,   0.245505,   0.446365,  -0.073124,  -0.090096,  -0.988845,  12.001800,  -0.163782,
        -0.000022,   0.213843,   0.122733,  -0.132993,   0.118263,  -0.058726,   0.141243,   0.108862,  -0.662757,  -1.301170,   1.287650,   0.300265,   0.123961,   0.171582,  -0.031164,   0.266411,   0.358555,   0.343084,   0.107462,  -0.055439,   0.194234,   0.315405,   0.139557,  -0.043661,   0.042923,   0.026757,   0.028459,  -0.163279,  -0.324664,  -0.139347,  -0.090684,   0.047722,
        -0.051756,  -0.127983,   0.092547,  -0.067856,   0.231966,   0.092327,  -0.972169,  -0.023438,   0.331071,   0.075586,   0.187463,   0.028948,  -0.302728,  -0.232759,  -0.035448,   0.273859,   0.369272,  -0.042347,   0.169455,   0.210388,   0.117233,  -0.039850,   0.011898,  -0.558262,   0.152509,  -0.365696,  -0.220740,   0.069182,   0.132470,   0.187847,   0.474651,   0.033561,
        -0.158116,   0.211713,   0.899826,   0.287835,  -1.087490,   0.493487,  -0.673290,   0.505847,   0.263478,  -0.457732,  -0.663906,   0.165832,   0.474893,   0.436837,  -0.252528,   0.187391,   0.620990,   0.043354,   0.362725,  -0.003528,  -0.147518,  -0.959994,   0.116223,   0.283054,   0.307202,   0.083932,   0.578299,   0.342341,   0.504426,   0.078714,  -0.185227,   0.013315,
         0.425005,   0.246847,  -1.554440,   2.526980,  -1.573880,   0.047978,   0.121390,  -0.812081,  -0.372474,   0.296394,  -0.023300,   0.049300,  -1.239180,  -0.207452,   0.407074,  -0.240218,  -0.090629,  -0.059772,   0.019488,  -0.387332,   0.129324,   0.333007,   0.304206,   0.143247,   0.560342,   0.218665,   0.179211,   0.010113,   0.166151,   0.213587,   0.170484,   0.219917,
        -1.397380,   4.619390,  -0.938636,   0.040416,   0.006322,  -0.617480,  -0.269622,   0.168008,  -0.177782,  -0.160524,  -0.809362,   0.036022,   0.329455,  -0.181319,  -0.585763,   0.202020,   0.270375,  -0.118188,  -0.173960,  -0.098392,  -0.171563,   0.144524,   0.226094,   0.154034,   0.215516,   0.003122,  -0.362615,   0.153389,  -0.078440,  -0.071969,  -0.506204,  -0.672571,
         0.293521,  -0.021996,   0.262856,   0.103871,  -0.001584,  -0.204827,  -0.008932,  -0.243569,   0.291438,   0.005224,  -0.094288,  -0.375727,   0.052496,   0.014277,   0.159654,   0.290966,   0.128371,   0.190392,  -0.203219,   0.464607,  -0.091035,   0.319132,   0.138552,   0.088264,  -0.207152,   0.308728,   0.214437,   0.088465,  -0.241641,  -0.132405,   0.218240,   0.055178,
         0.182728,   0.436342,  -0.234046,  -0.017547,   0.098993,  -0.284595,   0.769977,   0.152636,  -0.265093,  -0.266843,  -0.098654,  -0.077463,   0.087149,  -0.041555,   0.298076,  -0.026709,   0.088104,   0.052036,  -0.260141,   0.292382,   0.259444,  -0.006465,  -0.254553,  -0.156821,   0.189868,   0.112213,  -0.464214,  -0.137229,   0.188757,  -0.013949,  -0.106713,   0.503905,
         0.185655,  -0.034035,   0.037299,  -0.164813,   0.712827,   0.079864,  -0.260494,  -0.117883,  -0.180549,   0.233236,  -0.097529,   0.054040,   0.203692,   0.064281,   0.186804,   0.015041,  -0.413793,  -0.035399,  -0.069748,  -0.021043,   0.244981,  -0.037771,  -0.025306,   0.043040,  -0.631124,  -0.205775,  -0.341695,  -0.233744,  -0.070314,   0.371078,   0.221763,   0.159680,
         0.025294,  -0.151663,   0.483963,  -0.171442,  -0.182079,  -0.065041,  -0.071974,   0.173836,  -0.269455,   0.022773,  -0.072715,   0.003782,  -0.039105,  -0.079105,  -0.406551,  -0.082546,  -0.037316,  -0.166650,   0.002569,   0.326662,   0.002802,  -0.011430,   0.033102,   0.018906,  -0.554220,   0.024449,  -0.117557,   0.112269,   0.096725,  -0.037689,   0.165446,  -0.172783,
         0.388237,  -0.241691,  -0.367495,  -0.058090,  -0.200095,   0.117788,  -0.065018,  -0.140165,   0.052526,  -0.015848,  -0.118356,   0.009723,  -0.508699,   0.046211,   5.301060,   0.498167,
         0.216863,   0.046207,   0.016255,   0.223344,   0.075827,   0.085190,   0.212425,   0.224145,   0.422740,   0.299695,  -0.197700,   0.109963,   0.062826,   0.152263,  -0.112353,  -0.032791,  -0.068686,   0.164510,   0.176922,   0.080996,   0.071931,  -0.108208,   0.077281,   0.034879,   0.136471,  -0.234199,   0.061582,   0.086615,  -0.224112,   0.003969,   0.107169,   0.001337,
         0.069530,   0.113827,   0.024468,  -0.107296,   0.002152,   0.189772,   0.346562,   0.201928,   0.288096,  -0.106340,  -0.108000,   0.107654,   0.255050,  -0.066278,   0.069214,  -0.020137,   0.279553,   0.014162,   0.095554,  -0.115192,  -0.199243,   0.006630,   0.015399,  -0.264651,  -0.081352,  -0.102805,  -0.135318,  -0.034269,   0.088926,  -0.228323,  -0.163908,   0.091792,
         0.028431,   0.141959,  -0.061363,  -0.123138,   0.401425,   0.152724,   0.513961,  -0.007293,  -0.095078,   0.329832,   0.633508,  -0.220839,  -0.049433,  -0.182950,   0.273479,   0.241358,  -0.204023,   0.048696,  -0.044247,  -0.056093,   0.029648,  -0.203889,  -0.031303,  -0.058497,  -0.143435,  -0.112623,  -0.192809,  -0.106800,  -0.282902,   0.147807,   0.137160,   0.170594,
         0.129030,   0.001460,   0.367106,   0.018990,   0.802309,   0.116789,  -0.009551,   0.197784,   0.088399,   0.055442,  -0.030290,   0.123094,   0.327589,   0.116659,  -0.119077,  -0.088763,   0.126763,  -0.151800,  -0.017758,   0.686406,   0.061894,   0.053037,  -0.085847,  -0.085058,   0.001055,  -0.119608,  -0.279513,  -0.022514,   0.104327,   0.146084,   0.064781,  -0.061815,
         0.678772,   0.120803,   0.326444,   0.203331,   0.064504,  -0.035736,  -0.294338,  -0.073337,  -0.193900,   0.064148,   0.545843,   0.152004,  -0.179879,  -0.106081,  -0.035173,  -0.162962,   0.007056,  -0.112231,  -0.029082,  -0.048819,   0.065023,  -0.101168,  -0.047098,   0.109754,  -0.131090,   0.046227,   0.132084,  -0.022663,   0.176297,   0.224460,   0.690535,   1.067330,
        -1.130790,   0.165076,  -0.352253,  -0.044317,   0.262908,   0.284970,  -0.173745,   0.042234,  -0.584214,   0.174105,  -0.221360,  -0.220226,  -0.272532,  -0.165627,   0.020866,   0.144448,   0.311966,   0.267055,  -0.440512,   0.080419,   0.089569,  -0.045280,  -0.194275,  -0.118436,   0.370701,  -0.014109,  -0.071989,   0.123330,   0.483527,   0.211978,   0.131043,  -0.090398,
        -0.494932,   0.004109,   0.221195,   0.144239,  -0.536114,  -0.478648,   0.475549,  -0.219157,   0.026177,  -0.169780,  -0.472789,  -0.057493,   0.006956,   0.008933,   0.168588,  -0.134166,  -0.157659,  -0.075681,  -0.099295,  -0.100081,  -0.077822,   0.026923,   0.055829,   0.196272,  -0.060735,  -0.150287,   0.220349,   0.546568,  -0.166810,   0.019767,  -0.178525,   0.292026,
         0.420027,  -0.211877,  -0.307791,  -0.149059,   0.252572,   0.221551,  -0.087894,  -0.163243,  -0.129924,   0.025654,  -0.062948,  -0.661127,   0.139832,  -0.014119,  -0.085978,  -0.122818,  -0.165162,  -0.121029,  -0.131142,   0.116601,   0.033734,   0.066949,   0.093793,  -0.012823,   0.472084,   0.299176,   0.069979,   0.106607,  -0.178185,  -0.014202,   0.054958,  -0.048514,
        -0.305653,  -0.023555,   0.140689,  -0.016569,  -0.113794,  -0.192530,  -0.005909,   0.006595,  -0.176618,   0.162592,  -0.086712,   0.094611,  -0.109287,  -0.209811,  -0.108949,   0.109264,   0.118857,   0.083163,   0.116317,   0.006171,   0.066841,   0.179812,   0.210839,   0.087637,   0.279073,   0.007739,  -0.131531,   0.116213,   0.021692,  -0.188983,   0.072291,   0.026025,
        -0.043931,   0.028140,  -0.071961,  -0.143741,  -0.116135,   0.028927,  -0.017153,  -0.218472,   0.027276,   0.051121,  -0.119653,  -0.023498,   0.028538,   0.198209,  87.645700,  -0.890238,
         0.149300,  -0.039513,  -0.079863,   0.041129,   0.099204,   0.043379,   0.065883,   0.152610,   0.418171,   0.277419,  -0.068923,   0.013908,   0.098549,   0.294746,   0.137144,  -0.147462,   0.160476,   0.199128,   0.438304,   0.070450,  -0.062315,  -0.080009,   0.078940,  -0.139187,   0.067655,  -0.186647,   0.095495,  -0.003124,  -0.036603,  -0.017156,   0.073866,   0.114369,
        -0.007120,   0.261327,   0.159051,  -0.109410,   0.022078,   0.219391,   0.207628,   0.150702,   0.512549,  -0.163103,  -0.118394,   0.173732,   0.094533,  -0.102149,  -0.025451,  -0.018393,   0.148354,   0.049238,  -0.010383,  -0.064807,  -0.115341,  -0.011340,   0.042883,  -0.073572,  -0.161179,  -0.040875,  -0.155794,  -0.088549,   0.199976,  -0.116236,  -0.051954,   0.142681,
         0.009096,   0.111289,  -0.181085,   0.000724,   0.398969,   0.240659,   0.588853,  -0.097884,  -0.025435,   0.350361,   0.469276,  -0.056686,  -0.085502,  -0.155988,   0.292518,   0.221641,  -0.239204,   0.125134,  -0.049377,  -0.235127,  -0.035415,  -0.107518,  -0.068406,   0.022494,  -0.098231,  -0.146201,  -0.097898,  -0.102693,  -0.208461,   0.087099,   0.200905,   0.087157,
        -0.034018,   0.113554,   0.314421,   0.025190,   0.866288,   0.131164,   0.069168,   0.160028,   0.273114,   0.034059,  -0.029046,   0.118913,   0.326384,   0.178741,  -0.070325,  -0.148165,  -0.047765,  -0.165098,  -0.170778,   0.759868,  -0.003523,  -0.024795,  -0.096057,  -0.135171,  -0.023536,  -0.170472,  -0.389971,  -0.120957,   0.223315,   0.026578,  -0.010392,   0.119238,
         0.597741,   0.140349,   0.324583,   0.309192,   0.023886,  -0.029458,  -0.245318,  -0.129010,  -0.213936,  -0.060221,   0.455093,   0.070742,  -0.148565,  -0.009812,  -0.017694,  -0.053831,  -0.032586,  -0.036894,   0.072646,   0.082919,  -0.183637,  -0.188060,   0.017072,   0.118291,  -0.095526,  -0.052376,   0.311474,   0.122966,   0.056373,   0.222425,   0.602258,   0.826877,
        -0.995269,   0.125356,  -0.355491,  -0.055388,   0.070382,   0.156967,  -0.186756,  -0.117860,  -0.340589,  -0.022517,  -0.319544,  -0.347029,  -0.279289,  -0.225873,   0.166648,  -0.086151,   0.218855,   0.025529,  -0.205328,   0.027017,   0.023250,  -0.101611,  -0.269289,  -0.065630,   0.310106,  -0.074873,  -0.047447,   0.049016,   0.461202,   0.138846,   0.293196,   0.023914,
        -0.405519,   0.092971,   0.254105,   0.264306,  -0.599680,  -0.639278,   0.550354,  -0.325891,   0.032776,  -0.299623,  -0.448973,   0.068276,   0.005717,  -0.043291,   0.025687,  -0.035432,  -0.179492,  -0.068939,  -0.108701,  -0.166927,  -0.192514,   0.044615,   0.064105,   0.256543,  -0.149436,  -0.000408,   0.146371,   0.405073,  -0.143481,   0.004121,  -0.107061,   0.188661,
         0.390041,  -0.032007,  -0.386180,  -0.114105,   0.303339,   0.120666,  -0.197985,  -0.366457,  -0.017622,  -0.009119,  -0.058644,  -0.483761,  -0.019148,   0.085396,  -0.186433,  -0.274489,  -0.048951,  -0.045371,  -0.198844,   0.101883,   0.090525,   0.029425,   0.071520,   0.184445,   0.471857,   0.103862,   0.281783,   0.020949,  -0.133870,   0.156366,   0.158816,   0.000974,
        -0.246248,  -0.158049,  -0.103205,  -0.120156,  -0.220713,  -0.104673,   0.025766,  -0.006150,  -0.040030,   0.140537,  -0.038587,  -0.035156,  -0.062673,  -0.126014,   0.105597,   0.050267,   0.059748,   0.182774,   0.059310,   0.086420,   0.096103,   0.068128,   0.144125,  -0.209777,   0.585883,  -0.128101,   0.073066,   0.080924,   0.050167,  -0.180813,   0.023355,   0.049525,
         0.043282,   0.047332,  -0.102390,  -0.078231,   0.125745,   0.030322,  -0.068872,   0.001785,   0.137670,  -0.064979,  -0.088280,  -0.129507,   0.007936,   0.417675,  95.355800,  -1.062620,
         0.130752,  -0.135910,  -0.215418,  -0.024001,  -0.043212,   0.488812,   0.198601,  -0.138848,   0.115510,   0.209941,  -0.283801,   0.259434,  -1.361500,  -0.015762,  -0.233862,  -0.846308,  -0.045148,  -0.074008,   0.708874,  -0.141744,  -0.317350,  -0.150483,  -0.042866,   0.272887,   0.096683,   0.369902,   0.228203,  -0.390787,   0.569045,   0.005779,  -0.236941,  -0.153957,
         0.062486,   0.084596,   0.076477,  -0.041063,  -0.025601,   0.129169,  -1.367030,   0.068423,   0.507872,  -0.595546,   2.309410,   0.319662,  -0.012255,  -0.693005,   0.482036,  -0.425015,  -2.301770,  -0.234666,  -0.115136,  -0.015006,   1.096060,  -0.385227,   0.222254,   5.916100,   0.181849,   0.298142,  -0.583997,  -0.336016,  -0.172001,  -0.222159,   0.115485,   0.133462,
        -0.027181,   0.264655,  -0.586125,  -0.096879,   0.224998,   0.151255,  -2.633050,  -0.073687,   3.424800,   4.579410,  -3.727430,   0.473035,   0.511326,   3.725710,   0.708190,   4.790170,  -0.406249,  -0.734955,   0.007069,  -0.266084,  -0.376254,  -5.318830,   0.254326,   0.916614,   1.151850,   0.541890,  -0.232861,  -0.207428,  -0.550587,  -0.449577,  -0.129821,  -0.258222,
         0.229744,  -0.038751,   2.197540,  -0.098486,   0.065627,   0.691536,   1.659150,   0.092486,  -1.427640,  -0.203778,  -0.073890,   1.074370,   1.048850,  -0.341624,   1.115920,   0.230422,   0.939376,  -0.500863,  -0.148322,  -1.932150,  -0.188905,   0.352580,   2.831760,   0.143104,   0.663600,   0.414340,   0.475486,   0.284772,  -0.086294,   0.142489,   0.204278,   0.285320,
         0.040416,   1.058350,  -0.523346,   0.615603,   0.746391,   0.150036,   0.507031,  -0.921523,  -0.544913,  -0.189499,   0.167441,  -0.908985,   0.597678,  -0.255584,  -0.472879,  -0.532389,   0.913590,  -1.107310,   0.590428,   0.768348,   0.349302,   0.380111,   0.740414,  -0.031672,   0.113937,   0.181972,  -0.145118,   0.466280,   0.197153,   0.241438,   1.745220,  -0.086222,
         0.092068,   0.072132,  -0.575944,   0.122220,  -0.362899,  -0.072996,  -0.388744,   0.289681,   0.190882,  -0.252285,  -2.858380,   0.448798,  -0.706043,   0.050215,   0.181232,  -0.713141,  -0.186436,  -0.373802,  -1.508400,  -0.125880,  -0.398182,  -0.259923,   0.206034,  -0.002244,   0.117728,  -0.240369,   0.086379,  -0.086260,   0.654285,  -0.304222,   1.079940,   0.631816,
         0.137998,   0.232325,   2.666540,   2.107640,  -0.106743,  -0.382257,  -1.990620,   0.452205,  -0.215179,  -0.583227,  -0.239669,  -0.094711,   0.036375,  -2.332330,   0.792230,   0.333045,  -0.652153,  -0.246687,  -0.360662,   0.023971,   0.179679,   0.198881,  -0.035271,   0.092023,  -0.390320,   0.168931,  -0.772083,   0.234954,  -2.894140,   0.107706,   1.152110,   2.387790,
         4.435800,   0.462369,  -0.248997,  -5.253700,   7.156290,   3.008530,  -0.125100,  -0.505349,   2.864640,  -0.173772,  -0.101680,  -2.331310,   0.764909,   0.652967,  -0.329320,   0.054272,   0.018103,  -0.670483,   0.048499,  -0.036795,  -0.203925,   0.219938,   0.080196,  -0.273985,   1.234020,  -0.076064,  -0.299477,   0.966671,   0.858943,   0.814845,  -1.274510,   1.347450,
         0.241027,   0.647777,  -0.093452,  -1.089160,  -0.060124,  -0.070395,   1.281000,   0.560139,  -1.527770,  -3.185100,   0.166470,   1.239160,   0.347494,   0.041948,   0.428247,  -0.314834,   0.211797,   0.172061,  -0.019190,   0.111905,  -0.030222,   0.128129,  -0.207259,   0.324990,  -0.821164,  -0.673720,   0.788780,   0.752634,   0.231609,   0.213636,  -0.610496,   0.101581,
        -1.052380,   0.998131,   0.078643,  -0.063628,   1.163040,   0.011086,   0.772242,   0.232969,   0.244437,   0.672948,   0.056817,   0.971339,   0.685537,   0.065357,  25.659500,  -0.053917,
         0.137734,  -0.056481,   0.077950,   0.116201,  -0.109067,   0.169447,   0.114660,  -0.148396,   1.087210,   0.042723,   0.027280,  -0.289242,   0.127627,   0.511782,   0.045777,  -0.315997,   0.054534,   0.469152,   0.469044,   0.585265,  -0.774123,  -0.345901,   0.185049,  -0.144120,  -0.004695,  -0.049822,  -0.046008,  -0.083099,   0.020246,  -0.128661,  -0.377037,  -0.116772,
         0.149658,   0.101829,  -0.036681,  -0.148401,   0.255889,  -0.018020,  -0.050415,   0.033674,   0.465923,  -0.111698,   0.476382,   0.080952,  -0.299583,  -0.075740,   0.293371,   0.120926,  -0.169222,   0.028657,  -0.019386,   0.081104,   0.369001,   0.106355,   0.245943,   0.150098,   0.078411,   0.116377,   0.294770,   0.042200,   0.007940,  -0.097232,   0.190536,   0.106560,
        -0.136527,   0.216072,  -0.224116,  -0.114762,   0.597620,   0.144492,  -1.021650,  -0.375272,   1.838070,   0.456393,  -0.587891,  -0.368382,   0.441726,   0.932828,  -0.781257,   0.596582,  -0.698340,  -0.304072,   0.546565,  -0.345408,  -0.233828,  -0.429522,   0.007678,   0.375920,   0.305102,  -0.130848,  -0.209772,   0.083515,   0.065905,   0.007529,   0.183889,  -0.061114,
         0.088905,   0.154737,   0.551370,   0.175980,   0.242212,   0.065817,   0.424068,  -0.033561,  -0.162744,   0.058001,  -0.001153,   0.085786,   0.281601,  -0.127050,   0.407847,   0.241389,   0.341295,  -0.159032,  -0.010092,  -0.182485,  -0.059542,   0.347603,   0.219844,   0.371674,   0.334839,   0.126670,   0.272554,   0.206681,   0.076416,   0.156319,   0.280649,   0.113385,
        -0.416434,   0.369002,   0.077051,   0.130940,   0.241068,   0.020196,   0.115677,  -0.342202,   0.202563,  -0.130812,  -0.111170,   0.253264,   0.125830,  -0.023276,   0.072244,  -0.277099,   0.195657,  -0.288777,  -0.016030,   0.133730,  -0.055073,  -0.274525,   0.189993,  -0.349062,   0.124409,   0.171439,  -0.217522,   0.136651,   0.113152,  -0.120892,   0.924087,  -0.008019,
        -0.032863,  -0.007323,   0.116698,   0.513198,   0.158963,  -0.455263,   0.008981,   0.105455,   0.069467,  -0.100093,  -1.044870,  -0.282413,   0.171468,  -0.130679,   0.103712,   0.045720,   0.118374,  -0.020749,  -0.039854,  -0.385226,  -0.365427,   0.162292,   0.138020,   0.145366,  -0.045347,  -0.086599,   0.201293,   0.076395,   0.533180,  -0.063009,   0.553382,   0.004497,
         0.397009,  -0.006127,   0.532066,   0.100717,   0.280429,  -0.230910,  -0.250370,  -0.032630,  -0.108404,  -0.187238,   0.050291,   0.087764,  -0.020312,  -1.197080,   0.338121,   0.044766,   0.298049,   0.129282,  -0.105576,  -0.126525,   0.201562,   0.138309,  -0.084671,   0.262618,  -0.398684,  -0.194608,   0.148533,   0.165139,  -0.894782,  -0.519746,   0.762279,   0.739101,
         0.332523,  -0.034711,   0.242955,   0.637459,  -0.670083,   0.392132,  -0.406156,  -0.121394,   0.233052,  -0.283630,   0.036739,  -0.321075,   0.099782,   0.573380,   0.710125,   0.243922,  -0.289069,   0.100298,   0.264028,   0.117955,  -0.027922,   0.153841,   0.002624,   0.076019,   0.251065,  -0.020261,   0.208455,  -0.038831,   0.443157,  -0.124265,  -0.192298,  -0.229760,
         0.109756,  -0.133073,  -0.273014,  -0.283547,   0.130798,   0.065863,   0.411531,  -0.012092,  -0.047203,  -0.479049,  -0.133507,  -0.084021,   0.641965,   0.364990,   0.100278,   0.127373,   0.498166,   0.436007,  -0.036699,  -0.010140,  -0.045872,   0.113459,  -0.728125,   0.033586,   0.466000,  -0.312181,   0.086079,   0.291970,   0.168830,  -0.092417,   0.005218,  -0.066915,
        -0.378998,  -0.346719,   0.165526,   0.088517,   0.231605,  -0.236765,   0.166366,  -0.325130,   0.175133,   0.458070,   0.472580,   0.148215,   0.177919,  -0.043206,  12.048800,  -0.741820,
        -0.040667,  -0.032193,   0.037061,   0.006483,   0.175515,   0.118734,   0.038321,  -0.053113,  -0.524706,   0.016843,  -0.774959,   0.178091,   0.042793,   0.101218,  -0.014448,   0.091172,   0.117025,   0.321499,   0.045346,  -0.092086,  -0.240418,  -0.128915,  -0.055314,  -0.105615,   0.092938,  -0.011163,   0.040575,   0.241129,   0.061807,   0.135902,  -0.065326,  -0.066501,
         0.091196,   0.161312,   0.053447,   0.091262,  -0.050955,   0.041318,  -0.356531,  -0.119180,   0.110341,  -0.091419,   0.211644,   0.330383,   0.067966,   0.141437,   0.102916,   0.035453,   0.220497,   0.094934,  -0.382887,  -0.226733,   0.016859,   0.063757,  -0.112470,   0.395788,  -0.066490,  -0.100650,  -0.009171,   0.129755,  -0.243610,  -0.069991,   0.145221,   0.033847,
         0.039529,   0.091632,   0.048020,   0.038894,  -0.231849,  -0.202353,   1.306400,  -0.140051,   0.327673,   0.471059,   0.503989,   0.017301,   0.261733,  -0.049713,   0.372520,   0.066125,  -0.592508,  -0.331975,  -0.060111,  -0.021758,  -0.211993,   0.627577,   0.120096,  -0.100618,   0.014061,   0.134703,  -0.386476,   0.014287,   0.343513,   0.020591,   0.079942,  -0.006715,
        -0.037377,   0.062984,   0.014535,  -0.418845,   0.460465,   0.192963,   0.236774,   0.268215,  -0.012035,   0.053883,   0.063018,  -0.110016,   0.666669,   0.120313,  -0.342275,   0.014022,  -0.076784,  -0.039907,  -0.170812,   0.302377,   0.006553,  -0.098116,   0.039794,   0.111542,  -0.118155,   0.038971,   0.331935,   0.004118,  -0.004383,   0.170599,  -0.140269,  -0.065854,
         0.398753,   0.002424,  -0.874042,   0.268033,   0.132569,   0.103662,   0.373157,   0.099022,   0.074719,  -0.182009,   0.308944,   0.037240,  -0.353953,  -0.189549,  -0.094393,  -0.074823,  -0.038169,   0.109537,  -0.023160,   0.029035,  -0.099911,   0.058726,  -0.191372,  -0.024268,   0.206014,  -0.140096,   0.359240,  -0.048298,   0.191743,   0.010033,  -1.403850,  -0.947669,
         1.547280,   0.195173,   0.179611,  -0.398215,  -0.342734,   0.395047,   0.383595,   0.611712,   0.243075,  -0.055800,   0.587840,   0.532929,   0.203180,  -0.159879,   0.316657,  -1.060240,  -0.235895,  -0.105715,   0.019706,   0.113203,   0.092858,  -0.009315,   0.181559,  -0.100663,   0.101633,  -0.311353,   0.810544,   0.179885,  -1.166570,   0.403530,  -0.312759,   0.682346,
         0.465565,  -0.341084,  -0.840638,   0.501844,   0.326186,   0.528267,  -0.533240,   0.048999,   0.499826,   0.319550,   0.098028,   0.002371,   0.186896,  -0.886958,   0.164722,   0.056277,   0.172463,  -0.068412,   0.855031,   0.250512,   0.483579,  -0.042945,  -0.260326,   0.062624,   0.943305,   0.320533,  -1.147770,   0.288474,  -0.254735,   0.451126,   0.374213,  -0.764845,
        -1.678720,   0.421489,   0.477849,   0.550666,  -0.422628,  -0.041929,   0.777427,  -0.000264,  -0.032604,   0.123074,   0.328619,  -1.807090,   0.341510,  -0.028685,  -0.037135,   0.309361,   0.939447,   0.306731,   0.190735,   0.014853,   0.035143,   0.069722,   0.464382,   0.064847,  -2.740600,   2.270750,  -2.587820,   0.257602,   0.151602,  -0.764034,  -0.465194,   0.392855,
         0.208347,  -0.109467,  -0.783166,  -0.070495,   0.505049,  -0.268366,  -0.458759,   0.039230,   0.148415,  -0.711046,   0.496851,   0.032947,   0.163651,   0.177902,   0.437130,  -0.093201,   0.152977,   0.101002,   0.215961,   0.125928,   0.395576,   0.017094,  -2.063890,   2.537270,  -0.352640,   0.089323,  -0.107803,  -0.112351,  -0.910537,   0.272788,  -0.022422,  -0.223248,
        -0.462765,  -0.237580,   0.087884,   0.066859,  -0.575504,   0.055566,   0.152431,   0.248263,   0.034060,  -0.249285,   0.483055,   0.185134,  -0.112011,  -0.089088,  -0.497671,   0.087845,
         0.306532,   0.008840,  -0.174418,  -0.111388,   0.252135,  -0.084277,   0.081537,   0.171262,   0.189759,   0.576815,  -0.160017,  -0.029622,  -0.179510,  -0.117598,  -0.541561,  -0.116662,  -0.231791,  -0.049054,   0.094741,  -0.297125,   0.174238,  -0.031097,  -0.408031,  -0.006382,   0.123784,  -0.627456,  -0.150176,   0.021304,  -0.034885,   0.012891,   0.393732,  -0.147007,
         0.020626,   0.080214,   0.124625,  -0.044357,   0.014964,  -0.152313,  -0.190543,   0.181341,   0.438077,  -0.383545,   0.116134,   0.537901,   0.108642,  -0.153115,   0.095262,   0.177199,   0.548013,   0.125680,  -0.193293,  -0.130738,  -0.109513,   0.180830,   0.143592,  -0.290314,  -0.091353,   0.067314,  -0.282593,   0.051607,  -0.076777,  -0.195254,   0.036213,   0.151365,
         0.232006,   0.223995,  -0.103584,  -0.096454,   0.050469,   0.319241,   0.784622,  -0.142608,  -0.035803,  -0.218516,  -0.153807,   0.185238,  -0.096099,  -0.470655,   0.107134,  -0.258530,  -0.042118,   0.220874,  -0.044601,   0.173212,   0.354379,   0.213441,  -0.577875,  -0.035431,  -0.004130,   0.016876,   0.158273,  -0.048584,   0.158968,   0.313579,   0.031829,   0.088461,
         0.007258,   0.105147,   0.353219,   0.159440,   1.277480,  -0.071697,  -0.090471,   0.224203,   0.261703,  -0.074668,   0.052220,  -0.183284,  -0.071582,   0.498467,  -0.274954,  -0.199344,  -0.153227,  -0.123708,   0.145156,  -0.268946,  -0.056077,   0.218986,   0.377113,   0.164644,  -0.104023,  -0.065291,  -0.103076,   0.107291,   0.128658,   0.086399,   0.194695,   0.043384,
         0.348196,  -0.057563,   0.692770,   0.094092,   0.206214,   0.106417,  -0.087570,   0.051106,   0.078803,   0.074291,  -0.201146,   0.322831,   0.006184,  -0.088129,   0.249798,  -0.067290,   0.218233,   0.386183,  -0.035421,  -0.145246,  -0.267789,   0.214217,  -0.103455,   0.003767,  -0.034840,   0.026912,   0.141920,  -0.016582,   0.025171,   0.216286,   0.068969,   0.202910,
         0.122751,   0.047495,  -0.417466,  -0.151750,   0.150965,   0.059601,  -0.425627,  -0.347909,   0.095047,  -0.384841,   0.248136,  -0.003861,  -0.371808,  -0.223074,   0.180587,  -0.533187,  -0.019908,   0.040368,  -0.382702,  -0.152055,   0.794195,  -0.146997,   0.126143,   0.202832,  -0.201491,   0.175204,   0.043358,  -0.255719,   0.022444,   0.207877,   0.910811,  -0.248621,
         0.565738,   0.442689,   0.361091,  -0.141072,   0.239193,   0.271580,   0.461876,   0.612831,  -0.225417,  -0.539046,   0.403551,   0.083152,  -0.066888,   0.012223,   0.279179,  -0.037918,  -0.066726,  -0.116874,  -0.522264,  -0.149529,  -0.009978,   0.068417,   0.340730,   0.071270,  -0.180642,  -0.037722,   0.010752,   0.330272,   0.172228,   0.062753,  -0.560416,  -0.039124,
         0.216835,   0.203865,  -0.631450,  -0.700603,   0.406391,  -0.364405,   0.076290,   0.141802,  -0.494187,   0.289906,   0.064253,   0.562446,  -0.167490,  -0.047261,  -0.035723,  -0.351308,   0.327949,   0.016234,   0.083078,   0.204100,  -0.006145,   0.167202,  -0.125957,   0.201215,   0.353901,  -0.034238,   0.716344,  -0.228205,   0.039240,  -0.009683,   0.011517,  -0.194883,
         0.159165,   0.181987,  -0.103676,   0.363604,  -0.292636,  -0.324196,   0.205887,   0.013757,   0.007493,  -0.248766,  -0.054684,   0.092588,   0.116682,   0.107266,   0.166492,   0.172372,  -0.190130,  -0.086714,   0.282648,   0.162137,   0.074139,   0.119291,   0.173233,  -0.117598,   0.498815,   0.057930,  -0.142247,  -0.194481,  -0.003295,  -0.106734,  -0.194703,  -0.265124,
         0.276084,  -0.196804,  -0.067656,   0.171284,  -0.102516,   0.210686,   0.200639,   0.117676,   0.116539,  -0.103702,   0.024865,   0.161957,   0.066626,  -0.283578,  16.952300,   0.008355,
        -0.224026,   0.038526,   0.227421,   0.014085,  -0.083091,  -0.056064,  -0.058456,   0.012988,  -0.313063,  -0.741192,  -0.246684,  -0.071721,   0.040220,  -0.177565,  -0.100829,   0.015357,   0.121879,   0.166893,  -0.013198,   0.018269,   0.151468,   0.123572,   0.039476,  -0.002888,   0.101276,   0.171604,  -0.025701,  -0.047238,   0.217466,   0.053117,   0.027351,   0.140612,
        -0.002321,  -0.002845,   0.023225,  -0.096624,  -0.056593,  -0.047420,   0.422904,  -0.211110,  -0.790137,   0.120058,   0.086464,  -0.140668,  -0.133743,   0.104727,   0.077520,   0.125609,  -0.552016,   0.001482,   0.262105,   0.218847,  -0.140825,   0.105521,   0.169345,   0.017105,   0.060494,  -0.101833,   0.058535,   0.179704,   0.008031,   0.252921,   0.168790,  -0.153286,
        -0.148350,  -0.298061,   0.449174,   0.218903,   0.191530,  -0.362834,  -1.053230,   0.173830,   0.049662,   0.048610,  -0.250788,   0.067582,   0.250957,   0.450291,  -0.283191,   0.096694,   0.282449,   0.107496,   0.032814,   0.139034,   0.048828,  -0.345300,   0.324975,   0.106472,  -0.002838,   0.206896,   0.058317,   0.018285,  -0.064589,  -0.120416,  -0.088471,  -0.146755,
         0.061976,  -0.140322,   0.006438,  -0.233447,  -1.231030,   0.049478,  -0.072058,  -0.087443,  -0.177410,   0.052591,   0.005632,  -0.110871,   0.020700,  -0.025901,   0.118557,   0.158194,  -0.175962,   0.166784,   0.200940,  -0.217865,   0.243033,  -0.021632,   0.053005,   0.004167,  -0.095337,  -0.031804,   0.041651,  -0.174526,  -0.146314,  -0.132554,  -0.014502,  -0.060950,
        -0.264682,  -0.284153,  -0.815461,   0.151747,   0.182293,  -0.029625,  -0.166882,  -0.048383,   0.222414,   0.198304,   0.018914,  -0.010049,   0.062478,   0.093817,  -0.098408,   0.342500,  -0.102987,  -0.170365,  -0.043680,   0.109231,   0.254437,   0.136729,  -0.013051,   0.376122,   0.227670,   0.067719,  -0.133099,  -0.211238,   0.047853,  -0.104019,  -0.722866,  -0.837087,
        -0.184446,  -0.054031,   0.151068,   0.018257,  -0.099533,   0.000875,   0.085160,   0.120925,  -0.144067,  -0.061440,   0.332224,   0.169491,   0.005728,   0.263763,   0.059948,  -0.109103,  -0.013204,   0.060312,   0.231915,   0.196525,   0.075753,   0.240392,  -0.039844,  -0.177171,  -0.003741,  -0.195899,   0.182822,   0.088112,  -0.267214,   0.114368,  -1.104180,   0.383233,
         0.014196,  -0.215159,  -0.550196,   0.195316,   0.078335,   0.287820,  -0.643509,  -0.047650,   0.222253,   0.330588,  -0.217909,   0.097973,   0.294542,  -0.104087,   0.069927,   0.200399,   0.178106,   0.349815,   0.339729,   0.411566,   0.100565,  -0.001954,  -0.147591,  -0.177420,   0.334567,   0.147902,  -0.298170,  -0.348298,  -0.435519,   0.195086,  -0.089203,   0.180750,
        -0.195950,   0.118096,   0.172740,   0.197580,  -0.561278,  -0.078304,   0.296022,   0.161495,  -0.047122,   0.155335,   0.157328,  -0.826216,   0.226660,   0.210132,   0.050327,   0.310082,   0.305511,   0.126928,   0.167666,  -0.150904,  -0.051439,  -0.084513,   0.074684,  -0.044298,  -0.494982,  -0.192439,  -0.933910,  -0.167163,   0.070582,   0.031131,  -0.160668,   0.020209,
         0.103015,   0.087454,  -0.248997,  -0.012044,   0.241040,   0.193587,  -0.188819,   0.133156,   0.191871,  -0.437716,  -0.024660,  -0.055602,   0.308220,   0.251742,  -0.023989,  -0.067233,   0.108029,  -0.117127,  -0.106054,  -0.098728,   0.105476,  -0.142527,  -0.461155,   0.010555,  -0.778763,  -0.028119,   0.027597,   0.100532,  -0.277331,  -0.036468,   0.021823,  -0.059419,
        -0.276857,  -0.152900,   0.168215,  -0.014692,  -0.108308,   0.187715,   0.002526,  -0.049090,  -0.078717,   0.015796,   0.020688,   0.031904,   0.004620,   0.501877, -10.609400,  -0.404237,
        -0.161374,   0.086822,  -0.015434,  -0.013839,   0.148808,  -0.149270,   0.144833,   0.144476,  -0.795976,  -0.787168,   1.014710,   0.330676,  -0.021924,  -0.345392,  -0.307607,   0.295914,   0.046669,   0.029791,   0.105148,  -0.188578,  -0.033760,   0.419072,  -0.039299,   0.011208,   0.079449,  -0.238622,   0.027705,  -0.264266,  -0.076067,  -0.168341,   0.029203,   0.136192,
         0.007408,  -0.028199,   0.126443,   0.041548,   0.674918,   0.108019,  -0.986449,   0.593142,  -0.182368,   0.395398,   0.498554,  -0.434548,  -0.756837,  -0.226618,   0.052005,   0.266562,   0.037424,  -0.005676,   0.459172,   0.279473,   0.174316,   0.246453,   0.127729,  -0.542899,   0.192851,  -0.291117,   0.701732,   0.038038,   0.359862,   0.446659,   0.994857,   0.218239,
         0.053882,   0.195273,   0.833457,   0.391957,  -1.261590,   1.639830,  -1.011270,   0.440253,   0.446997,  -0.969971,  -1.261180,   0.692247,   0.499906,   0.422951,  -0.536912,  -0.264419,   0.817432,   0.202161,   0.148151,   0.097418,   0.115632,  -1.040180,   0.046694,   0.106930,   0.653821,   0.127509,   0.889320,   0.325207,   0.484875,   0.147795,  -0.353530,  -0.158682,
         0.471232,   0.285430,  -2.061610,   4.443220,  -1.606640,  -0.073755,  -0.052934,  -0.711735,  -0.344422,   0.289871,   0.287706,   0.078915,  -1.152880,   0.165900,   0.547004,  -0.115703,  -0.253022,   0.026322,   0.139885,  -0.908586,   0.204839,   0.298074,   0.086595,   0.346151,   0.682055,   0.185848,   0.220974,  -0.015863,   0.017416,  -0.014225,   0.308962,   0.147927,
        -2.223890,   2.433890,  -0.280859,  -0.106824,  -0.037261,  -0.390541,  -0.214066,  -0.007258,  -0.155326,   0.092309,  -0.928678,   0.076390,   0.364377,  -0.128957,  -0.177499,   0.053749,   0.158252,   0.630359,   0.004089,  -0.130998,  -0.218170,   0.029844,   0.119839,   0.020472,   0.329279,   0.201266,  -0.240405,   0.133646,   0.165527,  -0.013415,   0.013907,  -0.339297,
        -0.689335,   0.051509,   0.072654,  -0.086975,  -0.057655,   0.020286,  -0.243857,  -0.067477,  -0.131258,  -0.154884,  -0.106722,  -0.073534,  -0.147594,   0.191870,  -0.072761,   0.194225,   0.056835,   0.087386,   0.032827,   0.559167,  -0.098282,   0.266647,   0.154815,  -0.033550,  -0.045125,  -0.067090,   0.098506,   0.113227,   0.176634,  -0.175096,  -0.067111,   0.123902,
         0.127824,   0.250559,  -0.211324,   0.235335,   0.064369,  -0.226249,   0.410887,   0.070645,  -0.445704,  -0.256447,  -0.112341,  -0.068689,  -0.045278,  -0.254021,   0.142984,  -0.043897,   0.025515,  -0.011799,  -0.268901,   0.181901,   0.142316,   0.018397,  -0.154069,  -0.049348,   0.195062,   0.004141,  -0.261675,  -0.116848,   0.236352,  -0.046884,   0.033339,   0.542600,
         0.214734,   0.332746,   0.091315,  -0.167167,   0.607484,  -0.006371,  -0.500489,  -0.048774,   0.050059,   0.044115,  -0.241985,   0.250713,   0.012495,   0.188359,   0.206148,  -0.128451,  -0.586723,   0.010017,   0.084768,  -0.000395,   0.088979,   0.032891,  -0.001245,   0.074544,  -0.215994,  -0.386776,   0.159382,  -0.070491,   0.177923,   0.239460,   0.259041,  -0.068195,
         0.258927,  -0.114790,   0.500049,   0.160040,  -0.408526,  -0.061467,   0.376828,   0.023297,  -0.243852,   0.237101,  -0.097912,   0.018360,   0.072743,  -0.015845,  -0.431004,   0.152441,   0.235947,   0.100826,   0.012485,   0.330335,  -0.039215,  -0.009020,   0.185518,  -0.188768,  -0.288554,  -0.027179,   0.048242,   0.108480,  -0.028293,   0.034381,   0.353483,  -0.074875,
         0.290451,   0.161438,  -0.517587,   0.025542,   0.138566,   0.042365,  -0.109223,   0.165244,  -0.113960,   0.062502,   0.158794,   0.141409,  -0.366148,  -0.284597,   3.160930,  -0.143111,
        -0.281299,  -0.006561,   0.036973,  -0.199692,   0.095575,  -0.091180,  -0.150344,  -0.005043,   0.157037,  -0.933489,  -0.131437,   0.088807,  -0.087258,   0.145750,  -0.171973,  -0.166491,   0.147847,   0.134252,   0.521335,   0.105440,   0.003666,   0.010571,  -0.106024,  -0.160045,   0.014484,   0.105589,   0.131343,   0.051944,   0.106205,   0.033912,  -0.051948,   0.352421,
         0.168851,   0.140763,   0.293494,  -0.018893,   0.158293,   0.109333,   1.154120,  -0.198061,  -0.721676,   0.241752,  -0.140345,  -0.360857,  -0.388876,   0.168425,  -0.100143,  -0.043754,  -0.310466,  -0.015447,   0.183322,   0.293166,  -0.325573,  -0.031339,   0.218661,   0.011904,   0.139696,  -0.104750,  -0.045275,   0.378521,   0.022482,   0.422124,   0.304682,  -0.051797,
         0.073359,  -0.298343,   0.343329,   0.227006,   0.894695,  -0.312772,  -0.632907,   0.156611,  -0.195842,   0.057613,  -0.518628,   0.264291,   0.081383,  -0.077500,  -0.127376,  -0.007693,   0.122776,   0.033892,  -0.260476,   0.169595,   0.296670,  -0.718873,   0.270009,   0.130217,   0.118807,   0.184551,  -0.053275,  -0.047685,   0.162206,  -0.097828,  -0.029592,  -0.101494,
         0.016215,   0.083746,   0.615458,  -0.485485,  -0.921246,   0.084868,   0.028433,   0.257626,  -0.031143,  -0.212012,   0.167562,  -0.051356,   0.172074,   0.229712,  -0.194711,   0.025840,  -0.219482,  -0.111193,   0.049931,  -0.288035,   0.245568,   0.105462,  -0.091076,  -0.023110,  -0.254213,  -0.097878,  -0.035360,  -0.200398,  -0.015852,  -0.110506,   0.033963,  -0.044191,
         0.924467,  -0.258591,  -0.671682,  -0.023756,   0.013241,  -0.156469,  -0.304091,  -0.071585,  -0.071086,  -0.084754,   0.224128,  -0.180669,  -0.196715,  -0.093340,  -0.254068,   0.171081,  -0.024239,  -0.252692,   0.049843,   0.000023,   0.033901,  -0.130915,  -0.195773,   0.179232,   0.128790,  -0.064079,   0.016506,  -0.111491,  -0.100154,  -0.030905,  -0.849884,  -1.346140,
         0.788281,   0.137688,  -0.021214,   0.083188,   0.091056,  -0.081960,   0.076624,   0.058724,   0.535560,  -0.031969,   0.081119,  -0.025634,  -0.048527,   0.020818,   0.080830,  -0.192909,   0.061279,  -0.001361,   0.051785,   0.005867,  -0.115500,   0.295237,   0.009252,   0.045897,  -0.024886,   0.023726,   0.179642,   0.016852,  -0.365523,  -0.102476,  -0.531131,   0.250517,
         0.225560,   0.067154,   0.000948,   0.072650,   0.137132,   0.290897,  -0.144045,   0.057894,   0.109130,   0.009964,   0.094089,  -0.154842,   0.113549,   0.078136,   0.249298,   0.110934,   0.048162,   0.048931,   0.211603,   0.534025,   0.322002,   0.097533,  -0.278491,  -0.078031,   0.605147,   0.340444,  -0.295375,  -0.183762,  -0.414944,   0.088145,   0.178964,   0.133147,
        -0.106853,   0.163846,   0.272785,   0.352460,  -0.550566,  -0.010568,   0.154589,  -0.023544,   0.171297,   0.068605,   0.028248,  -0.564806,   0.400015,   0.410908,   0.014326,   0.289611,   0.288690,   0.142184,   0.226587,  -0.058048,   0.160932,  -0.123207,   0.194131,   0.187663,  -0.742981,   0.025905,  -0.846101,   0.025712,   0.039420,   0.105765,  -0.079089,   0.198622,
         0.045839,   0.010623,  -0.114428,  -0.098091,   0.212623,   0.018271,  -0.211084,   0.155203,   0.056616,  -0.529915,   0.225675,   0.066716,   0.203242,   0.185543,   0.039155,   0.198155,   0.142782,  -0.057852,   0.011131,  -0.093834,   0.125661,   0.101441,  -0.874059,   0.131259,  -0.861177,  -0.066185,  -0.184366,   0.152246,  -0.275628,  -0.054194,  -0.076954,   0.030010,
        -0.220138,  -0.129643,   0.157565,   0.013654,  -0.108522,   0.133376,   0.132067,  -0.131420,  -0.160317,   0.048937,   0.004411,   0.120882,   0.041790,   0.953594,  -0.706424,  -0.146763,
        -0.224518,   0.265232,   0.208360,  -0.002020,   0.073278,  -0.059616,  -0.050920,   0.086407,  -0.966636,  -1.275920,   0.992438,   0.152685,   0.139583,   0.276266,   0.167204,   0.027651,   0.258806,   0.421207,   0.261805,  -0.070224,   0.103280,   0.061884,   0.070232,  -0.110944,  -0.082762,   0.331400,   0.183056,  -0.033889,  -0.113225,   0.054178,   0.024335,   0.178651,
         0.103579,  -0.132714,   0.178110,  -0.079222,  -0.116105,   0.047456,  -0.540050,  -0.448129,   1.014950,  -0.032571,   0.337298,   0.719880,   0.754488,  -0.278969,   0.184178,   0.508169,   0.807849,   0.070604,  -0.080654,  -0.123654,   0.072679,   0.031160,  -0.034415,   0.484169,   0.222930,   0.075113,   0.104601,   0.194068,  -0.283239,   0.399475,   0.140406,  -0.042276,
        -0.111286,  -0.147011,   0.460438,   0.319647,  -0.605285,  -0.144382,   0.236490,   0.009905,   0.454188,   0.666478,   0.309631,  -0.073454,   0.440559,   0.336458,   0.367795,   0.186402,  -0.066336,  -0.036435,   0.508142,   0.068927,   0.023296,   0.224636,   0.556337,   0.209546,  -0.085867,   0.137131,  -0.194997,   0.208912,  -0.011895,  -0.130697,  -0.154590,  -0.042808,
        -0.059538,   0.096543,  -0.688514,  -0.045219,  -0.626087,   0.064348,   0.117739,   0.256758,   0.194783,  -0.023289,   0.230582,   0.011195,   0.079812,  -0.028160,  -0.013248,  -0.188248,  -0.010532,  -0.118205,  -0.048721,   0.297909,   0.181207,   0.294900,  -0.040897,  -0.133930,  -0.110695,   0.110513,   0.022752,  -0.210805,  -0.104696,   0.025128,  -0.098457,  -0.028057,
        -0.363816,   0.385634,  -0.798998,   0.054868,   0.001373,   0.065832,  -0.134640,  -0.015175,   0.083533,   0.026305,   0.213719,   0.009375,  -0.144116,  -0.175280,  -0.087167,   0.044839,  -0.135140,   0.117899,   0.039121,  -0.099187,   0.099506,  -0.020557,   0.029771,   0.172305,   0.099862,  -0.020464,  -0.084566,  -0.146625,  -0.074283,   0.030666,  -0.649815,  -0.689195,
         0.477989,  -0.112364,   0.259021,  -0.147118,  -0.160496,  -0.015536,   0.160702,  -0.118447,   0.352424,   0.231288,   0.351925,   0.192150,   0.020422,   0.136044,  -0.033926,  -0.226486,   0.036423,   0.316355,  -0.052222,   0.198901,   0.104821,   0.303561,   0.127232,  -0.023316,   0.019414,  -0.059270,   0.383971,   0.160645,   0.082094,  -0.048621,  -0.787037,   0.254908,
        -0.371708,  -0.277290,  -1.059800,   0.154760,  -0.060394,  -0.060609,  -0.451491,  -0.178901,   0.309639,   0.249477,  -0.422599,   0.122808,   0.233100,   0.016934,   0.062681,   0.018490,  -0.051389,   0.064288,   0.258257,   0.181593,   0.249832,  -0.074327,  -0.140494,  -0.228158,   0.731626,   0.253385,   0.199566,  -0.351042,  -0.486372,   0.227417,  -0.254663,  -0.181428,
        -0.521219,   0.362583,   0.147322,   0.146914,  -0.259536,  -0.052050,   0.216003,   0.097527,  -0.216906,   0.274866,  -0.005636,  -1.357990,   0.344849,   0.113716,   0.085554,   0.111638,   0.182135,   0.043280,   0.152514,  -0.119686,   0.053910,  -0.190329,   0.250407,   0.111022,  -0.050444,  -0.275560,  -0.914624,   0.049103,  -0.227758,   0.080444,  -0.225389,   0.079419,
        -0.119264,  -0.041744,   0.265803,  -0.041152,  -0.034882,  -0.073968,  -0.577530,   0.000363,   0.014838,  -0.645189,   0.190927,  -0.094755,  -0.173985,  -0.003143,  -0.195717,  -0.186654,  -0.010083,  -0.089036,  -0.041747,   0.077850,   0.050264,  -0.009428,  -0.003300,   0.054101,  -0.887406,   0.057632,  -0.277544,  -0.053275,  -0.583980,   0.279462,  -0.086696,  -0.220415,
         0.095657,  -0.222685,   0.063930,   0.084972,  -0.598213,   0.375499,   0.047898,  -0.857509,   0.008648,  -0.255603,  -0.235206,  -0.177106,  -0.158840,   0.353964,   2.285770,   0.353526
    };

    //! weights for layer 2
    const double ANN_WEIGHTS_CONTACT_STRAND_STRAND_LAYER_1[ 17] =
    {
        -0.259087,  -0.597632,  -0.693580,   0.959494,   0.735067,  -0.623302,   0.769421,  -0.850094,  -0.750768,  -0.595051,  -0.684021,   1.640890,  -0.620750,   0.651551,   0.787905,   0.490276,   0.673933
    };
    //! ANN CONTACT_STRAND_STRAND definition
    double ANN_CONTACT_STRAND_STRAND( const linal::Vector< double> &INP)
    {

      // declare variables
      int nora( 0), norb( 0), wei( 0);
      double *hid;

      // test net size
      BCL_Assert( INP.GetSize() == 303, "wrong input size!");

      // allocate memory
      linal::Vector< double> hidden[ 3];
      hidden[0] = linal::Vector< double>( 303);
      hidden[1] = linal::Vector< double>( 16);
      hidden[2] = linal::Vector< double>( 1);

      // normalize data
      hid = hidden[ 0].Begin();
      for( const double *inp = INP.Begin(); inp != INP.End(); inp++)
        ( *( hid++)) = ( *inp) * ANN_NORMALIZE_CONTACT_STRAND_STRAND_A[ nora++] + ANN_NORMALIZE_CONTACT_STRAND_STRAND_B[ norb++];

      // calculate network
      // calculate layer 1
      wei = 0;
      hid = hidden[1].Begin();
      for( size_t i = 0; i < 16; i++)
      {
        *hid = ANN_WEIGHTS_CONTACT_STRAND_STRAND_LAYER_0[ wei++];
        for( const double *inp = hidden[ 0].Begin(); inp != hidden[ 0].End(); inp++)
          ( *hid) += ANN_WEIGHTS_CONTACT_STRAND_STRAND_LAYER_0[ wei++] * ( *inp);
        *hid = double( 1.0) / ( double( 1.0) + exp( -( ( *hid)))); hid++;
      }

      // calculate layer 2
      wei = 0;
      hid = hidden[2].Begin();
      for( size_t i = 0; i < 1; i++)
      {
        *hid = ANN_WEIGHTS_CONTACT_STRAND_STRAND_LAYER_1[ wei++];
        for( const double *inp = hidden[ 1].Begin(); inp != hidden[ 1].End(); inp++)
          ( *hid) += ANN_WEIGHTS_CONTACT_STRAND_STRAND_LAYER_1[ wei++] * ( *inp);
        *hid = double( 1.0) / ( double( 1.0) + exp( -( ( *hid)))); hid++;
      }

      // denormalize data
      // end
      return hidden[ 2]( 0) * ANN_NORMALIZE_CONTACT_STRAND_STRAND_A[ nora] + ANN_NORMALIZE_CONTACT_STRAND_STRAND_B[ norb];

    }

    //! normalization values A*x+b
     const double ANN_NORMALIZE_CONTACT_SHEET_SHEET_A[ 304] =
    {
       0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 2, 0.0909091, 0.2, 0.1, 2, 2, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.000833333, 0.000833333, 0.000833333, 1.2
    };

    //! normalization values a*x+B
     const double ANN_NORMALIZE_CONTACT_SHEET_SHEET_B[ 304] =
    {
       0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0.0909091, 0.4, -0.2, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, -0.1
    };

    //! weights for layer 1
    const double ANN_WEIGHTS_CONTACT_SHEET_SHEET_LAYER_0[ 4864] =
    {
         0.122718,  -0.118309,  -0.147300,   0.014263,  -0.181864,   0.066530,   0.010351,  -0.037669,   1.087340,   0.144754,   0.035546,  -0.215915,   0.105573,  -0.182621,   0.104643,   0.053121,   0.140362,   0.251377,   0.124585,   0.165302,  -0.355437,  -0.223225,   0.155097,   0.030084,  -0.328546,  -0.308772,   0.131946,   0.009910,  -0.062275,  -0.051309,  -0.228476,  -0.030613,
         0.168859,   0.093109,   0.391972,   0.209663,  -0.000572,  -0.019519,   0.249348,   0.078009,   0.611711,   0.125114,  -0.160679,  -0.491908,  -0.112101,  -0.195430,  -0.182908,  -0.267822,   0.136741,  -0.280221,   0.063100,   0.202421,   0.014719,   0.209174,   0.056753,  -0.499096,  -0.375415,  -0.251113,   0.623320,  -0.237423,   0.093490,  -0.155204,   0.175256,   0.353772,
        -0.280925,   0.141355,   0.063361,   0.060201,   0.106729,   0.137995,   0.506555,   0.080015,   0.651756,  -0.094950,  -0.262336,   0.141306,   0.272638,   0.415911,   0.238234,   0.505538,  -0.285994,  -0.455326,   0.438962,   0.056848,  -0.084034,  -0.175487,   0.292554,   0.013021,   0.899302,   0.294849,   0.016500,   0.091405,  -0.068361,  -0.055819,   0.492057,   0.116813,
         0.053239,   0.127136,   0.143231,  -0.001637,   1.101860,   0.116889,   0.241640,  -0.277453,   0.313327,   0.201529,  -0.209047,  -0.236925,  -0.064270,   0.076948,  -0.209384,  -0.080815,   0.079663,  -0.126305,  -0.219694,  -0.034642,  -0.043551,  -0.070386,  -0.022052,  -0.090952,  -0.205973,   0.082754,   0.073838,   0.191788,   0.123899,   0.119535,  -0.081878,   0.050274,
         0.358287,  -0.124634,   0.997725,  -0.094943,   0.240790,  -0.448299,   0.155090,  -0.112818,  -0.211023,   0.012855,  -0.168679,   0.172577,  -0.062904,  -0.010081,   0.126513,  -0.196279,  -0.057592,   0.421429,  -0.401233,  -0.203023,   0.040449,   0.175990,  -0.012554,  -0.044640,  -0.155630,   0.012565,   0.064906,   0.317504,   0.085508,  -0.028329,   0.433831,   0.011980,
         0.185706,   0.045024,   0.615557,  -0.054364,  -0.257052,   0.053290,  -0.039314,   0.171635,   0.369262,  -0.045687,   0.118442,   0.134080,   0.381905,   0.073510,  -0.251408,  -0.001634,  -0.161746,  -0.239677,  -0.203136,  -0.307253,  -0.013862,   0.002751,   0.089421,   0.270436,   0.132414,   0.103077,   0.002548,   0.001109,   0.106146,   0.023269,   0.816986,  -0.107773,
         0.084457,  -0.297946,  -0.224926,  -0.065234,  -0.234498,  -0.133294,  -0.234999,   0.129679,  -0.001860,   0.342445,   0.072018,   0.132791,   0.119710,  -0.224735,  -0.392891,  -0.159304,  -0.066519,  -0.063165,   0.321512,   0.243614,  -0.026932,   0.266955,  -0.295149,   0.237951,   0.141464,   0.041859,   0.317058,   0.171910,   0.456101,   0.040369,   0.485827,   0.011494,
        -0.201613,   0.132421,   0.143732,   0.576202,   0.226348,   0.477902,  -0.151190,  -0.356203,   0.419406,  -0.074991,  -0.329187,   0.137390,   0.195794,   0.046295,  -0.185741,   0.063741,   0.001254,   0.113537,   0.180463,   0.076021,   0.223607,   0.067392,   0.108467,   0.117735,   0.264057,   0.130103,   0.446837,  -0.079405,   0.220188,   0.107729,   0.033809,  -0.038684,
        -0.319674,   0.039536,  -0.069595,   0.263597,   0.332632,  -0.039217,   0.116254,  -0.059285,  -0.150447,  -0.299521,  -0.134588,  -0.099102,   0.025598,  -0.095487,   0.050990,   0.006203,  -0.054307,   0.138390,  -0.021269,   0.111344,  -0.014295,   0.080075,   0.468690,  -0.201212,   1.233180,   0.130355,   0.091597,  -0.494745,   0.185264,  -0.124949,   0.094391,   0.182586,
        -0.215415,   0.056121,  -0.182194,  -0.025272,   0.122287,  -0.020536,  -0.252295,   0.406764,   0.236289,   0.261483,   0.381436,   0.108080,   0.027216,   0.547453,   1.299430,   1.243210,
         0.113892,  -0.166620,  -0.002125,   0.003904,  -0.060984,   0.108502,  -0.134155,  -0.140368,   0.310584,   0.186225,  -0.228620,  -0.265727,  -0.125297,  -0.005857,   0.024414,   0.120551,   0.008676,   0.091717,  -0.407501,  -0.121941,  -0.254033,   0.069864,  -0.076016,  -0.091533,   0.057211,  -0.404861,   0.067299,   0.097278,  -0.299835,  -0.052058,   0.056566,  -0.030657,
         0.099166,  -0.016183,   0.108029,  -0.029472,  -0.049422,   0.025088,   0.144413,   0.103495,   0.299651,   0.026640,  -0.130710,   0.803311,   0.814718,  -0.033410,  -0.008014,  -0.043746,  -0.743366,   0.312469,  -0.436106,  -0.101319,  -0.113362,  -0.272592,  -0.108235,  -0.699600,  -0.095513,   0.002338,   0.339518,   0.086506,  -0.388639,   0.313133,   0.791838,   0.493186,
        -0.285218,   0.020525,   0.506759,   0.165517,  -0.004382,   0.058450,  -0.017445,   0.030171,   0.764525,   0.179091,  -0.399591,   0.201662,   0.059886,   0.320922,  -0.789358,   0.176190,  -0.236124,   0.023755,   0.551591,  -0.251609,  -0.198882,  -0.089746,   0.289230,   0.037528,   0.314935,   0.167679,   0.016144,  -0.027595,   0.161883,   0.098772,  -0.027146,  -0.026170,
        -0.098101,   0.111103,   0.359249,  -0.180055,   0.893961,  -0.193350,   0.206980,   0.567402,   1.465710,  -0.149248,  -0.167589,   0.694094,   0.272119,   0.145760,  -0.246769,  -0.103950,   0.853974,  -0.149765,  -0.294408,  -0.304962,   0.044927,   0.096674,  -0.098150,  -0.112378,  -0.269178,   0.220227,   0.016771,   0.027667,   0.083986,   0.036213,  -0.025487,   0.129322,
         0.688036,  -0.179312,   0.459638,   0.018054,   0.063437,   0.067752,   0.070120,   0.130634,  -0.067352,   0.080089,   0.138439,  -0.757603,   0.955483,   0.296264,   0.374451,   0.056154,  -0.005245,   0.178209,  -0.219815,   0.436104,   0.389455,   0.145465,   0.572784,   0.043128,   0.082793,   0.085746,  -0.060450,   0.084195,   0.007227,   0.070916,   0.323102,   0.133447,
        -0.189315,  -0.206958,   0.195243,   0.120039,  -0.313627,   0.175837,   0.108307,  -0.284317,  -0.064936,  -0.035512,  -0.016238,   0.191386,   0.150126,  -0.025114,  -0.017313,  -0.985632,   0.004146,   0.036273,  -0.247275,   0.102692,   0.038266,  -0.148611,   0.017989,  -0.068720,  -0.157971,  -0.024125,  -0.145993,   0.086051,   0.038773,   0.035386,   0.662530,  -0.112019,
         0.223194,   0.080444,   0.696041,  -0.155619,   0.098158,   0.354176,  -0.427305,  -0.399808,  -0.480714,  -0.179265,   0.267809,  -0.032275,  -0.159854,  -0.710304,  -0.054289,  -0.020477,  -0.226996,  -0.150155,  -0.233319,   0.616708,   0.615717,   0.590639,  -0.198952,   0.146939,   0.480397,   0.046959,   0.052203,   0.101803,   0.423335,   0.338681,   0.834734,  -0.877082,
        -0.149640,  -0.093919,   0.049644,   0.251774,  -0.595605,   0.234181,  -0.384260,  -0.114258,   0.374185,  -0.212290,  -0.169050,  -0.044164,   0.243427,   0.452651,   0.018194,   0.290399,  -0.165623,  -0.073608,  -0.057252,   0.041423,  -0.088145,   0.138230,  -0.003342,   0.085129,   0.297889,   0.039418,   0.701557,  -0.097260,   1.310210,   0.667163,   0.366303,  -0.587019,
         0.364938,  -0.089020,   0.272604,   0.225590,  -0.419742,  -0.216318,   0.498006,  -0.291605,  -0.249682,   0.605564,   0.113044,  -0.072302,   0.055746,  -0.234568,  -0.226351,   0.104281,  -0.019152,   0.063796,   0.050178,   0.048140,  -0.051739,  -0.006585,   0.399368,  -0.050410,   0.606213,   0.127349,   0.390785,  -0.030716,  -0.202797,  -0.058083,  -0.120887,  -0.128285,
        -0.015666,  -0.217501,   0.329450,   0.663755,   0.147264,   0.125313,  -0.021680,  -0.097796,  -0.119072,   0.101333,  -0.038989,   0.035788,   0.707739,  -0.043722,   0.532145,  -0.047309,
         0.086531,   0.025165,  -0.085065,  -0.061302,   0.165429,   0.214122,   0.019577,   0.099869,   0.676443,   0.520203,  -0.455272,   0.200748,  -0.376923,  -0.076223,  -0.397541,   0.215072,  -0.357119,  -0.364345,   0.211734,  -0.016001,   0.081514,   0.231135,  -1.135930,   0.134705,  -0.237314,  -1.252740,  -0.047669,  -0.024698,   0.086415,  -0.149644,   0.189214,   0.070943,
        -0.035752,  -0.006282,   0.548948,   0.053558,   0.045080,   0.083946,   0.370819,   0.102943,   0.270768,   0.159317,  -0.401558,  -0.169430,   0.281559,   0.520958,  -0.204989,  -0.507001,   0.527212,  -0.091267,   0.225596,   0.158521,  -0.420294,   0.739435,   0.031825,  -1.585800,   0.034214,  -1.040170,  -0.490909,  -0.090936,   0.069913,   0.018912,   0.308446,   0.421012,
        -0.181730,   0.167410,   0.048985,   0.115445,   0.236539,   0.675607,  -0.157740,  -0.271171,   0.500667,   0.560010,   0.819485,  -0.150353,   0.303706,   0.778357,  -0.079872,   0.858612,  -0.296469,  -0.367959,   0.751480,  -0.039344,  -0.139128,  -0.001482,   0.432011,  -0.136241,   0.397322,   0.262539,  -0.035212,  -0.089059,  -0.228065,  -0.043806,   0.632225,  -0.022712,
        -0.154237,   0.172377,   1.066350,   0.353067,   0.245734,   0.194370,  -0.280857,  -0.077625,  -0.550971,   1.510380,  -0.406584,   0.009061,   0.269467,   0.029190,   0.161281,   0.274516,   0.032016,   0.230108,   0.153173,  -1.031990,  -0.057432,  -0.231810,  -0.389511,  -0.056163,   0.050170,   0.164656,   0.059556,   0.111394,   0.245976,   0.053996,   0.033241,   0.149718,
         0.495133,   0.787229,  -0.140400,   0.059833,  -0.304758,  -0.096888,  -0.580828,  -0.054306,   0.048701,   0.117968,  -0.014085,  -0.126584,   0.105292,  -0.098071,   0.062014,  -0.069623,   0.419869,  -0.451809,   0.072450,   0.245023,   0.284395,   0.359621,   0.148468,   0.187359,   0.022660,   0.180065,   0.151720,   0.160015,   0.239703,   0.308184,   0.914226,   0.413738,
        -0.557565,   0.240191,   0.064307,  -0.052184,  -0.070744,   0.600563,  -0.209309,  -0.131855,  -0.578809,  -0.091058,   1.470620,   0.510278,   0.088247,   0.194783,   0.866747,  -0.212280,   0.124451,  -0.002805,  -0.141293,  -0.005426,   0.626396,  -0.018534,  -0.113222,   0.021501,   0.275523,   0.139639,  -0.089101,   0.070347,   0.467236,   0.247755,   0.311162,  -0.125311,
         0.131499,   0.145067,   0.338519,  -2.754690,  -0.197020,  -0.041057,   0.533938,   0.147000,  -0.711577,   0.183636,  -0.100866,  -0.162468,  -0.487664,   0.084915,  -0.120294,  -0.114458,  -0.139195,  -0.164725,  -0.017525,   0.109898,  -0.645127,   0.217377,  -0.162293,   0.194617,   0.078174,   0.182511,   0.346837,   0.491833,  -0.020737,   0.060000,   0.200938,   0.570798,
        -1.182010,   0.208566,   0.193385,   0.646340,  -0.472866,  -0.037229,  -0.103885,  -0.234758,   0.318139,   0.039284,  -0.920294,   0.178623,   0.313520,   0.266544,  -1.734640,  -0.326826,   0.041551,   0.069548,   0.241462,   0.186183,   0.195550,   0.098967,   0.105902,  -0.056579,   0.708470,   0.318520,   0.381164,  -0.072152,   0.293934,   0.077715,   0.126633,  -0.320713,
        -0.705700,   0.592640,   0.154192,   0.909863,   0.301053,   0.009016,   0.290249,  -0.007114,  -0.068385,   0.155262,   0.078534,   0.043611,   0.573880,   0.024813,  -0.022045,   0.155494,   0.069896,   0.104820,  -0.033926,   0.249806,   0.053597,   0.043357,  -0.409945,   0.317262,   0.533045,  -0.017192,  -0.101537,   0.113058,  -0.350924,   0.160507,   0.090470,   0.033121,
         0.481918,  -0.145812,  -0.848828,  -0.772977,  -0.199443,   0.047381,   0.149259,   0.310655,   0.045703,   0.225728,   0.696279,   0.288975,  -0.598369,   0.512159,   1.423630,   0.912052,
        -0.143326,   0.634699,   0.155786,  -0.059172,  -0.002051,  -0.043126,   0.141272,   0.032673,  -0.738866,  -0.738802,   0.719357,  -0.101835,   0.170135,  -0.187143,  -0.424769,   0.193949,  -0.186173,  -0.024114,   0.269310,  -0.413594,   0.157551,   0.362903,   0.091873,   0.168692,   0.420642,  -0.405366,  -0.043125,   0.276506,  -0.076705,  -0.024378,   0.168846,   0.001490,
         0.217424,   0.135232,  -0.003907,   0.117659,   0.693385,   0.054653,  -0.513763,   0.702072,  -0.363442,   0.078571,   0.120234,  -0.491037,  -0.663762,   0.290406,   0.132538,   0.235639,  -0.138113,   0.223911,   0.521284,   0.363104,   0.410908,   0.028152,   0.180445,  -0.021083,   0.159930,  -0.002509,   0.013730,   0.161350,   0.256761,  -0.160374,   0.726925,  -0.135146,
         0.891986,   0.371536,   0.646573,  -0.135597,  -0.932178,   0.697987,  -0.483147,   0.436644,  -0.349081,  -1.175310,  -1.662190,   0.402898,  -0.193366,  -0.363772,  -0.390279,  -0.528058,   0.643937,   0.263274,  -0.382183,   0.190027,  -0.018610,  -0.992355,  -0.913448,  -0.348652,   0.293125,  -0.251161,   0.647462,   0.335616,   0.142877,   0.013956,  -0.137903,  -0.228103,
         0.518394,   0.172003,  -1.280850,   1.990330,  -1.064990,   0.362199,  -0.333923,  -0.655125,  -0.386930,   0.646147,  -0.031077,  -0.287930,  -0.773868,   0.119548,   0.467504,  -0.077462,  -0.372899,  -0.072437,  -0.011174,  -1.105650,   0.219867,   0.226280,  -0.008313,   0.150410,   0.761330,   0.228073,   0.229392,   0.073168,   0.175855,   0.057032,   0.219510,  -0.104993,
        -2.274560,   1.279490,  -0.082465,   0.054356,  -0.359827,  -0.314061,  -0.826565,   0.420730,  -0.445975,  -0.473687,  -0.345354,   0.119629,   0.163277,   0.346404,  -0.573355,   0.247447,   0.076045,   0.028426,  -0.198396,  -0.352979,   0.048842,  -0.297036,   0.303157,   0.184124,  -0.026089,  -0.003216,   0.122545,   0.176310,  -0.037283,   0.043659,  -0.584154,  -0.154353,
        -0.045430,   0.049805,   0.120502,   0.630623,  -0.109367,  -0.077121,   0.625897,   0.696678,   0.147735,  -0.030242,  -0.065788,  -0.185781,  -0.119650,   0.212909,   0.112568,  -0.043677,   0.274019,   0.130820,   0.047040,   0.142000,   0.054380,   0.086675,   0.318576,   0.271225,   0.003895,  -0.045218,   0.228412,   0.059239,   0.015463,  -0.098487,  -0.243557,  -0.050241,
        -0.360695,   0.235898,   0.629149,  -0.098449,   0.043742,   0.087611,  -0.160849,  -0.214652,  -0.120090,   0.015440,  -0.315651,   0.062968,  -0.189395,  -0.091745,   0.060849,   0.129418,  -0.071577,   0.152121,  -0.167843,   0.414571,   0.092836,  -0.110393,   0.086718,   0.025887,  -0.149231,   0.146788,  -0.091652,  -0.112425,   0.344779,  -0.229954,   1.221600,   1.096950,
         0.675850,  -0.003027,   0.651028,   0.509431,  -0.240225,   0.909574,  -0.433207,  -0.304706,   0.643004,   0.041114,  -0.307621,   0.371096,   0.103016,   0.375825,   0.247955,   0.583547,  -0.160749,  -0.082487,  -0.083975,   0.035784,  -0.218310,   0.292686,   0.128879,  -0.009873,  -0.318190,  -0.271916,   0.531516,  -0.055900,   0.121764,   0.056929,  -0.014020,   0.005480,
         0.019554,  -0.287522,   0.280307,  -0.301294,  -0.354488,   0.096044,  -0.109964,  -0.209750,  -0.089437,   0.335376,  -0.196848,  -0.076303,  -0.012859,  -0.123174,  -0.146202,  -0.055743,  -0.083387,  -0.092832,   0.023043,  -0.136993,  -0.031888,  -0.071138,  -0.314748,  -0.102568,   0.245048,   0.019365,   0.095763,   0.576813,   0.192353,  -0.228032,   0.268332,  -0.000358,
         0.034565,  -0.133190,  -0.138942,  -0.022366,   0.236459,   0.560292,  -0.050614,   0.089034,  -0.100030,   0.211808,   0.183568,   0.444001,  -0.137016,  -0.280070,   1.819860,   0.372559,
        -0.289516,  -0.008040,   0.043520,  -0.007427,  -0.053493,  -0.147815,   0.125041,  -0.054696,  -0.348890,  -0.902152,   0.038322,  -0.057638,   0.176901,  -0.003905,  -0.093296,   0.137462,  -0.005062,   0.065508,   0.323576,   0.100330,  -0.067853,  -0.167531,   0.148299,   0.014584,   0.192182,   0.010027,  -0.018999,  -0.003488,   0.109403,   0.045031,  -0.107720,   0.124185,
         0.095971,  -0.004743,  -0.243003,  -0.100991,   0.063412,  -0.020680,   0.111601,  -0.214626,  -0.453779,  -0.004932,   0.171711,   0.035370,  -0.417674,   0.013666,   0.079548,   0.449640,  -0.347889,   0.034935,   0.279005,   0.141403,   0.207249,  -0.057412,  -0.074031,   0.318646,   0.328367,   0.407998,   0.062876,   0.059032,   0.185555,  -0.197680,   0.108296,  -0.223863,
         0.588538,   0.092158,  -0.117296,  -0.083920,   0.162349,  -0.201002,  -0.133191,   0.103249,  -0.302183,  -0.342076,  -0.324779,  -0.027743,  -0.048462,  -0.194805,  -0.126408,  -0.098801,   0.360701,   0.253104,  -0.443100,   0.228046,   0.270691,  -0.128121,  -0.260572,  -0.001517,  -0.106436,  -0.011736,   0.152667,   0.109127,   0.013795,  -0.169161,  -0.313737,  -0.050038,
         0.174617,   0.001337,   0.064968,  -0.114206,  -0.585565,   0.044059,   0.110533,   0.127356,   0.005418,   0.031953,   0.252237,   0.070503,   0.012277,   0.170360,   0.172792,  -0.017417,   0.061006,  -0.072786,   0.011593,  -0.161231,   0.114909,   0.144103,   0.088359,   0.210598,   0.024955,   0.182822,   0.137559,   0.003892,  -0.053203,  -0.008931,   0.109990,  -0.070103,
        -0.303586,  -0.426087,  -0.051349,   0.056161,  -0.092402,   0.277227,  -0.209018,  -0.049998,  -0.055041,   0.051233,   0.030854,  -0.046999,  -0.067997,   0.092492,  -0.179109,   0.093633,  -0.061630,   0.108010,   0.111070,  -0.050664,  -0.094639,   0.071635,  -0.073594,   0.101648,   0.226431,  -0.002100,  -0.097819,  -0.076226,   0.121529,   0.010668,  -0.680431,  -0.840808,
         0.209721,   0.032127,   0.112083,   0.007422,  -0.091563,  -0.037871,   0.202256,   0.071175,   0.088813,  -0.253232,   0.110156,  -0.061830,  -0.008606,   0.060975,   0.137448,  -0.037631,  -0.085774,   0.194513,   0.001033,   0.126075,   0.322331,   0.170012,  -0.026089,   0.036827,  -0.101090,   0.017737,   0.024196,  -0.015680,  -0.165268,  -0.293746,  -0.242736,   0.213847,
         0.110645,  -0.128286,  -0.297462,  -0.018032,   0.181674,   0.110121,  -0.284256,  -0.176866,   0.127199,  -0.119124,   0.034523,  -0.027174,   0.011166,   0.219072,   0.268526,   0.111530,  -0.001414,   0.068803,   0.147602,  -0.136104,  -0.245284,  -0.417857,   0.701233,  -0.156488,  -0.112284,  -0.117758,   0.152803,  -0.183277,  -0.037741,   0.119020,  -0.256202,  -0.222096,
        -0.972831,   0.086233,   0.057937,  -0.251443,   1.284130,  -0.122685,   0.377660,   0.346870,  -0.386460,   0.046356,   0.190318,  -0.109334,  -0.389733,  -0.337933,  -0.077652,   0.134416,   0.505180,   0.035807,   0.037102,  -0.149649,  -0.250098,  -0.093544,   0.138933,  -0.028229,  -0.622676,  -0.108973,  -0.286758,   0.228677,  -0.093505,   0.181243,   0.051384,   0.196614,
         0.182154,   0.025652,   0.137389,  -0.032029,   0.016885,  -0.090726,  -0.026646,  -0.024823,  -0.145487,  -0.067705,   0.179986,   0.158073,   0.029906,  -0.020342,   0.033565,  -0.027156,  -0.003226,  -0.133369,  -0.056674,  -0.143931,   0.119544,  -0.173578,  -0.101098,  -0.315931,  -0.269096,   0.067710,  -0.051500,   0.434677,   0.020234,   0.014474,  -0.020783,  -0.083215,
        -0.064134,   0.086647,  -0.024870,  -0.075286,  -0.034492,  -0.049723,  -0.073601,   0.397645,   0.025388,  -0.031135,   0.109938,  -0.033271,  -0.219441,   0.228536,  -0.812980,   0.037701,
         0.002645,  -0.090422,  -0.140939,   0.037080,  -0.112221,   0.102773,  -0.133190,   0.021380,   0.298358,   0.165203,  -0.144975,  -0.357309,   0.063513,   0.127819,   0.188946,  -0.172210,   0.090238,   0.048266,  -0.136299,  -0.104118,  -0.120972,  -0.058413,  -0.093034,  -0.191443,  -0.106771,  -0.108013,   0.135796,   0.114834,  -0.307950,  -0.048838,  -0.138171,  -0.159446,
        -0.054628,  -0.068307,  -0.050362,  -0.026374,  -0.029907,  -0.009259,   0.020514,   0.075645,   0.271392,   0.006871,   0.274956,   1.191390,   1.050600,   0.101674,   0.170077,   0.146570,  -1.067030,   0.272101,  -0.380558,  -0.106651,   0.152506,  -0.250475,  -0.248878,  -0.480365,   0.009851,   0.141981,  -0.157515,  -0.088100,  -0.275980,   1.048420,   0.745840,   0.683178,
        -0.144970,   0.022936,   0.751199,   0.263960,  -0.011449,   0.156764,  -0.071619,  -0.001325,   0.485991,   0.194255,  -0.202096,   0.218626,  -0.020363,   0.109205,  -1.388530,   0.052950,  -0.123499,   0.206860,   0.434674,  -0.279152,  -0.064567,   0.049823,   0.176504,   0.159918,   0.184905,   0.189685,   0.158141,  -0.018026,   0.099066,   0.115148,  -0.008789,   0.058037,
        -0.053615,   0.090443,   0.446703,  -0.175185,   0.962793,  -0.359832,   0.260247,   0.408442,   0.889547,  -0.006022,  -0.241770,   0.583929,   0.168051,   0.099749,  -0.142833,  -0.023761,   0.621997,  -0.241698,  -0.161116,  -0.356706,   0.098813,  -0.112026,  -0.045344,  -0.077618,  -0.110581,   0.198996,  -0.058546,  -0.054945,   0.106576,  -0.009781,  -0.051788,   0.099840,
         0.671313,  -0.101275,   0.478221,  -0.027098,   0.034182,   0.203631,   0.487882,   0.041909,  -0.155440,   0.156926,   0.174254,  -0.576482,   0.734901,   0.283095,   0.412607,  -0.119958,   0.010733,   0.104436,  -0.118115,   0.270383,   0.720090,  -0.043421,   0.603588,  -0.163204,   0.059322,  -0.126127,  -0.052851,   0.145525,  -0.014219,   0.014972,   0.255258,   0.168249,
        -0.153093,  -0.232742,  -0.001186,   0.398301,  -0.039815,   0.176651,   0.031367,  -0.271366,  -0.054525,  -0.120865,  -0.186195,   0.086472,  -0.043092,   0.049131,   0.076345,  -0.546845,   0.093137,   0.055682,  -0.557429,   0.084064,  -0.155668,  -0.164382,  -0.088207,   0.064833,  -0.171008,  -0.058981,  -0.014762,  -0.077713,  -0.107826,   0.089662,   0.444223,  -0.124921,
         0.224062,   0.276110,   0.827445,  -0.070426,  -0.175760,   0.232169,   0.028783,  -0.451292,  -0.683193,  -0.288323,   0.246580,  -0.045074,  -0.293156,  -0.501329,  -0.040134,   0.062840,  -0.389639,  -0.202258,  -0.336944,   1.049810,   0.495628,   0.539256,  -0.149125,   0.126186,   0.450352,   0.172697,  -0.010709,   0.206156,   0.181147,   0.498733,   0.333884,  -1.038850,
         0.013610,   0.029768,  -0.130925,  -0.034161,  -1.692710,   0.092681,   0.005458,  -0.026663,   0.214854,  -0.167585,  -0.152733,  -0.433206,   0.248055,   0.242517,   0.123144,   0.254107,  -0.143018,  -0.070978,  -0.068711,   0.054323,  -0.003409,   0.160834,  -0.098650,   0.005358,   0.216943,  -0.045060,   0.592393,  -0.143087,   1.234870,   0.903399,   0.395649,  -0.370653,
         0.065621,   0.015549,   0.172930,   0.036922,  -0.434578,  -0.336130,   0.331833,  -0.435118,  -0.372436,   0.181565,   0.264616,  -0.011786,  -0.103437,  -0.304516,  -0.336751,   0.050074,  -0.116334,  -0.121250,   0.029105,   0.224980,  -0.110255,   0.023576,   0.168733,   0.059167,   0.545938,   0.196964,   0.356153,  -0.122727,   0.044041,  -0.093723,  -0.076734,  -0.167500,
        -0.089644,  -0.332469,   0.709446,   0.711087,   0.208682,  -0.021742,  -0.037474,   0.077366,  -0.042180,   0.176749,  -0.228190,   0.106728,   1.107570,   0.211422,   0.142349,  -0.183774,
         0.187336,  -0.115158,  -0.048793,  -0.030604,   0.053631,   0.265484,  -0.039272,  -0.069014,   0.683024,   0.109637,   0.257113,  -0.005864,   0.337621,  -0.043847,   0.143826,   0.024831,   0.096220,   0.191665,   0.248142,   0.398024,  -0.108925,  -0.143601,   0.230377,  -0.021784,  -0.170232,  -0.165712,   0.013116,   0.000608,  -0.320088,   0.006401,  -0.112846,   0.074877,
         0.131336,  -0.001154,   0.589270,   0.171814,  -0.024573,  -0.031922,   0.167363,  -0.067878,   0.874893,   0.065639,  -0.378159,  -0.302894,  -0.032058,   0.070697,  -0.310110,  -0.404861,   0.218217,  -0.197258,   0.026566,   0.249255,  -0.341583,   0.215241,   0.108545,  -0.442558,  -0.298922,  -0.318520,  -0.198962,  -0.025490,   0.067081,  -0.079502,   0.130369,   0.304380,
        -0.249620,   0.178490,   0.129515,   0.004519,   0.102010,   0.083665,   0.674687,  -0.150004,   0.366915,  -0.167555,  -0.304723,   0.120575,   0.078282,   0.326395,  -0.354519,   0.116820,  -0.249941,  -0.268199,   0.372439,   0.029786,  -0.096047,  -0.142767,   0.152632,  -0.158579,   0.356975,  -0.012400,   0.023047,   0.260909,   0.070807,  -0.020344,   0.649709,   0.087929,
         0.049846,   0.254009,   0.241175,  -0.026813,   1.013000,  -0.040336,  -0.038973,  -0.225459,   0.216090,   0.304313,  -0.427733,  -0.188830,  -0.071091,  -0.022513,   0.331890,   0.065966,  -0.247292,  -0.060480,   0.024765,  -0.160678,  -0.180214,  -0.131144,  -0.128128,   0.019400,   0.131696,  -0.007955,   0.008974,   0.086677,  -0.039681,   0.145469,  -0.089375,   0.110647,
         0.337797,   0.050233,   0.488045,  -0.111406,   0.229680,  -0.232981,   0.273662,  -0.063671,  -0.084327,   0.155273,  -0.203473,   0.127221,  -0.104751,  -0.212917,   0.123216,  -0.208282,  -0.064768,   0.074837,   0.033594,   0.039205,  -0.242153,  -0.061467,  -0.151763,   0.035783,  -0.167039,   0.034632,   0.029267,   0.079215,  -0.039827,   0.012327,   0.287278,   0.187651,
        -0.096626,   0.029758,   0.303828,  -0.008133,   0.004116,  -0.106899,  -0.066822,   0.093864,   0.433585,   0.108969,   0.005820,   0.040508,   0.271067,  -0.023816,  -0.306561,   0.038606,  -0.101330,  -0.079244,  -0.070341,  -0.363654,   0.002359,  -0.139646,  -0.003845,   0.173961,   0.211928,   0.168837,   0.074634,   0.045852,   0.134997,   0.129938,   0.829621,  -0.062451,
         0.044471,  -0.223132,   0.040448,  -0.105976,  -0.197255,  -0.110997,   0.210766,   0.194079,   0.004233,   0.334455,  -0.083418,  -0.002056,   0.062452,  -0.153031,  -0.389301,  -0.183028,  -0.315369,  -0.060389,   0.069623,   0.139162,   0.166803,   0.189794,  -0.126093,   0.213416,   0.126234,   0.208081,   0.033632,   0.208717,   0.490973,  -0.070010,   0.334718,  -0.232940,
        -0.146932,  -0.007887,   0.048442,   0.280502,  -0.400862,   0.091113,  -0.106179,  -0.393543,   0.332673,  -0.109563,  -0.260408,   0.032393,   0.111077,  -0.136129,  -0.107649,   0.021985,  -0.135112,   0.231349,   0.164247,   0.022658,   0.448291,   0.007840,   0.028467,   0.150905,   0.122708,   0.241463,   0.551041,  -0.054240,  -0.047787,  -0.134463,   0.004530,   0.027325,
        -0.229947,  -0.258533,  -0.011131,  -0.099814,   0.764177,   0.160939,  -0.245250,   0.222494,   0.023401,  -0.279079,  -0.219519,  -0.203044,  -0.180337,  -0.179719,   0.354181,  -0.055507,  -0.074340,   0.066509,  -0.148121,   0.136777,  -0.104495,  -0.002869,   0.075612,   0.045039,   0.875023,  -0.073083,   0.100453,  -0.125666,   0.409902,  -0.134268,   0.190034,   0.394684,
         0.220680,   0.226025,  -0.296871,  -0.304071,   0.059149,  -0.088414,  -0.145047,   0.201010,   0.193403,   0.186124,  -0.091936,  -0.070988,  -0.049893,   0.462674,   1.028460,   1.517760,
        -0.212001,  -0.015757,   0.075745,  -0.091334,   0.034895,  -0.180947,   0.029667,   0.045235,  -0.491004,  -0.094268,  -0.131704,   0.121571,   0.016137,  -0.091032,  -0.236343,   0.289825,  -0.085293,  -0.149716,   0.192791,  -0.164685,   0.276731,  -0.221966,  -0.249765,   0.119622,  -0.023491,  -0.362657,   0.019846,  -0.036502,  -0.028216,  -0.040739,   0.076021,   0.208777,
         0.263824,   0.131942,  -0.012516,  -0.125849,   0.240231,   0.103034,  -0.095203,   0.069678,  -0.337281,  -0.042811,  -0.177706,  -0.210512,  -0.522187,   0.394102,  -0.201413,   0.060835,  -0.474224,  -0.234143,   0.537527,   0.326251,  -0.088687,   0.172312,  -0.024623,   0.161751,   0.091749,  -0.002686,   0.007516,  -0.199442,   0.232134,  -0.165362,   0.349237,  -0.147586,
         0.273859,   0.205678,  -0.037312,  -0.235817,  -0.182540,  -0.155519,   0.093129,   0.006111,   0.039677,   0.380756,   0.300899,  -0.036637,   0.944725,   0.531114,   0.086634,   0.581311,  -0.121195,   0.046735,   0.108378,   0.012335,  -0.085267,   0.345239,   0.144300,   0.341473,   0.005656,   0.109728,  -0.080807,   0.033738,  -0.091009,  -0.185829,  -0.087435,  -0.202237,
         0.025043,  -0.107440,  -0.250141,  -0.033128,  -0.251717,   0.043453,  -0.114424,  -0.240003,  -0.033692,   0.092935,  -0.106087,  -0.069651,   0.092519,  -0.119016,   0.085837,  -0.119520,  -0.228126,  -0.036714,  -0.053242,  -0.128037,   0.052026,   0.067304,   0.087717,   0.093422,  -0.123162,   0.023370,  -0.103147,   0.038126,  -0.016184,  -0.019100,  -0.058586,  -0.008912,
        -0.641791,  -0.526757,   0.477250,   0.040754,   0.436338,   0.458604,   0.169401,  -0.006047,  -0.036718,   0.424229,   0.248178,   0.169843,  -0.177954,  -0.298416,  -0.086201,   0.038849,  -0.139487,   0.609584,   0.151913,   0.019592,  -0.150440,  -0.009248,  -0.276807,   0.057374,   0.127893,   0.031545,   0.183051,   0.028269,   0.083000,  -0.031754,  -0.740563,  -0.886280,
         0.806734,   0.127867,   0.226935,  -0.116314,   0.017852,   0.012462,   0.330161,   0.113798,   0.333022,  -0.089231,   0.201732,  -0.085822,   0.174676,   0.145695,   0.066183,  -0.275160,   0.081504,  -0.100861,   0.214871,   0.397886,   0.047916,   0.145971,   0.111887,   0.090670,  -0.073402,   0.084170,   0.298264,  -0.069826,   0.013940,  -0.329209,   0.004936,   0.072435,
         0.515802,   0.295059,  -0.163426,  -0.077768,   0.480837,   0.670718,  -0.035126,   0.433879,   0.151078,   0.093456,   0.548466,  -0.047492,   0.125329,   0.077970,   0.376835,   0.250432,   0.031477,   0.190429,  -0.105578,  -0.362561,  -0.092043,  -0.403394,   0.625163,  -0.031970,  -0.126369,  -0.045172,   0.048412,  -0.097526,  -0.138587,   0.129325,  -0.224638,  -0.300174,
        -0.863709,   0.229701,   0.143291,  -0.166253,  -0.168840,  -0.252749,   0.100816,   0.127909,  -0.401910,   0.104397,   0.093350,  -0.279432,  -0.312158,  -0.545716,  -0.233418,   0.095177,   0.281093,   0.096189,  -0.035717,   0.078108,  -0.091223,  -0.007483,   0.200729,  -0.055542,  -0.712000,   0.172560,  -0.634852,   0.105268,   0.241409,   0.176854,  -0.086747,   0.071074,
         0.435095,   0.329365,  -0.361762,  -0.114359,   0.063170,   0.036068,   0.255401,  -0.107475,  -0.112329,   0.107369,   0.193295,   0.254431,   0.099331,   0.245177,   0.130847,   0.144004,   0.072642,   0.021977,   0.365142,  -0.070886,   0.256877,  -0.030437,  -0.762127,   0.348821,  -0.363823,   0.056494,  -0.330920,  -0.248454,  -0.419823,   0.485111,  -0.162140,  -0.362091,
        -0.337808,  -0.266804,   0.713315,   0.160423,  -0.477172,   0.063070,   0.323249,  -0.084389,  -0.327118,  -0.192134,  -0.095482,  -0.105902,   0.009599,   0.540114,   0.712648,   0.551965,
        -0.161920,  -0.071053,  -0.098498,  -0.103379,   0.054634,   0.021274,  -0.003120,   0.001141,  -0.402708,  -0.251886,  -0.195645,  -0.020001,   0.434974,  -0.024036,   0.221338,   0.089511,   0.039517,   0.150197,  -0.112407,  -0.030996,   0.174297,  -0.229176,  -0.035184,   0.168195,  -0.066937,  -0.029084,   0.178465,  -0.082949,  -0.645047,   0.086803,   0.012415,   0.000507,
        -0.021982,   0.040660,  -0.147757,  -0.154914,   0.072268,  -0.091757,  -0.135824,  -0.150129,   0.110588,  -0.036616,   0.124939,   0.208770,   0.289914,   0.264029,  -0.104095,   0.061595,   0.859314,  -0.054520,  -0.403503,  -0.081494,  -0.183808,  -0.028668,  -0.035825,   0.034780,   0.338137,   0.071417,  -0.235649,   0.132832,  -0.189127,   0.013277,   0.573519,   0.194369,
         0.057017,   0.192454,  -0.051611,  -0.014151,  -0.090163,  -0.071541,   0.241839,   0.028787,   0.609517,   0.685590,   0.497278,   0.341951,   0.635893,   0.452865,  -0.045617,   0.521515,  -0.504896,  -0.293166,   0.349432,  -0.014916,  -0.204082,   1.006130,   0.445987,   0.421960,   0.028948,   0.262593,  -0.386897,   0.102445,   0.219848,   0.057128,  -0.117717,  -0.027607,
         0.057861,   0.078358,   0.516224,   0.019343,  -0.306668,  -0.027845,  -0.068468,  -0.107200,   0.376504,   0.203117,  -0.155664,  -0.032742,   0.348239,  -0.026283,  -0.154841,  -0.310345,   0.233553,   0.086601,  -0.162162,   0.385203,   0.131301,  -0.009066,  -0.355836,   0.270744,  -0.051850,   0.189312,   0.070699,   0.152338,  -0.144171,  -0.038821,   0.040155,   0.096068,
        -0.550514,  -0.219585,   0.026159,   0.107590,   0.480059,   0.075420,   0.027318,   0.005724,  -0.093398,   0.031029,   0.047709,  -0.234563,   0.014043,  -0.254700,  -0.027175,  -0.120918,  -0.291668,   0.111918,  -0.061535,  -0.221500,   0.362600,   0.434272,   0.052874,   0.149675,   0.290445,  -0.055840,   0.183608,   0.286064,   0.210569,   0.162015,  -1.663720,  -0.820610,
         0.980953,  -0.000198,   0.249531,  -0.415556,  -0.453123,  -0.051280,   0.354128,  -0.180723,   0.371876,  -0.299789,   0.321644,  -0.153188,   0.220397,   0.193265,   0.289419,  -1.175500,   0.085804,   0.052793,   0.286607,   0.074094,   0.177660,   0.097084,   0.243682,   0.212273,  -0.078932,   0.197576,   0.813025,  -0.008112,  -0.436483,   0.606566,  -0.369577,   0.296286,
         0.214814,  -0.301206,  -1.020800,   0.224772,   0.436384,   0.023738,  -0.374062,   0.030968,   1.055600,   0.378124,   0.065054,   0.269013,   0.607100,  -0.503605,   0.138403,   0.098451,   0.026875,  -0.036893,   0.402286,  -0.406150,   0.698059,  -0.173275,   0.584284,   0.383969,   0.077727,  -0.139699,  -0.868899,   1.134030,  -0.614023,   0.397833,   0.012488,  -0.885399,
        -1.827220,   0.388517,   0.332841,  -0.139920,  -0.440321,  -0.126179,   0.558900,   0.177365,  -0.488757,   0.086573,   0.210543,  -0.856252,  -0.517838,  -0.402729,   0.142688,  -0.097343,   1.165770,   0.244439,   0.394742,   0.206664,   0.116850,   0.074045,   0.468200,   0.092037,  -1.491660,   1.815160,  -0.959031,   0.060351,  -0.252646,  -0.268586,  -0.307586,   0.156375,
        -0.024576,  -0.225803,  -0.706723,  -0.333994,   0.185299,   0.039387,  -0.406024,  -0.007569,  -0.081577,  -0.825654,   0.198440,   0.156639,   0.268301,   0.296129,   0.365206,   0.022315,   0.393524,   0.110739,   0.249201,   0.358335,   0.413635,   0.109490,  -1.412350,   1.323160,  -0.302389,  -0.084755,  -0.183293,  -0.244071,  -0.827773,   0.864307,  -0.183781,  -0.498565,
        -0.351390,  -0.189676,   0.055920,  -0.118418,  -0.485331,  -0.244343,   0.011681,   0.478297,  -0.097880,  -0.300140,  -0.254393,  -0.102893,  -0.344171,   0.246005,  -0.890107,   1.070080,
        -0.265015,  -0.044871,   0.155691,  -0.017704,  -0.024167,  -0.076647,   0.065880,   0.061629,  -0.273517,  -0.742900,  -0.060682,  -0.006792,   0.241250,   0.017164,  -0.099282,   0.162655,  -0.123253,   0.200435,   0.014384,   0.057732,   0.066363,   0.005799,   0.030678,   0.056475,   0.012370,  -0.089529,  -0.116731,   0.004443,   0.031989,   0.208021,   0.090111,   0.129131,
         0.288568,   0.316073,  -0.090639,  -0.149225,   0.164173,  -0.007673,   0.152797,  -0.283892,  -0.441956,  -0.183153,   0.076288,  -0.135296,  -0.310059,  -0.056970,  -0.053166,   0.378449,  -0.575801,   0.038234,   0.256055,   0.151452,   0.100297,   0.046255,   0.073923,   0.123508,   0.086461,   0.234079,   0.222219,   0.344276,  -0.009185,  -0.303894,   0.100218,  -0.225950,
         0.356615,   0.164815,  -0.073993,  -0.189764,   0.022630,  -0.252712,  -0.126575,  -0.024891,  -0.033201,  -0.031523,  -0.098430,  -0.094507,   0.162544,   0.096582,  -0.054668,   0.097678,   0.021998,   0.195496,  -0.176576,   0.132108,  -0.070061,  -0.055060,  -0.270505,  -0.075858,   0.055334,  -0.039774,  -0.142839,  -0.003973,  -0.135091,  -0.042952,  -0.309275,  -0.155853,
         0.077816,   0.041266,   0.142008,  -0.077315,  -0.572707,   0.085650,   0.088858,  -0.005829,  -0.032703,   0.046355,   0.149674,   0.115064,  -0.416118,   0.086404,   0.077563,  -0.078902,  -0.020454,  -0.040180,   0.018158,  -0.142069,   0.204780,   0.026293,   0.024813,   0.114668,   0.178856,   0.080331,  -0.041381,  -0.058374,  -0.008711,  -0.115247,   0.123104,   0.010206,
        -0.567864,  -0.400919,   0.112648,   0.020766,  -0.048771,  -0.098034,  -0.341729,   0.127373,  -0.083594,   0.059473,  -0.150449,   0.068009,  -0.006352,  -0.020702,  -0.274191,   0.118802,  -0.051968,   0.324299,  -0.070959,  -0.228137,  -0.141550,  -0.007199,   0.057018,   0.060693,   0.165242,  -0.021015,  -0.111708,   0.087738,   0.211588,  -0.060811,  -0.502401,  -0.794530,
         0.346992,   0.105777,   0.245011,  -0.114735,  -0.249794,   0.102770,   0.291909,   0.065738,   0.034755,  -0.062878,   0.247472,  -0.060710,   0.112562,   0.028395,   0.034198,  -0.303646,  -0.073105,   0.061565,   0.049875,   0.081490,   0.297954,   0.223257,   0.063045,   0.033981,   0.040673,   0.076149,   0.249559,   0.072930,   0.363174,  -0.240833,  -0.303163,   0.055151,
         0.001872,  -0.089935,  -0.491532,   0.105312,   0.234298,   0.247395,  -0.323600,  -0.030739,   0.228968,   0.177791,   0.014137,   0.045063,   0.157244,   0.328591,   0.025932,   0.061138,   0.216919,   0.325943,   0.168900,  -0.354842,  -0.117927,  -0.320245,   0.240109,   0.034985,  -0.281250,  -0.066082,   0.585972,  -0.407867,   0.129322,  -0.051179,   0.344480,  -0.063212,
        -0.571161,   0.001485,   0.611312,   0.196270,  -0.322456,   0.071288,   0.041716,   0.019405,  -0.170932,   0.081445,  -0.048334,   0.053562,  -0.355067,  -0.137165,  -0.113943,   0.035091,   0.169966,   0.071954,   0.156880,  -0.003585,  -0.082862,  -0.132147,   0.192392,   0.014085,  -0.291060,  -0.082621,  -0.405065,   0.237067,   0.036103,  -0.041392,  -0.047708,   0.188226,
         0.204178,   0.183531,  -0.276441,   0.017749,  -0.012262,   0.091366,   0.071783,   0.021766,  -0.153906,  -0.142923,   0.181073,   0.048790,   0.052981,   0.100546,   0.077999,   0.012063,   0.003389,  -0.113243,   0.056564,  -0.198775,  -0.059712,  -0.070356,  -0.108008,  -0.070947,  -0.271319,   0.003011,  -0.222162,  -0.003582,  -0.187385,   0.038691,   0.055041,  -0.019501,
        -0.163911,   0.080605,   0.077964,  -0.142269,  -0.221642,  -0.155916,  -0.103204,   0.285404,   0.035849,  -0.016503,  -0.044956,  -0.098682,  -0.139764,   0.402047,  -0.814325,   0.591613,
         0.176325,  -0.042399,  -0.104134,  -0.096545,   0.285625,   0.058725,  -0.018270,   0.098494,   0.972569,   0.457919,  -0.147989,   0.169656,  -0.146725,  -0.139284,  -0.059995,  -0.086231,  -0.135176,  -0.122646,   0.392949,  -0.060886,  -0.059328,  -0.059641,  -0.347721,   0.053238,  -0.190742,  -0.132982,   0.000141,  -0.058839,  -0.148364,  -0.032937,   0.212743,  -0.066931,
        -0.050063,   0.037376,   0.319144,   0.209579,   0.056200,  -0.026804,   0.397186,   0.117281,   0.589920,   0.119386,  -0.388750,  -0.236956,   0.111607,  -0.152613,  -0.015353,  -0.357803,   0.553024,  -0.000159,  -0.051831,   0.070138,  -0.466720,   0.156912,   0.374506,  -0.796158,  -0.143045,  -0.244277,  -0.225442,  -0.153667,  -0.002701,  -0.124385,  -0.094729,   0.217918,
        -0.172749,   0.206481,  -0.036656,  -0.067218,   0.210420,   0.168927,   0.281439,   0.018402,   0.418129,   0.397607,   0.889577,  -0.055726,   0.176687,   0.519441,   0.364191,   0.132756,  -0.550843,  -0.571144,   0.343017,  -0.028562,  -0.162385,   0.215008,   0.439797,   0.047683,   0.097884,  -0.196691,  -0.191198,  -0.055112,  -0.090693,   0.044044,   0.821421,   0.008216,
         0.089814,   0.054937,   0.298361,  -0.021810,   0.609491,   0.218680,  -0.279673,  -0.245309,   0.181306,  -0.039822,  -0.175925,  -0.145715,   0.160454,  -0.019354,   0.212403,   0.092070,  -0.300036,   0.135563,   0.055373,  -0.309688,  -0.042689,  -0.259250,  -0.366334,  -0.169433,   0.032416,  -0.149048,   0.014353,   0.157433,   0.048929,   0.191501,  -0.119844,   0.045253,
         0.156199,   0.511173,   0.313138,   0.018114,   0.230843,  -0.144858,   0.126069,  -0.176823,   0.021985,   0.119470,   0.154171,   0.324675,  -0.224327,  -0.240804,   0.088721,  -0.198969,  -0.005058,  -0.189682,   0.050213,   0.030202,  -0.049466,   0.047184,  -0.208386,   0.031296,  -0.154086,  -0.019138,   0.062260,   0.070215,  -0.015579,   0.058971,   1.234650,   0.210641,
        -0.221062,   0.209510,   0.160404,  -0.139419,  -0.028367,  -0.131930,  -0.178881,  -0.025264,   0.243186,   0.207030,   0.108091,   0.114013,  -0.055545,   0.012636,  -0.074042,  -0.038075,   0.069475,  -0.147112,  -0.131156,  -0.200342,  -0.068106,   0.108214,   0.100649,   0.032153,   0.522780,   0.149238,   0.117886,   0.111783,   0.560963,   0.238594,   0.264711,  -0.019644,
        -0.189187,  -0.280658,  -0.044662,   0.101446,  -0.308057,  -0.445999,   0.294074,   0.050568,   0.154950,   0.370777,  -0.442612,   0.147151,   0.004642,  -0.799886,  -0.264995,  -0.136238,  -0.092679,  -0.063631,   0.240554,  -0.019364,  -0.096163,   0.062876,  -0.300565,   0.337372,  -0.004896,  -0.092674,   0.239718,   0.317722,   0.176730,  -0.121987,   0.415160,   0.248551,
         0.527852,   0.041541,   0.341369,   0.778873,   1.287650,   0.474757,  -0.459714,  -0.644224,   0.322306,  -0.147214,  -0.266770,   0.129469,   0.305106,   0.092575,  -0.052358,  -0.109867,  -0.149232,   0.065395,   0.012980,   0.193731,   0.180148,   0.192459,   0.141976,   0.009334,   0.483033,   0.295520,   0.358794,   0.098955,  -0.220080,  -0.291971,   0.001662,  -0.170934,
        -0.179749,  -0.386159,   0.045318,  -0.113451,   0.430614,   0.077706,  -0.406066,   0.150966,  -0.012113,  -0.545318,  -0.122996,  -0.018857,  -0.112254,  -0.084816,   0.100516,   0.050009,  -0.060133,   0.044188,   0.052443,   0.120001,   0.177094,   0.073632,   0.270926,   0.434867,   0.299083,  -0.065775,  -0.089366,  -0.167868,  -0.030602,  -0.196345,   0.088754,   0.073220,
        -0.055321,   0.106731,  -0.145107,  -0.209440,  -0.071776,  -0.063058,  -0.070205,  -0.331587,   0.318953,   0.106032,  -0.073974,  -0.002130,  -0.065293,   0.269745,   2.416530,   1.622040,
         0.256072,  -0.086739,   0.018478,  -0.001661,  -0.097009,   0.224562,  -0.068507,  -0.012125,   0.247312,   0.159926,   0.121178,  -0.097314,   0.193739,  -0.099047,   0.107577,   0.034855,   0.055720,   0.206004,   0.160616,   0.148608,  -0.068940,  -0.055786,   0.132677,   0.038477,  -0.250206,  -0.293005,  -0.014761,  -0.002091,  -0.278406,   0.048871,  -0.074912,   0.037217,
         0.111450,   0.038324,   0.459240,   0.149224,  -0.011441,   0.100466,   0.108404,   0.011409,   0.818354,   0.162515,  -0.154494,  -0.340065,   0.042118,  -0.105880,  -0.100238,  -0.394904,   0.233933,  -0.096090,  -0.003153,   0.346615,  -0.146924,   0.143928,   0.017440,  -0.413597,  -0.207728,  -0.387262,  -0.198551,  -0.211035,  -0.008056,  -0.055225,   0.188577,   0.523969,
        -0.071543,   0.178324,   0.157952,  -0.010963,   0.163921,   0.040446,   0.535088,  -0.180604,   0.389370,  -0.106824,  -0.452557,  -0.192833,   0.141492,   0.380821,  -0.361395,   0.455182,  -0.225237,  -0.277708,   0.406748,  -0.051619,  -0.001786,  -0.021087,   0.048828,  -0.026910,   0.965789,   0.361526,  -0.006696,   0.159933,   0.069959,   0.031210,   0.537216,  -0.013762,
         0.005247,   0.108964,   0.400293,  -0.003498,   0.809441,   0.083045,   0.020899,  -0.099748,   0.143382,   0.384479,  -0.274026,  -0.126867,  -0.386578,  -0.065872,   0.174926,   0.119250,  -0.061149,   0.125904,  -0.249630,  -0.185811,  -0.156872,  -0.108045,  -0.238383,  -0.032456,   0.052545,  -0.000651,   0.202313,   0.087550,  -0.040457,   0.190561,  -0.081840,  -0.081882,
         0.312005,  -0.032006,   0.539081,  -0.091865,   0.340269,  -0.173571,   0.288784,  -0.197457,  -0.066230,   0.131279,  -0.215263,   0.303669,  -0.144170,  -0.259941,   0.309430,  -0.117402,  -0.014424,   0.274291,  -0.205608,  -0.193459,   0.079894,   0.029385,  -0.164414,   0.117533,  -0.047188,  -0.017963,   0.090173,   0.275205,   0.045970,   0.018908,   0.301842,   0.132088,
        -0.043011,   0.098253,   0.215581,  -0.083977,   0.018583,   0.147548,  -0.310878,   0.025402,   0.243640,   0.107276,   0.236822,   0.074996,   0.336231,   0.128785,  -0.175538,  -0.049503,  -0.077857,  -0.115670,  -0.190851,  -0.503484,   0.103854,   0.003565,   0.161271,   0.087710,   0.149660,   0.180376,   0.032576,  -0.024304,   0.087471,   0.142299,   0.689560,  -0.028933,
        -0.004936,  -0.128017,   0.085796,  -0.282773,  -0.175199,  -0.128442,   0.082706,  -0.000115,   0.083949,   0.593646,  -0.062010,   0.070657,  -0.049327,  -0.075375,  -0.502083,  -0.230759,  -0.257271,  -0.075871,   0.323220,   0.231262,   0.045749,   0.141518,  -0.191975,   0.310953,  -0.001717,   0.096506,   0.118423,   0.159455,   0.378078,   0.010921,   0.425186,   0.073453,
        -0.355846,  -0.178308,  -0.051572,   0.443641,  -0.129249,   0.185225,  -0.090761,  -0.250460,   0.333132,  -0.189066,  -0.179079,   0.064756,   0.174822,  -0.033752,  -0.220183,  -0.098486,  -0.008841,   0.278092,   0.195079,   0.198073,   0.357715,   0.151262,  -0.006517,   0.205507,  -0.005477,   0.120076,   0.374162,  -0.128907,  -0.074439,  -0.166421,  -0.050383,  -0.294905,
        -0.189206,  -0.084398,  -0.284561,   0.152635,   0.730941,   0.037603,   0.017868,   0.137125,   0.083815,  -0.282092,  -0.149530,  -0.076570,   0.018936,   0.035890,   0.349872,  -0.069601,  -0.121664,  -0.072686,  -0.096623,   0.157428,   0.049113,   0.160694,   0.233679,   0.035111,   1.105510,   0.000368,   0.122255,  -0.272457,   0.190101,  -0.099789,   0.165533,   0.418179,
         0.058498,   0.042972,  -0.134210,  -0.162407,   0.111963,   0.074848,  -0.385291,   0.656795,   0.174468,   0.163095,  -0.126335,  -0.155469,   0.022893,   0.285146,   1.000140,   1.669580,
         0.102710,  -0.173897,  -0.069438,  -0.010982,  -0.138856,   0.262096,  -0.050098,   0.013484,   0.804394,   0.053503,  -0.070603,  -0.118865,   0.319169,   0.120969,   0.197543,   0.031125,   0.131292,   0.247792,   0.168199,  -0.044465,  -0.118737,  -0.118483,   0.186541,   0.035878,  -0.185392,   0.034844,   0.238442,   0.023545,  -0.298064,   0.066968,  -0.167135,  -0.023211,
         0.104404,  -0.004535,   0.043292,   0.163478,  -0.111877,   0.028761,   0.220769,   0.027398,   0.644480,   0.136653,  -0.092773,  -0.163308,   0.310429,  -0.174744,  -0.028732,  -0.206090,   0.335543,  -0.041585,  -0.260317,   0.075035,   0.232357,   0.020138,   0.002271,  -0.290678,  -0.061712,  -0.021304,   0.381532,  -0.018983,  -0.175012,   0.205419,   0.911951,   0.825705,
        -0.141529,   0.117559,   0.209066,   0.249510,   0.010114,   0.006192,   0.478113,   0.053412,   1.013440,   0.188981,  -0.457725,  -0.107210,   0.227648,   0.280333,  -1.004660,   0.545659,  -0.209260,  -0.222132,   0.349248,  -0.075517,  -0.108360,  -0.697690,   0.284754,  -0.048604,   0.874916,   0.820893,  -0.025182,   0.109957,   0.178889,   0.059910,   0.143672,   0.127909,
        -0.073520,  -0.030621,   0.152085,   0.023463,   0.923764,  -0.127003,   0.240120,   0.057993,   0.420884,  -0.198022,  -0.087866,  -0.003670,   0.203207,   0.144721,  -0.168286,  -0.203095,   0.239992,  -0.004404,  -0.153073,  -0.051599,  -0.058400,  -0.095880,   0.311353,   0.107794,  -0.041018,  -0.053412,  -0.030733,   0.134213,  -0.060645,   0.085775,  -0.043423,   0.002420,
         0.414855,  -0.112071,   0.925970,  -0.085297,   0.492317,  -0.009683,   0.071326,  -0.272721,  -0.147372,   0.120667,  -0.010386,   0.018062,  -0.039918,  -0.225091,   0.118917,  -0.139726,  -0.102833,   0.337646,  -0.080104,   0.078663,   0.557789,   0.095801,  -0.170839,  -0.130018,  -0.070269,   0.082679,   0.026248,   0.274999,  -0.106178,  -0.062604,   1.125780,   0.066008,
         0.604424,  -0.102841,   0.324788,   0.083359,  -0.309221,   0.028438,  -0.093248,  -0.060059,   0.514940,  -0.025226,  -0.042122,   0.012111,   0.362441,   0.061764,  -0.040723,  -0.072489,  -0.125947,  -0.187157,  -0.449278,  -0.107116,  -0.188028,  -0.079639,   0.025105,  -0.011177,  -0.082034,   0.144201,   0.063248,  -0.121912,   0.189300,  -0.076333,   1.593830,  -0.090379,
         0.047316,  -0.046714,  -0.120132,  -0.008607,   0.034286,  -0.165325,  -0.748608,  -0.331460,  -0.202648,   0.181828,   0.115988,   0.054401,  -0.136854,  -0.318117,  -0.167272,   0.041948,  -0.408269,  -0.179994,   0.123637,   0.493793,   0.335171,   0.568513,  -0.256122,   0.164214,   0.246846,   0.102658,   0.194635,   0.022864,   0.370460,   0.011064,   0.590405,  -0.761473,
         0.223971,  -0.097849,   0.338232,   0.473384,  -0.838704,   0.271480,  -0.333943,  -0.302232,   0.364891,   0.068765,  -0.211403,  -0.184010,   0.114061,   0.166028,   0.060583,   0.512024,  -0.119213,   0.076791,   0.084565,   0.080089,   0.143127,   0.142894,   0.046295,   0.020605,   0.335221,  -0.042485,   0.647806,  -0.109036,   0.216717,   0.128678,   0.409400,  -0.540019,
        -0.076743,  -0.183661,   0.019366,   0.076603,  -0.132429,  -0.133854,  -0.083561,  -0.184644,   0.153838,   0.119078,  -0.108007,  -0.028261,   0.321383,   0.128944,  -0.303423,   0.002971,  -0.132437,  -0.125599,  -0.031267,   0.009486,  -0.151219,  -0.015321,   0.597842,  -0.071850,   1.546820,   0.096223,   0.103656,  -0.140569,  -0.010001,  -0.334376,   0.122038,   0.009481,
         0.203953,  -0.064040,  -0.055738,  -0.106144,   0.077283,  -0.145215,  -0.327321,   0.170655,   0.135007,   0.107429,   0.050810,  -0.006196,  -0.031966,   0.170364,   1.872900,   0.594007,
        -0.127394,   0.099169,  -0.003649,  -0.160721,   0.039692,  -0.147112,   0.159413,   0.093973,  -0.429689,  -0.457825,  -0.009343,   0.069160,  -0.085459,  -0.137891,  -0.349965,   0.372415,  -0.102341,  -0.022574,   0.123694,  -0.210303,   0.330779,  -0.109487,  -0.180046,   0.131548,   0.081550,  -0.216221,  -0.088844,  -0.004046,  -0.137374,   0.082432,   0.011903,   0.091337,
         0.160840,   0.223754,  -0.128753,  -0.129003,   0.233900,  -0.027809,   0.006895,  -0.158924,  -0.364714,   0.023224,   0.082725,  -0.172792,  -0.425588,   0.094684,  -0.006742,   0.137196,  -0.554202,  -0.199160,   0.189946,   0.207837,   0.121923,  -0.004526,  -0.022943,   0.145077,   0.099158,   0.082729,   0.177778,   0.060189,   0.196294,  -0.390045,   0.205447,  -0.061788,
         0.269385,   0.228586,   0.008794,  -0.147951,   0.081057,  -0.205004,  -0.257902,  -0.037216,  -0.083729,   0.003109,  -0.226918,  -0.078326,   0.326426,   0.034820,   0.023887,   0.106187,   0.052679,   0.086210,  -0.113356,   0.183577,  -0.043399,  -0.110276,  -0.173728,   0.036189,   0.022316,   0.094523,  -0.115655,   0.171496,   0.021154,   0.004244,  -0.168250,  -0.116392,
         0.150515,  -0.070525,  -0.117103,   0.099426,  -0.754601,   0.117150,  -0.065829,  -0.034778,  -0.181990,   0.164135,  -0.026556,  -0.123057,  -0.307232,   0.067464,   0.099642,  -0.050185,  -0.023053,   0.016324,   0.015825,  -0.194915,   0.194122,   0.237858,  -0.017435,   0.071556,   0.070160,   0.107078,   0.019155,   0.081491,   0.041492,  -0.079502,   0.025829,  -0.008023,
        -0.538584,  -0.361784,  -0.096134,   0.040549,   0.035991,   0.161596,  -0.175328,   0.040177,  -0.052664,   0.065136,  -0.156518,   0.089896,   0.035852,  -0.004756,  -0.264772,   0.068279,  -0.003275,   0.481463,   0.114340,  -0.167267,   0.027629,   0.088849,  -0.140997,  -0.001050,   0.054568,  -0.123101,  -0.041654,  -0.008039,   0.045545,  -0.115099,  -0.674304,  -0.981779,
         0.552519,   0.157445,   0.273881,  -0.138994,   0.112195,  -0.035012,   0.403332,   0.146736,   0.537944,  -0.157922,   0.078671,  -0.190964,   0.219991,  -0.015641,   0.029079,  -0.337799,   0.007769,  -0.066973,  -0.053812,   0.058340,   0.219480,   0.152065,   0.210264,   0.229670,  -0.182935,   0.034935,   0.344778,   0.007514,   0.212348,  -0.224979,  -0.257341,   0.128493,
         0.248514,   0.247566,  -0.112290,  -0.023470,   0.565589,   0.672208,  -0.444735,   0.230668,   0.311288,  -0.051848,   0.465530,   0.034889,   0.151735,   0.180304,   0.230144,   0.275118,   0.042696,   0.188414,  -0.008061,  -0.279888,  -0.014268,  -0.371880,   0.277666,  -0.012650,  -0.176037,  -0.151171,   0.258507,  -0.125243,  -0.033130,  -0.018005,   0.000942,  -0.160276,
        -0.463703,   0.035767,   0.528668,  -0.051113,  -0.198479,   0.005077,   0.105622,   0.105006,  -0.214278,   0.059886,  -0.055482,  -0.102483,  -0.202307,  -0.299422,  -0.136561,   0.134384,   0.217549,   0.146414,   0.146481,  -0.060845,  -0.164369,  -0.088982,   0.154304,  -0.042164,  -0.473053,   0.145394,  -0.784314,   0.194926,   0.011067,   0.027440,   0.054470,   0.093153,
         0.147889,   0.112821,  -0.209953,  -0.001140,   0.075344,  -0.033482,   0.094316,   0.036890,  -0.090577,  -0.146707,   0.175324,   0.166918,   0.204862,   0.200658,   0.273080,   0.177054,   0.072379,  -0.131137,   0.264068,  -0.060362,   0.045321,   0.095256,  -0.475895,   0.120670,  -0.479653,   0.109239,  -0.231699,  -0.025597,  -0.189104,   0.341464,  -0.069297,  -0.199797,
        -0.224235,  -0.164981,   0.291178,  -0.139271,  -0.352500,  -0.095924,  -0.033788,   0.089836,  -0.225369,  -0.116745,  -0.137652,  -0.077492,  -0.063084,   0.707873,  -0.601748,   0.914449,
        -0.246579,   0.191297,   0.132012,  -0.149625,   0.053676,  -0.130894,   0.157640,  -0.000277,  -0.381241,  -0.790083,  -0.035682,  -0.037752,   0.162984,   0.106895,   0.151339,   0.000378,  -0.090414,   0.127310,   0.128348,   0.098683,   0.090300,  -0.159294,  -0.009442,   0.032747,   0.326588,  -0.056930,  -0.073289,   0.070142,   0.063020,   0.189491,  -0.037400,   0.043032,
         0.063059,  -0.008890,  -0.199437,  -0.108858,   0.172358,  -0.064654,  -0.031448,  -0.148529,  -0.495210,  -0.008685,   0.118765,   0.114419,  -0.164421,   0.047815,  -0.035802,   0.198059,  -0.262259,   0.081390,   0.139000,  -0.030128,   0.219174,  -0.139923,  -0.045890,   0.416717,   0.247598,   0.334066,   0.103847,   0.183244,   0.078330,  -0.257685,  -0.036672,  -0.358025,
         0.634194,  -0.003009,  -0.257054,  -0.140756,  -0.146151,  -0.354688,  -0.165090,   0.124889,  -0.650982,  -0.142084,  -0.470207,   0.061485,  -0.181383,  -0.571300,   1.821710,  -0.160991,   0.186229,   0.108434,  -0.579298,   0.063008,   0.203183,  -0.339348,  -0.354207,  -0.027903,  -0.004363,   0.014834,   0.153384,  -0.021055,  -0.089417,  -0.102041,  -0.299074,  -0.254792,
         0.138839,   0.016611,  -0.205357,   0.015286,  -0.613774,  -0.024338,  -0.059848,   0.070668,  -0.203307,   0.012903,   0.136771,  -0.072395,  -0.162100,   0.016703,   0.262552,  -0.052697,  -0.061694,   0.078314,  -0.005385,  -0.140883,   0.147318,   0.009266,   0.430548,   0.199572,   0.131566,   0.026295,   0.151738,  -0.034176,   0.061595,  -0.060873,  -0.002757,  -0.104965,
        -0.521339,  -0.124080,  -0.232783,   0.027242,  -0.267444,   0.219366,  -0.366715,   0.194220,  -0.087862,  -0.127151,   0.112193,  -0.099570,   0.193363,   0.200799,  -0.237552,   0.229937,   0.058163,  -0.022034,   0.138566,  -0.035565,  -0.195047,  -0.040921,   0.051559,   0.139884,   0.185163,   0.040181,  -0.009547,  -0.091882,  -0.013914,   0.095011,  -0.668616,  -0.654269,
         0.187621,  -0.035470,  -0.058235,  -0.178911,  -0.125925,   0.004413,   0.086122,  -0.080336,  -0.022182,  -0.294598,   0.210180,  -0.024302,  -0.131706,  -0.048607,   0.233948,   0.169674,   0.059137,   0.099224,  -0.126823,   0.065104,   0.371333,   0.120501,   0.060225,  -0.047195,  -0.148267,   0.032933,   0.113151,  -0.091534,  -0.105765,  -0.281697,  -0.367012,   0.308484,
         0.106922,  -0.057442,  -0.177980,   0.043446,   0.267744,   0.246141,   0.330544,   0.134263,   0.178914,  -0.068275,   0.179611,  -0.135425,  -0.041303,   0.287206,   0.262735,   0.266863,   0.075705,   0.029051,   0.100935,  -0.208036,  -0.280955,  -0.463597,   0.676214,  -0.102109,  -0.271699,  -0.074582,  -0.113721,  -0.408830,   0.109137,  -0.035295,  -0.311395,  -0.037698,
        -0.930007,   0.093695,  -0.097819,  -0.265904,   1.660960,  -0.267720,   0.230430,   0.248184,  -0.410038,   0.003695,   0.215475,  -0.086233,  -0.408609,  -0.098302,  -0.116042,  -0.050924,   0.155151,   0.049134,   0.060946,  -0.098090,  -0.144596,  -0.149103,  -0.023771,  -0.017503,  -0.442276,  -0.133142,  -0.406230,   0.220088,   0.021005,   0.307201,  -0.017386,   0.112865,
         0.051375,  -0.026169,   0.133465,  -0.004653,  -0.064428,   0.045983,   0.006496,  -0.152070,   0.079297,   0.072736,   0.242771,   0.230785,   0.136711,   0.140105,   0.035470,   0.088285,   0.100711,  -0.076990,   0.104857,  -0.123388,   0.023143,  -0.013560,  -0.139989,  -0.306997,  -0.147520,  -0.030882,   0.054887,   0.394472,  -0.178868,   0.043345,   0.049353,  -0.152995,
         0.006975,   0.202674,   0.260836,  -0.050199,  -0.079336,  -0.027133,   0.105017,   0.086403,   0.098231,  -0.146836,   0.091657,   0.133007,  -0.163032,  -0.022666,  -0.943535,   0.032579,
        -0.127892,   0.078160,   0.160149,  -0.147463,  -0.010660,  -0.189444,   0.106893,   0.061279,  -0.291443,  -0.881388,   0.040639,   0.122439,   0.265180,  -0.027430,  -0.074823,   0.113869,  -0.080123,   0.173489,  -0.005362,   0.057500,   0.033952,  -0.166403,   0.131787,   0.080769,  -0.013398,  -0.051356,   0.063397,   0.018757,   0.198601,   0.302946,  -0.021174,   0.101972,
         0.116047,   0.069079,  -0.143815,  -0.226293,   0.033776,   0.108875,   0.287428,  -0.186925,  -0.495806,  -0.037176,   0.088371,  -0.133794,  -0.336659,   0.090090,  -0.073534,   0.221660,  -0.537627,  -0.044455,   0.273795,   0.137724,   0.218762,   0.036011,   0.105358,   0.097848,   0.117908,   0.208870,   0.367785,   0.352448,   0.173176,  -0.301019,   0.082356,  -0.097429,
         0.256861,   0.122628,  -0.149013,  -0.126139,   0.140491,  -0.364053,  -0.108427,  -0.006657,  -0.202086,  -0.136506,  -0.091906,  -0.265758,  -0.060797,  -0.195387,   0.068235,   0.077606,   0.218847,   0.072652,  -0.320342,   0.223329,   0.128282,  -0.160983,  -0.226912,  -0.030312,  -0.018763,   0.159087,   0.128365,   0.000280,  -0.120063,  -0.075391,  -0.234090,  -0.160632,
         0.049344,  -0.010826,   0.078841,  -0.065092,  -0.631680,  -0.001609,   0.152298,   0.041660,  -0.084512,   0.074978,   0.328667,   0.203292,  -0.198444,   0.060570,   0.058292,  -0.099259,   0.170690,  -0.010858,  -0.048473,  -0.198288,   0.084947,   0.072686,   0.168707,   0.130983,  -0.001611,   0.060619,  -0.017371,   0.029234,  -0.155526,  -0.209461,   0.039662,  -0.097320,
        -0.202544,  -0.406219,  -0.186073,  -0.010227,  -0.013328,   0.032968,  -0.132273,   0.162531,  -0.053585,  -0.021018,  -0.108592,  -0.220153,  -0.005958,  -0.026781,  -0.327067,   0.166761,   0.073031,   0.172414,  -0.059402,  -0.089859,  -0.125563,   0.095675,  -0.040814,   0.093590,   0.096994,   0.001774,   0.064980,   0.073480,   0.012685,  -0.065792,  -0.521767,  -0.895810,
         0.327304,   0.146418,   0.219667,  -0.103059,  -0.086806,   0.060665,   0.174447,   0.062466,   0.036478,  -0.096393,   0.225941,  -0.021202,   0.188961,   0.046267,   0.241907,  -0.163861,  -0.011517,   0.004073,   0.043845,   0.155630,   0.219079,   0.288418,   0.078059,   0.167857,   0.008340,  -0.000012,   0.155940,   0.042974,  -0.001457,  -0.175315,  -0.403800,   0.161692,
         0.001969,  -0.166789,  -0.371788,  -0.024306,   0.083745,   0.180126,  -0.432619,  -0.099525,   0.382097,  -0.043893,   0.169282,   0.216055,   0.173089,   0.145290,   0.089916,   0.176881,   0.091656,   0.283546,   0.133508,  -0.263400,  -0.162016,  -0.231624,   0.596596,  -0.032254,  -0.254398,  -0.134465,   0.150090,  -0.372614,  -0.026040,  -0.097711,   0.074892,  -0.246507,
        -0.867676,  -0.107155,   0.306097,   0.037499,   1.159280,  -0.089192,   0.275206,   0.289614,  -0.204567,   0.175035,  -0.022786,   0.042561,  -0.328816,  -0.322688,  -0.138656,   0.109464,   0.312649,   0.006847,   0.100085,  -0.093131,  -0.296489,  -0.212946,   0.159779,  -0.068719,  -0.392325,  -0.005596,  -0.473150,   0.107687,   0.089320,   0.096612,   0.150721,   0.034282,
         0.055748,   0.145167,  -0.012705,   0.119568,  -0.090075,   0.064892,   0.071833,   0.044568,  -0.212518,  -0.197021,   0.046167,   0.063582,   0.029293,   0.080612,   0.152416,  -0.090510,   0.086917,  -0.216244,   0.032081,  -0.101467,   0.119026,   0.025767,  -0.010587,  -0.279034,  -0.282055,   0.119321,  -0.105281,   0.080744,  -0.044808,  -0.072140,  -0.068520,  -0.070184,
        -0.051865,   0.031807,   0.183543,  -0.167711,  -0.081170,  -0.016794,  -0.088853,   0.396774,  -0.066148,  -0.168092,  -0.090771,  -0.044904,  -0.224594,   0.125600,  -0.783272,   0.251036
    };

    //! weights for layer 2
    const double ANN_WEIGHTS_CONTACT_SHEET_SHEET_LAYER_1[ 17] =
    {
        -0.137359,  -0.386343,  -0.502366,  -0.474267,   1.492680,   0.299190,  -0.554219,  -0.454494,   0.416579,   1.469120,   0.380457,  -0.617655,  -0.364990,  -0.370207,   0.264691,   0.477347,   0.258101
    };
    //! ANN CONTACT_SHEET_SHEET definition
    double ANN_CONTACT_SHEET_SHEET( const linal::Vector< double> &INP)
    {

      // declare variables
      int nora( 0), norb( 0), wei( 0);
      double *hid;

      // test net size
      BCL_Assert( INP.GetSize() == 303, "wrong input size!");

      // allocate memory
      linal::Vector< double> hidden[ 3];
      hidden[0] = linal::Vector< double>( 303);
      hidden[1] = linal::Vector< double>( 16);
      hidden[2] = linal::Vector< double>( 1);

      // normalize data
      hid = hidden[ 0].Begin();
      for( const double *inp = INP.Begin(); inp != INP.End(); inp++)
        ( *( hid++)) = ( *inp) * ANN_NORMALIZE_CONTACT_SHEET_SHEET_A[ nora++] + ANN_NORMALIZE_CONTACT_SHEET_SHEET_B[ norb++];

      // calculate network
      // calculate layer 1
      wei = 0;
      hid = hidden[1].Begin();
      for( size_t i = 0; i < 16; i++)
      {
        *hid = ANN_WEIGHTS_CONTACT_SHEET_SHEET_LAYER_0[ wei++];
        for( const double *inp = hidden[ 0].Begin(); inp != hidden[ 0].End(); inp++)
          ( *hid) += ANN_WEIGHTS_CONTACT_SHEET_SHEET_LAYER_0[ wei++] * ( *inp);
        *hid = double( 1.0) / ( double( 1.0) + exp( -( ( *hid)))); hid++;
      }

      // calculate layer 2
      wei = 0;
      hid = hidden[2].Begin();
      for( size_t i = 0; i < 1; i++)
      {
        *hid = ANN_WEIGHTS_CONTACT_SHEET_SHEET_LAYER_1[ wei++];
        for( const double *inp = hidden[ 1].Begin(); inp != hidden[ 1].End(); inp++)
          ( *hid) += ANN_WEIGHTS_CONTACT_SHEET_SHEET_LAYER_1[ wei++] * ( *inp);
        *hid = double( 1.0) / ( double( 1.0) + exp( -( ( *hid)))); hid++;
      }

      // denormalize data
      // end
      return hidden[ 2]( 0) * ANN_NORMALIZE_CONTACT_SHEET_SHEET_A[ nora] + ANN_NORMALIZE_CONTACT_SHEET_SHEET_B[ norb];
    }

  } // namespace contact
} // namespace bcl

//////////////////////////////////////////////////////////////////////////////////
// end of automatically built code                                              //
//////////////////////////////////////////////////////////////////////////////////
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
#include "contact/bcl_contact_calculate_correlations_mi.h"

// includes from bcl - sorted alphabetically
#include "align/bcl_align_alignment_interface.h"
#include "biol/bcl_biol_aa_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {

    //! @brief ConstAATypePair used to make code more readable
    typedef storage::Pair< biol::AAType, biol::AAType> AATypePair;

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CalculateCorrelationsMI::s_Instance
    (
      GetObjectInstances().AddInstance( new CalculateCorrelationsMI())
    );

    const size_t CalculateCorrelationsMI::s_LogBase = 20;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    CalculateCorrelationsMI::CalculateCorrelationsMI()
    {
    }

    //! @brief Clone function
    //! @return pointer to new CalculateCorrelationsMI
    CalculateCorrelationsMI *CalculateCorrelationsMI::Clone() const
    {
      return new CalculateCorrelationsMI( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CalculateCorrelationsMI::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    CorrelationMatrix CalculateCorrelationsMI::operator ()( const align::AlignmentInterface< biol::AABase> &ALIGNMENT_INTERFACE) const
    {
      // Acquire N length of MSA and number of sequences from ALIGNMENT_INTERFACE
      const size_t msa_length = ALIGNMENT_INTERFACE.GetSize();
      const size_t num_sequences = ALIGNMENT_INTERFACE.GetDepth();

      // Create a vector to store all the probabilities for base type x at position i
      storage::Vector< storage::Map< biol::AAType, double> > probability_vector( msa_length);

      // Create a matrix to hold each map of pair probabilities
      storage::SymmetricMatrix
      <
        storage::Map< AATypePair, double>
      >
        pair_probability_map_matrix( msa_length);

      // Get the list of assignments
      const util::ShPtrList< align::Assignment< biol::AABase> > assignments( ALIGNMENT_INTERFACE.GetAssignments());

      // Create vector to hold iterators pointing to each assignment list in the Alignment Interface
      // Using std::vector because storage::Vector tries to write out iterators and can't/not designed to
      std::vector< util::SiPtrList< const biol::AABase>::const_iterator> assignments_vector;
      assignments_vector.reserve( msa_length);

      // Iterate through assignment iterator vector setting each to beginning of corresponding members list
      util::ShPtrList< align::Assignment< biol::AABase> >::const_iterator
        assignments_itr( assignments.Begin()),
        assignments_itr_end( assignments.End());
      // Save end iterator for first position to keep track of how far to go
      util::SiPtrList< const biol::AABase>::const_iterator
        first_members_itr_end( ( *assignments_itr)->GetMembers().End());
      for
      (
        ;
        assignments_itr != assignments_itr_end;
        ++assignments_itr
      )
      {
        // Add the iterator for the current member list to the assignments_vector
        assignments_vector.push_back( ( *assignments_itr)->GetMembers().Begin());
      }

      // Iterate through each row represented by a vector of iterators pointing each to members in the same sequence
      // TODO: Use IS EMPTY?
      while( assignments_vector.front() != first_members_itr_end)
      {
        // Set sequence position index i which keeps track of positions in the alignment
        size_t i_index( 0);

        // Iterate through each position i
        for
        (
          std::vector< util::SiPtrList< const biol::AABase>::const_iterator>::iterator
            assignments_vector_itr( assignments_vector.begin()),
            assignments_vector_itr_end( assignments_vector.end());
          assignments_vector_itr != assignments_vector_itr_end;
          ++assignments_vector_itr, ++i_index
        )
        {
          // Save frequency data for type at this i in this row
          // Get current SiPtr and check to see if is NULL and set type appropriately
          util::SiPtr< const biol::AABase> aa_siptr( **assignments_vector_itr);
          biol::AAType temp_type_i( aa_siptr.IsDefined() ? aa_siptr->GetType() : biol::GetAATypes().GAP);

          // Check for key in map and insert if not present otherwise, increment count for type
          if( probability_vector( i_index).Has( temp_type_i))
          {
            ++( probability_vector( i_index)[ temp_type_i]); // #################### Not sure if this works!!!!!!
          }
          else
          {
            probability_vector( i_index).InsertElement
            (
              storage::Pair< biol::AAType, double>( temp_type_i, double( 1.0))
            );
          }

          // TODO: Change temp type stuff to a symmetric matrix using enums

          // Set sequence position index j which keeps track of positions in the alignment
          size_t j_index( 0);
          // Iterate through each position j where j is always > i
          std::vector< util::SiPtrList< const biol::AABase>::const_iterator>::iterator
            j_assignments_vector_itr( assignments_vector_itr);
          ++j_assignments_vector_itr;
          for( ; j_assignments_vector_itr != assignments_vector_itr_end; ++j_assignments_vector_itr, ++j_index)
          {
            // Get current SiPtr and check to see if is NULL and set type appropriately
            util::SiPtr< const biol::AABase> aa_siptr( **j_assignments_vector_itr);
            const biol::AAType temp_type_j( aa_siptr.IsDefined() ? aa_siptr->GetType() : biol::GetAATypes().GAP);

            // Create type ordered type pair to search map from types i and j
            AATypePair temp_pair;
            if( temp_type_i <= temp_type_j)
            {
              temp_pair = AATypePair( temp_type_i, temp_type_j);
            }
            else
            {
              temp_pair = AATypePair( temp_type_j, temp_type_i);
            }

            // Save a new map into pair_probability_map_matrix to keep track of x,y pair frequencies at position i,j
            if( pair_probability_map_matrix( i_index, j_index).IsEmpty())
            {
              pair_probability_map_matrix( i_index, j_index) = storage::Map< AATypePair, double >();
              pair_probability_map_matrix( i_index, j_index).InsertElement
              (
                storage::Pair< AATypePair, double>( temp_pair, double( 1.0))
              );
            }
            else
            {

              // Check to see if map at this i,j already has temp_pairing, increment x,y (x must always be <= y) or
              // add with count of 1 if not present yet
              if( pair_probability_map_matrix( i_index, j_index).Has( temp_pair))
              {
                ++( pair_probability_map_matrix( i_index, j_index)[ temp_pair]);
              }
              else
              {
                pair_probability_map_matrix( i_index, j_index).InsertElement
                  (
                    storage::Pair< AATypePair, double>( temp_pair, double( 1.0))
                  );
              }
            }
          } // End of iterator for position j
          // Increment position i iterator within the vector once done with all comparisons
          ( *assignments_vector_itr) = ++( *assignments_vector_itr);

        } // End of iterator for position i
      } // End of while loop which iterates through all rows in the ALIGNMENT_INTERFACE

      ////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////
      // Iterate through all counts and change all at current position into probabilities by dividing them by
      // num_sequences for all probabilities in the probability vector
      for( size_t pos_index( 0); pos_index < msa_length; ++pos_index)
      {
        for
        (
          storage::Map< biol::AAType, double>::iterator
            count_itr( probability_vector( pos_index).Begin()),
            count_itr_end( probability_vector( pos_index).End());
          count_itr != count_itr_end;
          ++count_itr
        )
        {
          probability_vector( pos_index)[ ( *count_itr).first] /= num_sequences;
        }
      }

      // Iterate through all positions in the matrix of pair probabilities and change all the probabilities to
      // decimal percentages by dividing by total num_sequences
      for( size_t i_index( 0); i_index < msa_length; ++i_index)
      {
        size_t j_index( i_index);
        ++j_index;
        for( ; j_index < msa_length; ++j_index)
        {
          for
          (
            storage::Map< AATypePair, double>::iterator
              map_itr( pair_probability_map_matrix( i_index, j_index).Begin()),
              map_itr_end( pair_probability_map_matrix( i_index, j_index).End());
            map_itr != map_itr_end;
            ++map_itr
          )
          {
            pair_probability_map_matrix( i_index, j_index)[ map_itr->first] /= num_sequences;
          }
        }
      }

      // Create correlation matrix to be filled in and returned
      CorrelationMatrix correlation_matrix;
      // Iterate through SymmetricMatrix with Maps of all probability pairs xi,yj- skips duplicate half across diagonal
      for( size_t i_index( 0); i_index < pair_probability_map_matrix.GetSize(); ++i_index)
      {
        // Set double sum_ij that keeps track of the sum at position i,j
        double sum_ij( 0.0);
        for( size_t j_index( 0); j_index <= i_index; ++j_index)
        {
          // Iterate through all frequency pairs at this i,j positionPredictorCorrelations
          for
          (
            storage::Map< AATypePair, double>::const_iterator
              pair_freq_itr( pair_probability_map_matrix( i_index, j_index).Begin()),
              pair_freq_itr_end( pair_probability_map_matrix( i_index, j_index).End());
            pair_freq_itr != pair_freq_itr_end;
            ++pair_freq_itr
          )
          {
            // Get types x and y at positions i and j respectively
            const biol::AAType &x_type( pair_freq_itr->first.First());
            const biol::AAType &y_type( pair_freq_itr->first.Second());

            //###### Should be possible to skip this for the sake of a slight speedup assuming that the above code
            //###### ensuring pairs are stored such that i <= j is not tampered with...
            // Create temp_pair for lookup such that in pair x,y x <= y
//            AATypePair temp_type_pair;
//            if( x_type <= y_type)
//            {
//              temp_type_pair = AATypePair( x_type, y_type);
//            }
//            else
//            {
//              temp_type_pair = AATypePair( y_type, x_type);
//            }
            // Get frequency for (xi,yj)
            const double prob_xi_yj( pair_freq_itr->second);

            // Get frequency for xi and for yj independently
            const double prob_xi( probability_vector( i_index).GetValue( x_type));
            const double prob_yj( probability_vector( j_index).GetValue( y_type));

            // Calculate for log_b(a) by calculating log(a)/log(b)
            const double log_result( log( ( prob_xi_yj) / ( prob_xi * prob_yj)) / log( double( s_LogBase)));
            sum_ij += prob_xi_yj * log_result;
          }
          // Save -sum_ij to correlation matrix
          correlation_matrix( i_index, j_index) = 0 - sum_ij;
        }
      }

      //!**********************************************************

//      // Iterate through probabilityies for all pairings xi,yj
//      for
//      (
//          storage::Vector< storage::Map< biol::AAType, double> >::const_iterator
//            i_prob_itr( probability_vector.Begin()),
//            prob_itr_end( probability_vector.End());
//          i_prob_itr != prob_itr_end;
//          ++i_prob_itr, ++i_index
//      )
//      {
//        // Set variable keeping track of the sum at position i,j
//        double sum_ij( 0.0);
//        // Set j_index as one position advanced from i_index
//        size_t j_index( i_index);
//        ++j_index;
//        // Set j position iterator one advanced from the i position iterator
//        storage::Vector< storage::Map< biol::AAType, double> >::const_iterator j_prob_iter = i_prob_itr;
//        ++j_prob_iter;
//        for
//        ( ; j_prob_iter != prob_itr_end; ++j_prob_iter, ++j_index)
//        {
//          // Set the first iterator for all the types in the probability map
//          for
//          (
//            storage::Map< biol::AAType, double>::const_iterator
//              k_type_itr( probability_vector(i_index).Begin()),
//              type_itr_end( probability_vector( i_index).End());
//            k_type_itr != type_itr_end;
//            ++k_type_itr
//          )
//          {
//            // Set the second iterator for all the types present in the probability map and start it one position advanced
//            storage::Map< biol::AAType, double>::const_iterator l_type_itr( k_type_itr);
//            ++l_type_itr;
//            for
//            ( ; l_type_itr != type_itr_end; ++l_type_itr)
//            {
//              double prob_xi( probability_vector( i_index).GetValue( k_type_itr->first));
//              double prob_yj( probability_vector( j_index).GetValue( l_type_itr->first));
//              // Create temp_pair for lookup such that in pair x,y x <= y
//              AATypePair temp_type_pair;
//              if ( k_type_itr->first <= l_type_itr->first)
//              {
//                temp_type_pair = AATypePair( k_type_itr->first, l_type_itr->first);
//              }
//              else
//              {
//                temp_type_pair = AATypePair( l_type_itr->first, k_type_itr->first);
//              }
//              double prob_xi_yj( pair_probability_map_matrix( i_index, j_index)[ temp_type_pair]);
//              // Calculate for log_b(a) by calculating log_2(a)/log_2(b)
//              double log_result( log2(( prob_xi_yj) / ( prob_xi * prob_yj)) / log2( double( 20.0)));
//              sum_ij +=  prob_xi_yj * log_result;
//            }
//          }
//        }
//
//        // Save the negative of the sum as the correlation for pair i,j in the correlation matrix
//        correlation_matrix( i_index, j_index) = sum_ij;
//      }

      return correlation_matrix;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CalculateCorrelationsMI::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &CalculateCorrelationsMI::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace contact
} // namespace bcl
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
#include "contact/bcl_contact_calculate_correlations_sm.h"

// includes from bcl - sorted alphabetically
#include "align/bcl_align_alignment_interface.h"
#include "biol/bcl_biol_aa_base.h"
#include "math/bcl_math_running_average_sd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {

  //////////
  // data //
  //////////

    //! Default is e_BLOSUM_62
    const score::AAAssignmentBLOSUM CalculateCorrelationsSM::s_SimilarityMatrix( score::AAAssignmentBLOSUM::e_BLOSUM_62);

    // TODO: Use assign score so you have a general score with type set as parameter

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CalculateCorrelationsSM::s_Instance
    (
      GetObjectInstances().AddInstance( new CalculateCorrelationsSM())
    );

    const double CalculateCorrelationsSM::s_GapCutoff( 0.1);

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    CalculateCorrelationsSM::CalculateCorrelationsSM() :
      m_WeightMap() // holds all the weights for all pairs of sequences k and l
    {
    }

    //! @brief Clone function
    //! @return pointer to new CalculateCorrelationsSM
    CalculateCorrelationsSM *CalculateCorrelationsSM::Clone() const
    {
      return new CalculateCorrelationsSM( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CalculateCorrelationsSM::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief calculates a correlation matrix if given an AlignmentInterface of type bio::AABase
    //! @details
    //! @param ALIGNMENT_INTERFACE is the MSA representation from which the CM is calculated
    //! @return Returns a correlation matrix of dimensions N by N where N is the length of the MSA, containing doubles

    // Currently returns nan for perfectly conserved positions
    CorrelationMatrix CalculateCorrelationsSM::operator()( const align::AlignmentInterface< biol::AABase> &ALIGNMENT_INTERFACE) const
    {
      // Acquire N length of MSA, number of sequences K from AlignmentInterface
      const size_t msa_length( ALIGNMENT_INTERFACE.GetSize());
      const size_t num_sequences( ALIGNMENT_INTERFACE.GetDepth());

      // Create vector to hold matrices with similarity values for each column of the MSA
      storage::Vector< linal::Matrix< double> > similarity_matrix_vector;
      similarity_matrix_vector.AllocateMemory( msa_length);

      // Create vector to hold all statistics as saved via RunningAverageSD< double> indexed by position pos_index
      storage::Vector< math::RunningAverageSD< double> > stat_mean_sd_vector;

      // Get the list of assignments
      const util::ShPtrList< align::Assignment< biol::AABase> > &assignments( ALIGNMENT_INTERFACE.GetAssignments());

      // Set sequence position index that keeps track of position in alignment
      size_t pos_index( 0);

      // Itr Alignment Interface across all nodes and get assignments
      for
      (
        util::ShPtrList< align::Assignment< biol::AABase> >::const_iterator
          assignment_itr( assignments.Begin()),
          assignment_itr_end( assignments.End());
        assignment_itr != assignment_itr_end;
        ++assignment_itr, ++pos_index
      )
      {
        // TODO: MOVE TO CONSTRUCTOR and SKIP ALLOCATION OF MEMORY to line 98
        // Create matrix of doubles to store similarity values to be calculated at this "column" of assignments and pushback into vector of matrices to keep track
        // TODO: Initialize all members of the matrix to zero, do I just add to the third field to accomplish that?
        similarity_matrix_vector.PushBack( linal::Matrix< double>( num_sequences, num_sequences));

        // Get SiPtrLst of members
        const util::SiPtrList< const biol::AABase> &members( ( *assignment_itr)->GetMembers());

        // Create DataSetStatisc
        math::RunningAverageSD< double> stat_mean_col_current;

        // Set sequence position k_index, and l_index
        size_t k_index( 0);
        size_t l_index( 0);

        for
        (
          util::SiPtrList< const biol::AABase>::const_iterator
            k_assign_itr( members.Begin()),
            assign_itr_end( members.End());
          k_assign_itr != assign_itr_end;
          ++k_assign_itr, ++k_index
        )
        {
          // Set l index and l iterator such that sequences aren't compared to self and all pairs are accounted for once
          l_index = k_index + 1;
          // Set l_assign_itr such that it is != to k_assign_itr
          util::SiPtrList< const biol::AABase>::const_iterator l_assign_itr( k_assign_itr);
          ++l_assign_itr;
          for( ; l_assign_itr != assign_itr_end; ++l_assign_itr, ++l_index)
          {
            // Calculate similarity value for Sikl using s_SimilarityMatrix
            // this code is ridiculous
            // Check for null pointers since GAPs are stored in alignments as null pointers and not AABases of type GAP
            const util::SiPtr< const biol::AABase> &aa_k_siptr( *k_assign_itr);
            const biol::AAType type_k( aa_k_siptr.IsDefined() ? aa_k_siptr->GetType() : biol::GetAATypes().GAP);
            // Check for null pointers since GAPs are stored in alignments as null pointers and not AABases of type GAP
            const util::SiPtr< const biol::AABase> &aa_l_siptr( *l_assign_itr);
            const biol::AAType type_l( aa_l_siptr.IsDefined() ? aa_l_siptr->GetType() : biol::GetAATypes().GAP);
            const double similarity_kl
            (
              s_SimilarityMatrix.Probability
              (
                s_SimilarityMatrix.e_BLOSUM_62,
                type_k,
                type_l
              )
            );

            // for current stat_mean_col_current add the value so that later it can computer mean and standard deviation
            stat_mean_col_current += similarity_kl;

            // add similarity value for Sikl to position k,l in the current similarity matrix
            similarity_matrix_vector( pos_index)( k_index, l_index) = similarity_kl;

            // add to weight if non-identical elements present in k and l in position k,l in the weights map keyed by
            // k,l
            const storage::Pair< size_t, size_t> temp_pair( k_index, l_index);
            // check to see if this pair is already in the m_WeightMap if no, insert pair along with current count for non-identical elements
            double non_identical_count( ( type_k != type_l));
            // TODO: Change map to symmetric matrix
            // TODO: check logic
            if( m_WeightMap.Has( temp_pair) && non_identical_count) // Only executed if non_identical_count is != 0
            {
              m_WeightMap[ temp_pair] += non_identical_count;
            }
            else
            {
              m_WeightMap.InsertElement( storage::Pair< storage::Pair< size_t, size_t>, double>( temp_pair, non_identical_count));
            }
          }
        }
        // pushback most recently finished statistics object into statistics vector
        stat_mean_sd_vector.PushBack( stat_mean_col_current);
      } // end alignment iterator for loop

      // TODO: Sum to 1 function in math matrix, change this and maybe make a linal::symmetric matrix

      // Normalize weights in m_WeightsMap such that they sum to 1
      // First iterate through Map and sum up all values in weights_sum
      double weights_sum( 0.0);
      for
      (
        storage::Map< storage::Pair< size_t, size_t>, double>::const_iterator
          map_itr( m_WeightMap.Begin()),
          map_itr_end( m_WeightMap.End());
        map_itr != map_itr_end;
        ++map_itr
      )
      {
        weights_sum += map_itr->second;
      }

      // Iterate through Map again, now saving the values as value/total sum of values to normalize
      for
      (
        storage::Map< storage::Pair< size_t, size_t>, double>::iterator
          map_itr( m_WeightMap.Begin()),
          map_itr_end( m_WeightMap.End());
        map_itr != map_itr_end;
        ++map_itr
      )
      {
        map_itr->second /= weights_sum;
//        m_WeightMap[ map_itr->first] = ( m_WeightMap[ map_itr->first] / weights_sum);
      }

      // Create CorrelationMatrix to fill up with values
      CorrelationMatrix correlation_matrix( msa_length);

      // Iterate through all pairs of k l summing values of formula for all pairings i j
      for( size_t pos_i( 0); pos_i < msa_length; ++pos_i)
      {
        for( size_t pos_j( pos_i + 1); pos_j < msa_length; ++pos_j)
        {
          // Keep track of sum of all pairs of kl
          double sum_for_ij( 0.0);

          // Iterate through all possible pairs of k and l applying the formula below
          // For all pairs of i and j where i != j
          // Correlation at Rij =
          // 1/N^2 *
          // Sum(for all pairs of sequences k and l where k != l)
          // (Weight(k,l)  *  ( similarity at i between k and l - average similarity at i) *
          // ( similarity at j between k and l - average similarity at j)) divided by
          // (standard deviation of similarities at column i)*(standard deviation of similarities at j)
          for( size_t seq_k( 0); seq_k < num_sequences; ++seq_k)
          {
            for( size_t seq_l( seq_k + 1); seq_l < num_sequences; ++seq_l)
            {
              // TODO: statistics function for this score, if not there implement and call here
              sum_for_ij += ( GetWeight( seq_k, seq_l)
                  * ( similarity_matrix_vector( pos_i)( seq_k, seq_l) - stat_mean_sd_vector( pos_i).GetAverage())
                  * ( similarity_matrix_vector( pos_j)( seq_k, seq_l) - stat_mean_sd_vector( pos_j).GetAverage()))
                  / ( stat_mean_sd_vector( pos_i).GetStandardDeviation() * stat_mean_sd_vector( pos_j).GetStandardDeviation());
            }
          }
          // Divide the sum of all pairs kl by the total sequence length squared and save that to r_ij
          correlation_matrix( pos_i, pos_j) = ( sum_for_ij)/( math::Sqr( msa_length));
          // Correlation matrix is diagonally symmetrical
          correlation_matrix( pos_j, pos_i) = correlation_matrix( pos_i, pos_j);
        }
      }
      // Return filled in correlation_matrix
      return correlation_matrix;
    }

    //! @brief Checks to see if param CalculateCorrelationsSM object is equal to this one
    //! @param CALCULATE_CORRELATIONS reference to a CalculateCorrelationsSM object
    //! @return bool true for all members are equal and false otherwise
    bool CalculateCorrelationsSM::operator==( const CalculateCorrelationsSM &CALCULATE_CORRELATIONS) const
    {
      if
      (
        s_GapCutoff == CALCULATE_CORRELATIONS.s_GapCutoff &&
          s_SimilarityMatrix.GetTableType() == CALCULATE_CORRELATIONS.s_SimilarityMatrix.GetTableType() &&
          m_WeightMap.GetSize() == CALCULATE_CORRELATIONS.m_WeightMap.GetSize()
      )
      {
        for
        (
          storage::Map< storage::Pair< size_t, size_t>, double>::const_iterator
            map_itr_a( m_WeightMap.Begin()),
            map_itr_b( CALCULATE_CORRELATIONS.m_WeightMap.Begin()),
            map_itr_end( m_WeightMap.End());
          map_itr_a != map_itr_end;
          ++map_itr_a, ++map_itr_b
        )
        {
          if
          (
            map_itr_a->first.First() != map_itr_b->first.First() ||
              map_itr_a->first.Second() != map_itr_b->first.Second() ||
              map_itr_a->second != map_itr_b->second
          )
          {
            return false;
          }
        }
        return true;
      }
      // else
      return false;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CalculateCorrelationsSM::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_WeightMap, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &CalculateCorrelationsSM::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_WeightMap, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Takes two sequences from the MSA given to CalculateWeights and returns the weight value for that pair
    //! @param K_SEQ_INDEX from the MSA for sequence k
    //! @param L_SEQ_INDEX from the MSA for sequence l
    //! @return Returns a double weight value
    double CalculateCorrelationsSM::GetWeight( const size_t K_SEQ_INDEX, const size_t L_SEQ_INDEX) const
    {
      // check to make sure that K_SEQ_INDEX is < L_SEQ_INDEX and if not swap them and use them as the keys for map lookup
      if( K_SEQ_INDEX < L_SEQ_INDEX)
      {
        return m_WeightMap.GetValue( storage::Pair< size_t, size_t>( K_SEQ_INDEX, L_SEQ_INDEX));
      }

      return m_WeightMap.GetValue( storage::Pair< size_t, size_t>( L_SEQ_INDEX, K_SEQ_INDEX));
    }

  } // namespace contact
} // namespace bcl
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
#include "contact/bcl_contact_correlation_matrix.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_running_average.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CorrelationMatrix::s_Instance
    (
      GetObjectInstances().AddInstance( new CorrelationMatrix())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Constructor taking in dimension
    CorrelationMatrix::CorrelationMatrix( const size_t DIMENSIONS, const double VALUE, const bool REPLACE) :
      storage::SymmetricMatrix< double>( DIMENSIONS),
      m_Value( VALUE)
    {
      // TODO: make another class to handle this case        m_ReplaceValues( REPLACE)
    }

    //! @brief Clone function
    //! @return pointer to new CorrelationMatrix
    CorrelationMatrix *CorrelationMatrix::Clone() const
    {
      return new CorrelationMatrix( *this);
    }

  /////////////////
  // data access //
  /////////////////

    // TODO: Make example test for GetAverage

    //! @brief Returns the average for the symmetric matrix
    //! @return double corresponding to average correlation value of matrix excluding NaN and the diagonal
    double CorrelationMatrix::GetAverage() const
    {
      math::RunningAverage< double> average;
      for( size_t i( 0), size( GetSize()); i < size; ++i)
      {
        for( size_t j( i + 1); j < size; ++j)
        {
          const double value( storage::SymmetricMatrix< double>::operator()( i, j));
          if( util::IsDefined( value))
          {
            average += value;
          }
        }
      }

      return average;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CorrelationMatrix::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! Returns value but substitutes in m_Value if NaN is encountered
    //! TODO: Will this allow changing of NaN?  Should I try and stop this?
    const double &CorrelationMatrix::operator()( const size_t I_POS, const size_t J_POS) const
    {
      const double &value( storage::SymmetricMatrix< double>::operator()( I_POS, J_POS));
      return util::IsDefined( value) ? value : m_Value;
    }

    //! operator( POS) return reference to changeable element at POS
    double &CorrelationMatrix::operator()( const size_t I_POS, const size_t J_POS)
    {
      // Call const method and return
      return storage::SymmetricMatrix< double>::operator()( I_POS, J_POS);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CorrelationMatrix::Read( std::istream &ISTREAM)
    {
      // read members
      storage::SymmetricMatrix< double>::Read( ISTREAM);
      io::Serialize::Read( m_Value, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &CorrelationMatrix::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      storage::SymmetricMatrix< double>::Write( OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Value, OSTREAM, INDENT);              // write members

      // return the stream
      return OSTREAM;
    }

  } // namespace contact
} // namespace bcl
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
#include "align/bcl_align_alignment_node.h"
#include "align/bcl_align_handler_pir.h"
#include "contact/bcl_contact_correlation_matrix.h"
#include "contact/bcl_contact_correlation_storage_file.h"
#include "io/bcl_io_directory_entry.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {
  //////////
  // data //
  //////////

    //! format to convert Key to string                                                       // CHECK FORMAT
    const util::Format CorrelationStorageFile::s_KeyToString
    (
      util::Format().W( 6).R().Fill( '0')
    );

    //! default file extension for objects
    const std::string CorrelationStorageFile::s_FileExtension( "bcl");

    //! default prefix for filename of stored model information
    const std::string CorrelationStorageFile::s_FilePrefix( "ampair");

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CorrelationStorageFile::s_Instance
    (
      GetObjectInstances().AddInstance( new CorrelationStorageFile())
    );

    //! Default delimiter to be used for reading/writing alignments to file
    const char CorrelationStorageFile::s_Delim( '$');

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  /////////////////
  // data access //
  /////////////////
    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CorrelationStorageFile::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief initialize the correlation storage
    //! @param INITIALIZER encodes where data is stored
    //! @param INIT_FLAG flag for type of initialization
    //! @return true if initialize was successful
    bool CorrelationStorageFile::Initialize( const std::string &INITIALIZER, const InitializerType INIT_FLAG)
    {
      BCL_MessageStd( "Initializer: " + INITIALIZER);
      util::ShPtr< io::Directory> sp_dir( new io::Directory( INITIALIZER));

      // check if directory exists
      const bool dir_exists( sp_dir->DoesExist());

      switch( INIT_FLAG)
      {
        // attach to existing storage - directory needs to exist
        case e_Attach:
        {
          if( dir_exists)
          {
            m_Directory = sp_dir;
          }
          return dir_exists;
        }
        // create storage, directory does need to not exist
        case e_Create:
        {
          // if it exists, it will not be overwritten
          if( dir_exists)
          {
            BCL_MessageCrt( "cannot create " + INITIALIZER + " since the directory already exists!");
            return false;
          }
          const bool make_success( sp_dir->Make());
          if( make_success)
          {
            m_Directory = sp_dir;
          }
          return make_success;
        }
        // overwrite directory, when it exist, will be overwritten
        case e_Overwrite:
        {
          if( dir_exists)
          {
            // clear out any existing model files; ignore other files
            m_Directory = sp_dir;

            // get the keys of all current models
            storage::Vector< std::string> keys( GetAllKeys());

            bool clear_success( true);

            // iterate over the keys and delete them
            for
            (
              storage::Vector< std::string>::const_iterator itr( keys.Begin()), itr_end( keys.End());
              itr != itr_end && clear_success;
              ++itr
            )
            {
              io::DirectoryEntry current_entry( sp_dir, s_FilePrefix + *itr + io::File::GetExtensionDelimiter() + s_FileExtension);
              clear_success &= current_entry.Remove();
            }
            if( !clear_success)
            {
              BCL_MessageCrt( "cannot clear " + INITIALIZER + " to overwrite!");
              return clear_success;
            }
            m_Directory = sp_dir;

            // initalize as create
            return clear_success;
          }
          else // directory does not exist
          {
            const bool make_success( sp_dir->Make());
            if( make_success)
            {
              m_Directory = sp_dir;
            }
            return make_success;
          }
        }
        default:
        {
          BCL_Exit( "unknown InitializerType supplied", -1);
          break;
        }
      }

      // end
      return true;
    }

    //! @brief number of data items in source
    //! @param SOURCE source of data file
    //! @return number of MSA correlation matrix pairs in storage List at SOURCE
    size_t CorrelationStorageFile::GetSize() const
    {
      return GetAllKeys().GetSize();
    }

    //! @brief get all keys for given source
    //! @return all keys of given source
    storage::Vector< std::string> CorrelationStorageFile::GetAllKeys() const
    {
      // keys
      storage::Vector< std::string> keys;

      // Check that directory exists
      if( m_Directory.IsDefined())
      {
        // directory content
        const storage::List< io::DirectoryEntry> dir_entries
        (
          m_Directory->ListEntries( io::Directory::e_File, s_FilePrefix, "." + s_FileExtension)
        );

        // iterate over all entries
        storage::List< io::DirectoryEntry>::const_iterator itr( dir_entries.Begin()), itr_end( dir_entries.End());
        for( ; itr != itr_end; ++itr)
        {
          // filename without extension
          std::string key( io::File::RemoveLastExtension( itr->GetName()));

          // extract the key and add it to keys vector
          key.erase( 0, s_FilePrefix.length());
          keys.PushBack( key);
        }
      }

      return keys;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get AlignmentMatrixPair
    //! @param SOURCE source of molecule eg. given AlignmentMatrixPairs
    //! @param KEY key identifier for specific molecule in given source
    //! @return shptr to AlignmentMatrixPair of interest and shptr to NULL if not found
    util::ShPtr< CorrelationStorageFile::AlignmentMatrixPair> CorrelationStorageFile::Retrieve( const std::string &KEY) const
    {
      // Create filename using key
      std::string filename( StorageFilename( KEY));

      // Check for file
      if( !io::DirectoryEntry( filename).DoesExist())
      {
        // BCL Message for verbose
        BCL_MessageVrb( "File :" + filename + " does not exist");

        // Return ShPtr to NULL if failed
        return util::ShPtr< AlignmentMatrixPair>( 0);
      }

      // BCL Message for verbose
      BCL_MessageVrb( "File :" + filename + " found");
      // Open file, create object, and read in
      io::IFStream read;
      io::File::MustOpenIFStream( read, filename);
      std::string alignment;
      CorrelationMatrix correlation_matrix;
      BCL_MessageDbg( "Reading alignment portion " + filename);
      io::Serialize::Read( alignment, read);
      BCL_MessageDbg( "Reading alignment portion " + filename + " done");
      BCL_MessageDbg( "Reading correlation matrix " + filename);
      io::Serialize::Read( correlation_matrix, read);
      BCL_MessageDbg( "Reading correlation matrix " + filename + " done");

      // Create alignment handler to read in string portion into alignment object
      align::HandlerPIR< biol::AABase> pir_handler( s_Delim);
      std::stringstream stream( alignment);
      BCL_MessageDbg( "Reading alignment " + filename);
      BCL_MessageDbg( "Read alignment :" + alignment);
      util::ShPtr< align::AlignmentNode< biol::AABase> > alignment_ptr
      (
        pir_handler.ReadAlignment( stream, biol::AASequence())
      );
      BCL_MessageDbg( "Reading alignment " + filename + " done");

      // Translate input into AlignmentMatrixPair object
      util::ShPtr< CorrelationStorageFile::AlignmentMatrixPair> am_pair_ptr( new AlignmentMatrixPair( *alignment_ptr, correlation_matrix));

      // Close file
      io::File::CloseClearFStream( read);

      // Return ShPtr to AMPair
      return am_pair_ptr;
    }

    //! @brief get ensemble of AlignmentMatrixPair
    //! @return shptr list of AlignmentMatrixPair from given source directory
    util::ShPtrList< CorrelationStorageFile::AlignmentMatrixPair> CorrelationStorageFile::RetrieveEnsemble() const
    {
      // Get keys
      storage::Vector< std::string> key_vector( GetAllKeys());

      // Return retrieve for all keys
      return RetrieveEnsemble( key_vector);
    }

    //! @brief get ensemble of AlignmentMatrixPair for given keys
    //! @param KEYS vector of keys
    //! @return shptr list of molecules from given source
    util::ShPtrList< CorrelationStorageFile::AlignmentMatrixPair> CorrelationStorageFile::RetrieveEnsemble
    (
      const storage::Vector< std::string> &KEYS
    ) const
    {
      // Retrieve given keys and save them in the list
      util::ShPtrList< AlignmentMatrixPair> am_pair_ptr_list;
      storage::Vector< std::string>::const_iterator iter( KEYS.Begin());
      storage::Vector< std::string>::const_iterator iter_end( KEYS.End());
      for( ; iter != iter_end; ++iter)
      {
        am_pair_ptr_list.PushBack( Retrieve( *iter));
      }

      // Return list of results
      return am_pair_ptr_list;
    }

//    //! @brief get ensemble of AlignmentMatrixPair for given keys which does not map to keys, following plain numbering
//    //! @param RANGE range of AlignmentMatrixPair
//    //! @return shptr list of AlignmentMatrixPair from given source
//    util::ShPtrList< AlignmentMatrixPair> CorrelationStorageFile::RetrieveEnsemble
//    (
//      const math::Range< size_t> &RANGE
//    ) const
//    {
//        // IMPLEMENT
//    }

    //! @brief store AlignmentMatrixPairs
    //! @param ALIGNMENT_MATRIX_PAIR to store
    //! @return key string associated with target sequence - pdb sequence ID, otherwise empty string
    std::string CorrelationStorageFile::Store( const AlignmentMatrixPair &ALIGNMENT_MATRIX_PAIR)
    {
      // key of stored AlignmentMatrixPair
      const std::string key( ALIGNMENT_MATRIX_PAIR.First().GetSequenceIds().FirstElement());

      // Store with key and return key
      if( CorrelationStorageFile::Store( ALIGNMENT_MATRIX_PAIR, key))
      {
        return key;
      }
      else
      {
        // else return empty string
        return std::string();
      }
    }

    //! @brief store AlignmentMatrixPairs as pairs of alignment string in wrapper and correlation matrix
    //! @param ALIGNMENT_MATRIX_PAIR to store
    //! @param KEY key under which AlignmentMatrixPair was stored
    //! @return true if store was successful false if unsuccessful or file already exists
    bool CorrelationStorageFile::Store( const AlignmentMatrixPair &ALIGNMENT_MATRIX_PAIR, const std::string &KEY)
    {
      // OFstream for writing into bcl file
      io::OFStream write;

      // filename for writing into
      std::string filename( StorageFilename( KEY));

      // Write out filename if verbose
      BCL_MessageVrb( "File path for storage: " + filename);

      // Check for file already of that pdb_id and prefix
      if( io::DirectoryEntry( filename).DoesExist())
      {
        return false;
      }

      // open file for reading
      io::File::MustOpenOFStream( write, filename);

      // translate Alignment given by node to string stream for storage
      align::HandlerPIR< biol::AABase> pir_handler( s_Delim);
      std::stringstream alignment_string;
      pir_handler.WriteAlignment( alignment_string, ALIGNMENT_MATRIX_PAIR.First());

      // write AlignmentMatrixPair
      io::Serialize::Write( alignment_string.str(), write) << '\n';
      io::Serialize::Write( ALIGNMENT_MATRIX_PAIR.Second(), write);
      // close file stream
      io::File::CloseClearFStream( write);

      return true;
    }

    //! @brief store ensemble of AlignmentMatrixPair
    //! @param ENSEMBLE shptr list of AlignmentMatrixPair
    //! @return vector of keys
    storage::Vector< std::string> CorrelationStorageFile::Store( const util::ShPtrList< AlignmentMatrixPair> &ENSEMBLE)
    {
      // Create vector of return keys
      storage::Vector< std::string> key_vector;

      // Iterate through ShPtrList of AlignmentMatrixPairs and store each
      util::ShPtrList< AlignmentMatrixPair>::const_iterator itr( ENSEMBLE.Begin());
      util::ShPtrList< AlignmentMatrixPair>::const_iterator itr_end( ENSEMBLE.End());
      for( ; itr != itr_end; ++itr)
      {
        // Store AlignmentMatrixPair and store returning key in vector to be returned
        key_vector.PushBack( Store( **itr));
      }

      // return vector of keys some of which may be empty strings
      return key_vector;
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CorrelationStorageFile::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Directory, ISTREAM);
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &CorrelationStorageFile::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Directory, OSTREAM, INDENT);
      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace contact
} // namespace bcl
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
#include "contact/bcl_contact.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_range.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    const std::string &GetNamespaceIdentifier()
    {
      static const std::string *s_namespace_name( new std::string( util::ExtractNamespaceIdentifier( __PRETTY_FUNCTION__)));
      return *s_namespace_name;
    }

    //! @brief return default sequence separation range for contacts
    //! @return default sequence separation range contacts
    const math::Range< size_t> &GetDefaultSequenceSeparationRange()
    {
      // initialize static range
      static const math::Range< size_t> s_range( 6, std::numeric_limits< size_t>::max());

      // end
      return s_range;
    }

    //! @brief return sequence separation range for short-range contacts
    //! @return sequence separation range for short-range contacts
    const math::Range< size_t> &GetDefaultSequenceSeparationShortRange()
    {
      // initialize static range
      static const math::Range< size_t> s_range_short( 6, 11);

      // end
      return s_range_short;
    }

    //! @brief return sequence separation range for mid-range contacts
    //! @return sequence separation range for mid-range contacts
    const math::Range< size_t> &GetDefaultSequenceSeparationMidRange()
    {
      // initialize static range
      static const math::Range< size_t> s_range_mid( 12, 23);

      // end
      return s_range_mid;
    }

    //! @brief return sequence separation range for long-range contacts
    //! @return sequence separation range for long-range contacts
    const math::Range< size_t> &GetDefaultSequenceSeparationLongRange()
    {
      // initialize static range
      static const math::Range< size_t> s_range_long( 24, std::numeric_limits< size_t>::max());

      // end
      return s_range_long;
    }

  } // namespace contact
} // namespace bcl
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
#include "contact/bcl_contact_data.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {
  //////////
  // data //
  //////////

    //! single instance
    const util::SiPtr< const util::ObjectInterface> Data::s_Instance( GetObjectInstances().AddInstance( new Data()));

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Data::Data() :
      m_Data( Types::s_NumberValidTypes, util::GetUndefined< double>()), // set all probabilities to undefined
      m_Merged( util::GetUndefined< double>())
    {
    }

    //! @brief Clone function
    //! @return pointer to new Data
    Data *Data::Clone() const
    {
      return new Data( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Data::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the merged contact prediction
    //! @return a const reference to the prediction value
    const double &Data::MergedPrediction() const
    {
      return m_Merged;
    }

    //! @brief return the merged contact prediction
    //! @return a changeable reference to the prediction value
    double &Data::MergedPrediction()
    {
      return m_Merged;
    }

    //! @brief returns if all contact prediction data is defined
    //! @return if data is defined
    bool Data::IsDefined() const
    {
      return m_Data.IsDefined() && util::IsDefined( m_Merged);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief swaps prediction values at positions CONTACT_TYPE_A and CONTACT_TYPE_B
    //! @param CONTACT_TYPE_A first type position to swap
    //! @param CONTACT_TYPE_B second type position to swap
    void Data::Swap( const Type CONTACT_TYPE_A, const Type CONTACT_TYPE_B)
    {
      std::swap( m_Data( CONTACT_TYPE_A), m_Data( CONTACT_TYPE_B));
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief returns the prediction for contact type CONTACT_TYPE
    //! @param CONTACT_TYPE the contact type for which to return the prediction
    //! @return a const reference to the prediction value
    const double &Data::operator[]( const Type &CONTACT_TYPE) const
    {
      return m_Data( CONTACT_TYPE);
    }

    //! @brief returns the prediction for contact type CONTACT_TYPE
    //! @param CONTACT_TYPE the contact type for which to return the prediction
    //! @return a changeable reference to the prediction value
    double &Data::operator[]( const Type &CONTACT_TYPE)
    {
      return m_Data( CONTACT_TYPE);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Data::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Data, ISTREAM);
      io::Serialize::Read( m_Merged, ISTREAM);

      return ISTREAM; // return the stream
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Data::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Data, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_Merged, OSTREAM, INDENT);

      return OSTREAM; // return the stream
    }

  } // namespace contact
} // namespace bcl
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
#include "contact/bcl_contact_map.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_geometry_packing.h"
#include "biol/bcl_biol_atom.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {

  //////////
  // data //
  //////////

    //! identifier string to be used for prediction maps
    const std::string Map::s_Identifier = "CHAINS";

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Map::Map() :
      storage::ObjectNDHashMap< 2, biol::AAData, storage::Pair< linal::Vector< size_t>, bool> >(),
      m_Chains(),
      m_Boundary( 8)
    {
    }

    //! @brief construct contact map from a chain
    //! @param CHAIN Chain for which sequence contacts are going to be deteced
    //! @param BOUNDARY number of residues to be excluded from both sides
    Map::Map
    (
      const util::ShPtr< assemble::Chain> &CHAIN,
      const int BOUNDARY
    ) :
      storage::ObjectNDHashMap< 2, biol::AAData, storage::Pair< linal::Vector< size_t>, bool> >(),
      m_Chains( 1, CHAIN),
      m_Boundary( BOUNDARY)
    {
      FillMap();
    }

    //! @brief construct contact map from a protein model( intra-sequence contacts and inter-sequence contacts),
    //! @param PROTEIN_MODEL Protein Model for which inter and intra sequence contacts are going to be deteced
    //! @param BOUNDARY number of residues to be excluded from both sides
    Map::Map
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      const int BOUNDARY
    ) :
      storage::ObjectNDHashMap< 2, biol::AAData, storage::Pair< linal::Vector< size_t>, bool> >(),
      m_Chains( PROTEIN_MODEL.GetChains()),
      m_Boundary( BOUNDARY)
    {
      FillMap();
    }

    //! @brief virtual copy constructor
    Map *Map::Clone() const
    {
      return new Map( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Map::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns stored chains
    //! @return stored chains
    const util::ShPtrVector< assemble::Chain> &Map::GetChains() const
    {
      return m_Chains;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief iterates over the stored chains and finds the one with matching CHAIN_ID
    //! @param CHAIN_ID chain id of the chain that is beings searched
    //! @return Chain  with the searched chain id
    const util::ShPtr< assemble::Chain> &Map::GetChain( const char CHAIN_ID) const
    {
      // iterate over every chain stored
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator chain_itr( m_Chains.Begin()),
          chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // if the sequence id matches return it
        if( ( *chain_itr)->GetSequence()->GetChainID() == CHAIN_ID)
        {
          return *chain_itr;
        }
      }

      // else Exit
      BCL_Exit( "No chain with the provided chain id \'" + util::Format()( CHAIN_ID) + "\' is not stored!!", -1);
      return m_Chains.FirstElement();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Returns contacts vector for the provided AA Data pointers
    //! @param AA_DATA_POINTERS pair of pointers to AAData of amino acids of interest
    //! @return the pair of vector of contact vectors and the merged prediction
    const storage::Pair< linal::Vector< size_t>, bool> &Map::GetContactVector
    (
      const storage::VectorND< 2, util::SiPtr< const biol::AAData> > &AA_DATA_POINTERS
    ) const
    {
      // initialize undefined predictions vector
      static const storage::Pair< linal::Vector< size_t>, bool> s_undefined_contact_vector
      (
        linal::Vector< size_t>( Types::s_NumberValidTypes, util::GetUndefined< size_t>()),
        false
      );

      // search in the hash map
      storage::HashMap< size_t, storage::Pair< linal::Vector< size_t>, bool> >::const_iterator itr
      (
        Find( AA_DATA_POINTERS)
      );

      // if not found return undefined
      if( itr == End())
      {
        return s_undefined_contact_vector;
      }

      // end
      return itr->second;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Map::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &Map::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

    //! @brief helper function to read predictions from a file
    //! @param ISTREAM input stream to be read from
    //! @return input stream from which the map was read
    std::istream &Map::ReadMap( std::istream &ISTREAM)
    {
      // header that is going to be searched
      std::string header;

      // while reading header and not end of file yet
      while( ISTREAM >> header && !ISTREAM.eof())
      {
        // if header matches the identifier
        if( header == s_Identifier)
        {

          // read the two chain ids
          char chain_id_a, chain_id_b;
          ISTREAM >> chain_id_a >> chain_id_b;

          // check the chain_id_a exists
          BCL_Assert
          (
            GetChain( chain_id_a).IsDefined(),
            "No chain is stored in the map with given chain id: " + std::string( 1, chain_id_a)
          );

          // check the chain_id_b exists
          BCL_Assert
          (
            GetChain( chain_id_b).IsDefined(),
            "No chain is stored in the map with given chain id: " + std::string( 1, chain_id_b)
          );

          // read the contact map for two chains with the provided chain ids
          ReadContacts
          (
            ISTREAM,
            *GetChain( chain_id_a),
            *GetChain( chain_id_b)
          );
        }
      }

      // end
      return ISTREAM;
    }

    //! @brief helper function to write predictions to a file
    //! @param OSTREAM output stream to write to
    //! @param WRITE_ONLY_CONTACTS boolean value that determines whether non contacting residue couples are written
    //! @return output stream which was written to
    std::ostream &Map::WriteMap
    (
      std::ostream &OSTREAM,
      const bool WRITE_ONLY_CONTACTS
    ) const
    {
      // iterate over every sequence
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator chain_itr_a( m_Chains.Begin()),
        chain_itr_end( m_Chains.End());
        chain_itr_a != chain_itr_end;
        ++chain_itr_a
      )
      {
        // versus every other sequence
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator chain_itr_b( m_Chains.Begin());
          chain_itr_b != chain_itr_end;
          ++chain_itr_b
        )
        {
          // output the header
          OSTREAM << s_Identifier << ' '
                  << util::Format()( (*chain_itr_a)->GetChainID()) << ' '
                  << util::Format()( (*chain_itr_b)->GetChainID()) << '\n';

          // write predictions
          WriteContacts( OSTREAM, **chain_itr_a, **chain_itr_b, WRITE_ONLY_CONTACTS);
        }
      }
      return OSTREAM;
    }

    //! @brief reads contacts for the provided chains
    //! @param ISTREAM input stream
    //! @param CHAIN_A first chain
    //! @param CHAIN_B second chain
    //! @return istream from which was read
    std::istream &Map::ReadContacts
    (
      std::istream &ISTREAM,
      const assemble::Chain &CHAIN_A,
      const assemble::Chain &CHAIN_B
    )
    {
      // initialize sequence id and AAType pairs and the vector hold predictions
      int seq_id_a, seq_id_b;
      char type_a, type_b;
      storage::Pair< linal::Vector< size_t>, bool> contact_vector
      (
        linal::Vector< size_t>( Types::s_NumberValidTypes, util::GetUndefined< size_t>()),
        false
      );

      // while it is possible to read
      while( ISTREAM >> seq_id_a >> type_a >> seq_id_b >> type_b && !ISTREAM.eof())
      {
        // make the sure sequences match
        BCL_Assert
        (
          CHAIN_A.GetSequence()->GetAA( seq_id_a - 1)->GetType()->GetOneLetterCode() == type_a &&
          CHAIN_B.GetSequence()->GetAA( seq_id_b - 1)->GetType()->GetOneLetterCode() == type_b,
          " The provided contact pair does not match the provided sequences" +
          util::Format()( CHAIN_A.GetSequence()->GetAA( seq_id_a - 1)->GetType()->GetOneLetterCode()) +
          " vs " + util::Format()( type_a) + " and " +
          util::Format()( CHAIN_B.GetSequence()->GetAA( seq_id_b - 1)->GetType()->GetOneLetterCode()) +
          " vs " + util::Format()( type_b)
        );

        // read the predictions vector
        ISTREAM >> contact_vector.First()( GetTypes().HELIX_HELIX)
                >> contact_vector.First()( GetTypes().HELIX_SHEET)
                >> contact_vector.First()( GetTypes().SHEET_HELIX)
                >> contact_vector.First()( GetTypes().STRAND_STRAND)
                >> contact_vector.First()( GetTypes().SHEET_SHEET);

        // read the is_in_contact from the file
        ISTREAM >> contact_vector.Second();

        // since default constructed conteact vectors have 0, but when reading all values are set to undefined
        // HELIX_STRAND and STRAND_HELIX contacts cannot be read and will remain undefined
        contact_vector.First()( GetTypes().HELIX_STRAND) = 0;
        contact_vector.First()( GetTypes().STRAND_HELIX) = 0;

        // Insert amino acid pair
        Insert
        (
          storage::VectorND< 2, util::SiPtr< const biol::AAData> >
          (
            CHAIN_A.GetSequence()->GetAA( seq_id_a - 1)->GetData(),
            CHAIN_B.GetSequence()->GetAA( seq_id_b - 1)->GetData()
          ),
          contact_vector
        );
      }

      // end
      return ISTREAM;
    }

    //! @brief Writes predictions for the provided chains
    //! @param OSTREAM output stream
    //! @param CHAIN_A first chain
    //! @param CHAIN_B second chain
    //! @param WRITE_ONLY_CONTACTS boolean value that determines whether non contacting residue couples are written
    //! @return ostream which was written to
    std::ostream &Map::WriteContacts
    (
      std::ostream &OSTREAM,
      const assemble::Chain &CHAIN_A,
      const assemble::Chain &CHAIN_B,
      const bool WRITE_ONLY_CONTACTS
    ) const
    {
      // iterate over first sequence
      for
      (
        biol::AASequence::const_iterator aa_itr_a( CHAIN_A.GetSequence()->Begin()),
          aa_itr_a_end( CHAIN_A.GetSequence()->End());
        aa_itr_a != aa_itr_a_end;
        ++aa_itr_a
      )
      {
        // iterate over second sequence
        for
        (
          biol::AASequence::const_iterator aa_itr_b( CHAIN_B.GetSequence()->Begin()),
            aa_itr_b_end( CHAIN_B.GetSequence()->End());
          aa_itr_b != aa_itr_b_end;
          ++aa_itr_b
        )
        {
          // get the contact vector for the two amino acids
          storage::Pair< linal::Vector< size_t>, bool> contact_vector
          (
            GetContactVector
            (
              storage::VectorND< 2, util::SiPtr< const biol::AAData> >
              (
                ( *aa_itr_a)->GetData(),
                ( *aa_itr_b)->GetData()
              )
            )
          );

          // if these residues are not in contact and WRITE_ONLY_CONTACTS flag is set skip
          if( WRITE_ONLY_CONTACTS && contact_vector.Second() == 0)
          {
            continue;
          }

          // if these residues are not found
          if( !contact_vector.First().IsDefined())
          {
            // if write_only_contacts flag is true skip
            if( WRITE_ONLY_CONTACTS)
            {
              continue;
            }
            else
            {
              contact_vector.First() = linal::Vector< size_t>( size_t( Types::s_NumberValidTypes), size_t( 0));
              contact_vector.Second() = false;
            }
          }

          // output values
          OSTREAM << util::Format().W( 4)( ( *aa_itr_a)->GetSeqID()) << ' '
                  << ( *aa_itr_a)->GetType()->GetOneLetterCode() << ' '
                  << util::Format().W( 4)( ( *aa_itr_b)->GetSeqID()) << ' '
                  << ( *aa_itr_b)->GetType()->GetOneLetterCode() << ' '
                  << util::Format()( contact_vector.First()( GetTypes().HELIX_HELIX)) << ' '
                  << util::Format()( contact_vector.First()( GetTypes().HELIX_SHEET)) << ' '
                  << util::Format()( contact_vector.First()( GetTypes().SHEET_HELIX)) << ' '
                  << util::Format()( contact_vector.First()( GetTypes().STRAND_STRAND)) << ' '
                  << util::Format()( contact_vector.First()( GetTypes().SHEET_SHEET)) << ' '
                  << util::Format()( contact_vector.Second()) << '\n';
        }
      }

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Given 2 AAs, this function returns a bool of whether these two residues
    //! @brief are in contact according to rules of the specified contact type
    //! @param AMINO_ACID_A first amino acid
    //! @param AMINO_ACID_B second amino acid
    //! @return boolean indicating whether AMINO_ACID_A and AMINO_ACID_B are in contact with any contact type
    bool Map::IsInContact( const biol::AABase &AMINO_ACID_A, const biol::AABase &AMINO_ACID_B)
    {
      for( Types::const_iterator itr( GetTypes().Begin()), itr_end( GetTypes().End()); itr != itr_end; ++itr)
      {
        bool in_contact( IsInContact( AMINO_ACID_A, AMINO_ACID_B, *itr));
        if( in_contact)
        {
          return true;
        }
      }
      return false;
    }

    //! @brief Given 2 AAs and a contact type, this function returns a bool of whether these two residues
    //! @brief are in contact according to rules of the specified contact type
    //! @param AMINO_ACID_A first amino acid
    //! @param AMINO_ACID_B second amino acid
    //! @param CONTACT_TYPE contact type
    //! @return boolean indicating whether AMINO_ACID_A and AMINO_ACID_B are in contact with CONTACT_TYPE
    bool Map::IsInContact
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B,
      const Type &CONTACT_TYPE
    )
    {
      // if it is one of the valid contacttypes
      if( CONTACT_TYPE->IsValid())
      {
        // calculate the distance between CB atoms of two residues
        const double distance( Distance( AMINO_ACID_A.GetFirstSidechainAtom(), AMINO_ACID_B.GetFirstSidechainAtom()));

        // if  distance is smaller than the supplied threshold for that given contacttype
        if( distance < CONTACT_TYPE->GetResidueDistanceCutoff())
        {
          return true;
        }
      }
      // flag is supplied, check the unknown contacts
      if( CONTACT_TYPE == GetTypes().e_Undefined)
      {
        // calculate the distance between CB atoms of two residues
        const double distance( Distance( AMINO_ACID_A.GetFirstSidechainAtom(), AMINO_ACID_B.GetFirstSidechainAtom()));

        // if distance is smaller than the supplied threshold for undefined contacts
        if( distance < Types::GetUnknownResidueDistanceCutoff())
        {
          return true;
        }
      }
      // else return false
      return false;
    }

    //! @brief checks if the two SSEs are in contact in the given contact map
    //! @param CONTACTS the contact map the defines the contacts between amino acids
    //! @param FIRST first SSE
    //! @param SECOND second SSE
    //! @return true if the SSEs are in contact
    bool Map::IsInContact( const Map &CONTACTS, const assemble::SSE &FIRST, const assemble::SSE &SECOND)
    {
      util::SiPtrVector< const biol::AABase> first_aas( FIRST.GetData()), second_aas( SECOND.GetData());
      for
      (
        util::SiPtrVector< const biol::AABase>::const_iterator
          itr_aa_a( first_aas.Begin()), itr_aa_a_end( first_aas.End());
        itr_aa_a != itr_aa_a_end; ++itr_aa_a
      )
      {
        for
        (
          util::SiPtrVector< const biol::AABase>::const_iterator
            itr_aa_b( second_aas.Begin()), itr_aa_b_end( second_aas.End());
          itr_aa_b != itr_aa_b_end; ++itr_aa_b
        )
        {
          if( CONTACTS.IsInContact( **itr_aa_a, **itr_aa_b))
          {
            return true;
          }
        }
      }
      return false;
    }

    //! @brief This function fills m_Map with real contact information ( 0|1) for each residue couple from
    //! @brief each contacttype for the given set of sses
    void Map::FillMap()
    {
      BCL_MessageDbg( "Initializing the contactmap");

      // make sure all chains have SSE information
      CheckChainsForSSEs();

      // iterate over every chain stored
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator chain_itr_a( m_Chains.Begin()),
          chain_itr_end( m_Chains.End());
        chain_itr_a != chain_itr_end;
        ++chain_itr_a
      )
      {
        FillMap( **chain_itr_a);

        // iterate over every other chain
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator chain_itr_b( chain_itr_a + 1);
          chain_itr_b != chain_itr_end;
          ++chain_itr_b
        )
        {
          FillMap( **chain_itr_a, **chain_itr_b);
        }
      }
    }

    //! @brief This function fills m_Map with real contact information ( 0|1) for each residue couple from
    //! @brief each contact type for the given chain
    //! @param CHAIN a chain
    void Map::FillMap( const assemble::Chain &CHAIN)
    {
      // iterate over every sse element( sse_itr_a)
      for
      (
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
          sse_itr_a( CHAIN.GetData().Begin()), sse_itr_end( CHAIN.GetData().End());
        sse_itr_a != sse_itr_end;
        ++sse_itr_a
      )
      {
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
        sse_itr_b( sse_itr_a);
        ++sse_itr_b;

        // iterate over every later sse element( sse_itr_b)
        for
        (
          ;
          sse_itr_b != sse_itr_end;
          ++sse_itr_b
        )
        {
          // precalculate the contact type between sses, since this will stay same for all aa in these two sses
          Type contact_type
          (
            assemble::SSEGeometryPacking( **sse_itr_a, **sse_itr_b).GetContactType()
          );

          // HELIX_STRAND and STRAND_HELIX not handled yet - set to undefined
          if( contact_type == GetTypes().HELIX_STRAND)
          {
            contact_type = GetTypes().UNDEFINED_HELIX_STRAND;
          }
          else if( contact_type == GetTypes().STRAND_HELIX)
          {
            contact_type = GetTypes().UNDEFINED_STRAND_HELIX;
          }

          BCL_MessageDbg( "Determined contact type :" + util::Format()( contact_type));

          // iterate over every residue in sse_itr_a
          for
          (
            biol::AASequence::const_iterator aa_itr_a( ( *sse_itr_a)->Begin()),
              aa_itr_end_a( ( *sse_itr_a)->End());
            aa_itr_a != aa_itr_end_a;
            ++aa_itr_a
          )
          {
            // vs every residue in sse_itr_b
            for
            (
              biol::AASequence::const_iterator aa_itr_b( ( *sse_itr_b)->Begin()),
                aa_itr_end_b( ( *sse_itr_b)->End());
              aa_itr_b != aa_itr_end_b;
              ++aa_itr_b
            )
            {
              // if lastseqid was set
              if( CheckBoundaryCondition( **aa_itr_a, CHAIN) && CheckBoundaryCondition( **aa_itr_b, CHAIN))
              {
                InsertAminoAcidPair( **aa_itr_a, **aa_itr_b, contact_type, true);
              }
            }
          }
        }
      }
    }

    //! @brief This function fills m_Map with real contact information ( 0|1) for each residue couple from
    //! @brief each contact type for the given pair of chains
    //! @param CHAIN_A first chain
    //! @param CHAIN_B second chain
    void Map::FillMap( const assemble::Chain &CHAIN_A, const assemble::Chain &CHAIN_B)
    {
      // iterate over every sse element( sse_itr_a) in first chain
      for
      (
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
          sse_itr_a( CHAIN_A.GetData().Begin()), sse_itr_end_a( CHAIN_A.GetData().End());
        sse_itr_a != sse_itr_end_a;
        ++sse_itr_a
      )
      {
        // iterate over every later sse element( sse_itr_b) in the second chain
        for
        (
          storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
            sse_itr_b( CHAIN_B.GetData().Begin()), sse_itr_end_b( CHAIN_B.GetData().End());
          sse_itr_b != sse_itr_end_b;
          ++sse_itr_b
        )
        {
          // precalculate the contactype between sses, since this will stay same for all aa in these two sses
          Type contact_type
          (
            assemble::SSEGeometryPacking( **sse_itr_a, **sse_itr_b).GetContactType()
          );

          // HELIX_STRAND and STRAND_HELIX not handled yet - set to undefined
          if( contact_type == GetTypes().HELIX_STRAND || contact_type == GetTypes().STRAND_HELIX)
          {
            contact_type = GetTypes().e_Undefined;
          }

          BCL_MessageDbg( "Determined contact type :" + util::Format()( contact_type));

          // iterate over every residue in sse_itr_a
          for
          (
            biol::AASequence::const_iterator aa_itr_a( ( *sse_itr_a)->Begin()),
              aa_itr_end_a( ( *sse_itr_a)->End());
            aa_itr_a != aa_itr_end_a;
            ++aa_itr_a
          )
          {
            // vs every residue in sse_itr_b
            for
            (
              biol::AASequence::const_iterator aa_itr_b( ( *sse_itr_b)->Begin()),
                aa_itr_end_b( ( *sse_itr_b)->End());
              aa_itr_b != aa_itr_end_b;
              ++aa_itr_b
            )
            {
              // if lastseqid was set
              if( CheckBoundaryCondition( **aa_itr_a, CHAIN_A) && CheckBoundaryCondition( **aa_itr_b, CHAIN_B))
              {
                InsertAminoAcidPair( **aa_itr_a, **aa_itr_b, contact_type, true);
              }
            }
          }
        }
      }
    }

    //! @brief checks the boundary conditions of a residue to determine whether it should be included in the map
    //! @param AMINO_ACID amino acid for which the conditions are going to be checked
    //! @param CHAIN
    bool Map::CheckBoundaryCondition( const biol::AABase &AMINO_ACID, const assemble::Chain &CHAIN) const
    {
      return
      (
        AMINO_ACID.GetSeqID() > m_Boundary
      )
      &&
      (
        AMINO_ACID.GetSeqID() < CHAIN.GetSequence()->GetLastAA()->GetSeqID() - m_Boundary
      );
    }

    //! @brief insert the provided amino acid pair and contact type and optional reverse contact type,
    //! @brief inserts the related data into hashmap and residue pair list
    //! @param AMINO_ACID_A first amino acid
    //! @param AMINO_ACID_B second amino acid
    //! @param CONTACT_TYPE contact type
    //! @param INSERT_REVERSE also insert reverse of CONTACT_TYPE
    void Map::InsertAminoAcidPair
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B,
      const Type &CONTACT_TYPE,
      const bool INSERT_REVERSE
    )
    {
      //! create the vectors for storage
      storage::Pair< linal::Vector< size_t>, bool> contact_vector
      (
        linal::Vector< size_t>( Types::s_NumberValidTypes, size_t( 0)),
        false
      );

      // determine if these two residues are in contact
      const bool is_in_contact( IsInContact( AMINO_ACID_A, AMINO_ACID_B, CONTACT_TYPE));

      // if in contact, update the vector
      if( is_in_contact && CONTACT_TYPE->IsValid())
      {
        contact_vector.First()( CONTACT_TYPE) = true;
      }

      // set the boolean value accordingly
      contact_vector.Second() = is_in_contact;

      // create simple pointers to amino acid data
      util::SiPtr< const biol::AAData> sp_a( AMINO_ACID_A.GetData());
      util::SiPtr< const biol::AAData> sp_b( AMINO_ACID_B.GetData());

      // insert h,i
      Insert
      (
        storage::VectorND< 2, util::SiPtr< const biol::AAData> >( sp_a, sp_b),
        contact_vector
      );

      // if reverse flag is set
      if( INSERT_REVERSE)
      {
        storage::Pair< linal::Vector< size_t>, bool> contact_vector_swapped( contact_vector);
        linal::Vector< size_t> &contact_vector( contact_vector_swapped.First());
        std::swap( contact_vector( GetTypes().HELIX_SHEET), contact_vector( GetTypes().SHEET_HELIX));

        // insert i,h
        Insert
        (
          storage::VectorND< 2, util::SiPtr< const biol::AAData> >( sp_b, sp_a),
          contact_vector_swapped
        );
      }
    }

    //! @brief this functions checks that the chains stored also have SSE information
    //! @brief so that the contact map can be generated correctly
    void Map::CheckChainsForSSEs() const
    {
      // iterate over every chain stored
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator chain_itr( m_Chains.Begin()),
          chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        //check that there is at least 1 sse stored with this chain
         BCL_Assert
         (
           !( *chain_itr)->GetData().IsEmpty(),
           "The following chains does not have any SSE information: " + util::Format()( ( *chain_itr)->GetChainID())
         );
      }
    }

  } // namespace contact
} // namespace bcl
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
#include "contact/bcl_contact_order.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container_generator_protein_model.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_function_cached.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {
    math::FunctionCached< assemble::ProteinModel, double>
      Order::s_AbsoluteInstanceCached
      (
        util::ShPtr< Order>( new Order( e_Absolute, "raw")),
        &assemble::ProteinModel::GetDestructorSignal,
        &assemble::ProteinModel::GetChangeSignal
      );

    //! @brief NormalizationType as string
    //! @param NORMALIZATION_TYPE the NormalizationType
    //! @return the string for the NormalizationType
    const std::string &Order::GetNormalizationTypeDescriptor( const NormalizationType &NORMALIZATION_TYPE)
    {
      static const std::string s_descriptors[] =
      {
        "Absolute",
        "RelativeAAUsed",
        "RelativeSeqeunceLength",
        "RelativeSqrsAAUsed",
        "RelativeSqrSeqeunceLength",
        GetStaticClassName< NormalizationType>()
      };

      return s_descriptors[ NORMALIZATION_TYPE];
    }

  //////////
  // data //
  //////////

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &Order::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "co");

      // end
      return s_default_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Order::Order() :
      m_Normalization( e_Absolute),
      m_Scheme( GetDefaultScheme()),
      m_CacheNeighborLists( false),
      m_NeighborGenerator
      (
        assemble::AANeighborListContainerGeneratorProteinModel::AANeighborListGenerator
        (
          g_ContactCbDistanceCutoff, s_ContactOrderSequencesDistanceCutoff, false, false
        )
      )
    {
    }

    //! @brief constructor from normalization type, scheme string and whether to use cache generating neighbor lists
    //! @brief NORMALIZE type of normalization
    //! @brief SCHEME scheme to report
    //! @brief CACHE whether to use cache generating neighbor lists
    Order::Order
    (
      const NormalizationType NORMALIZE,
      const std::string &SCHEME,
      const bool CACHE
    ) :
      m_Normalization( NORMALIZE),
      m_Scheme( SCHEME),
      m_CacheNeighborLists( CACHE),
      m_NeighborGenerator
      (
        assemble::AANeighborListContainerGeneratorProteinModel::AANeighborListGenerator
        (
          g_ContactCbDistanceCutoff, s_ContactOrderSequencesDistanceCutoff, false, CACHE
        )
      )
    {
    }

    //! @brief copy constructor
    Order *Order::Clone() const
    {
      return new Order( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Order::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief contact order for aa sequence
    //! @param AA_SEQUENCE sequence of interest
    //! @return contact order of sequence
    double Order::ContactOrder( const biol::AASequence &AA_SEQUENCE) const
    {
      // all amino acids
      const util::SiPtrVector< const biol::AABase> amino_acids( AA_SEQUENCE.GetMembers());

      // return normalized contact order
      return Normalize
      (
        ContactOrder
        (
          assemble::AANeighborListContainer
          (
            amino_acids, g_ContactCbDistanceCutoff, s_ContactOrderSequencesDistanceCutoff, false
          )
        ),
        amino_acids.GetSize(),
        amino_acids.GetSize()
      );
    }

    //! @brief contact order for the given chain
    //! @param CHAIN chain of interest
    //! @return contact order for given chain
    double Order::ContactOrder( const assemble::Chain &CHAIN) const
    {
      // all amino acids
      const util::SiPtrVector< const biol::AABase> amino_acids( CHAIN.GetAminoAcids());

      // return normalized contact order
      return Normalize
      (
        ContactOrder
        (
          assemble::AANeighborListContainer
          (
            amino_acids, g_ContactCbDistanceCutoff, s_ContactOrderSequencesDistanceCutoff, false
          )
        ),
        amino_acids.GetSize(),
        CHAIN.GetSequence()->GetSize()
      );
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that calculates the contact order for a given protein model
    //! @param PROTEIN_MODEL protein model
    //! @return contact order for the given protein model
    double Order::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // currently, only one chain is defined
      if( PROTEIN_MODEL.GetChains().GetSize() != 1)
      {
        BCL_MessageVrb( "only single chain protein contact orders can be calculated!");
        return util::GetUndefined< double>();
      }
      // test whether this is the static instance. If so, it is necessary to compute the contact order
      if( this == &*s_AbsoluteInstanceCached.GetFunction())
      {
        return ContactOrder( m_NeighborGenerator->operator ()( PROTEIN_MODEL));
      }

      // return normalized contact order
      return Normalize
      (
        s_AbsoluteInstanceCached( PROTEIN_MODEL),
        PROTEIN_MODEL.GetChains().FirstElement()->GetNumberAAs(),
        PROTEIN_MODEL.GetChains().FirstElement()->GetSequence()->GetSize()
      );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer Order::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        "Computes the contact order of a protein model allowing for various means of normalization"
      );
      serializer.AddInitializer
      (
        "normalization",
        "method of normalization",
        io::Serialization::GetAgent( &m_Normalization),
        GetNormalizationTypeDescriptor( e_Absolute)
      );
      serializer.AddInitializer
      (
        "cache neighbor lists",
        "whether to cache neighbor lists. This may speed up calculations if many methods make use of neighbor lists on "
        "the same protein model. If contact order is the only neighbor list being used, it may slow the computation "
        "down by ~10%",
        io::Serialization::GetAgent( &m_CacheNeighborLists),
        "False"
      );
      return serializer;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool Order::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      m_NeighborGenerator =
        assemble::AANeighborListContainerGeneratorProteinModel::AANeighborListGenerator
        (
          g_ContactCbDistanceCutoff, s_ContactOrderSequencesDistanceCutoff, false, m_CacheNeighborLists
        );
      return true;
    }

  ////////////
  // helper //
  ////////////

    //! @brief normalize the contact order
    //! @param CONTACT_ORDER the unnormalized contact order
    //! @param NR_AAS_USED nr of amino acids used
    //! @param SEQUENCE_LENGTH the length of the amino acid sequence
    //! @return normalized contact order, depending of m_Normalization
    double Order::Normalize( const double CONTACT_ORDER, const size_t NR_AAS_USED, const size_t SEQUENCE_LENGTH) const
    {
      switch( m_Normalization)
      {
        case e_Absolute                 : return CONTACT_ORDER;
        case e_RelativeAAsUsed          : return NR_AAS_USED == 0 ? 0 : CONTACT_ORDER / double( NR_AAS_USED);
        case e_RelativeSequenceLength   : return SEQUENCE_LENGTH == 0 ? 0 : CONTACT_ORDER / double( SEQUENCE_LENGTH);
        case e_RelativeSqrAAsUsed       : return NR_AAS_USED == 0 ? 0 : math::Sqr( CONTACT_ORDER) / double( NR_AAS_USED);
        case e_RelativeSqrSequenceLength: return SEQUENCE_LENGTH == 0 ? 0 : math::Sqr( CONTACT_ORDER) / double( SEQUENCE_LENGTH);
        default                         : return util::GetUndefined< double>();
      }

      // should never reach this
      return util::GetUndefined< double>();
    }

    //! @brief calculate contact order form NeighborListContainer
    //! @param NEIGHBOR_LIST_CONTAINER container of neighbors
    //! @return the contact order
    double Order::ContactOrder( const assemble::AANeighborListContainer &NEIGHBOR_LIST_CONTAINER)
    {
      size_t seq_distance_sum( 0);
      size_t contacts( 0);

      // iterate over all aa neighbot list
      for
      (
        assemble::AANeighborListContainer::const_iterator
          itr( NEIGHBOR_LIST_CONTAINER.Begin()), itr_end( NEIGHBOR_LIST_CONTAINER.End());
        itr != itr_end;
        ++itr
      )
      {
        // sequence id of the center amino acid
        const int center_seq_id( itr->second.GetCenterAminoAcid()->GetSeqID());

        // iterate over neighbors
        for
        (
          assemble::AANeighborList::const_iterator neigh_itr( itr->second.Begin()), neigh_itr_end( itr->second.End());
          neigh_itr != neigh_itr_end;
          ++neigh_itr
        )
        {
          // calculate sequence distance
          const int seq_distance( neigh_itr->First()->GetSeqID() - center_seq_id);

          // assure that only one of the two pairs (stored redundant in AANeighborList) are used and test the distance
          if
          (
               seq_distance < int( s_ContactOrderSequencesDistanceCutoff)
            || neigh_itr->Second() >= g_ContactCbDistanceCutoff
          )
          {
            continue;
          }

          seq_distance_sum += size_t( seq_distance);
          contacts++;
        }
      }

      // if no contacts are found
      if( contacts == 0)
      {
        return double( 0.0);
      }

      // otherwise normalize
      return double( seq_distance_sum) / double( contacts);
    }

  } // namespace contact
} // namespace bcl
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
#include "contact/bcl_contact_prediction_map.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "contact/bcl_contact_ann.h"
#include "linal/bcl_linal_matrix.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {

  //////////
  // data //
  //////////

    //! How many residues to exclude from begin and end
    static const size_t s_NumberSkippedResiduesInEnds = size_t( 4);

    //! The minimum number of residues a sequence should have contact prediction
    static const size_t s_MinimumNumberResidues = size_t( 20);

    //! minimal sequence distance for contacts
    static const size_t s_MinSequenceDistance = size_t( 4);

    //! The number of values used for describing each residue for prediction (7 properties + 3 JUFO + 20 Blast)
    static const size_t s_NumberOfValuesPerResidue = size_t( 30);

    //! identifier string to be used for prediction maps
    static const std::string s_Identifier = "CHAINS";

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PredictionMap::PredictionMap() :
      storage::ObjectNDHashMap< 2, biol::AAData, storage::Pair< linal::Vector< double>, double> >(),
      m_Sequences()
    {
    }

    //! @brief construct prediction map ( from ANNs) from one chain( intra-sequence contacts)
    //! @param CHAIN chain
    PredictionMap::PredictionMap( const assemble::Chain &CHAIN) :
      storage::ObjectNDHashMap< 2, biol::AAData, storage::Pair< linal::Vector< double>, double> >(),
      m_Sequences( 1, CHAIN.GetSequence())
    {
      // make sure the sequence matches the size
      BCL_Assert
      (
        CHAIN.GetSequence()->GetSize() >= s_MinimumNumberResidues,
        "The given sequence needs to be at least " + util::Format()( s_MinimumNumberResidues) +
        " residues long, but it has only " + util::Format()( CHAIN.GetSequence()->GetSize()) + " residues"
      );
      // fill the map
      FillMap();
    }

    //! @brief construct predion map from ProteinModel
    //! @param PROTEIN_MODEL Protein Model for which prediction map is going to be constructed
    PredictionMap::PredictionMap( const assemble::ProteinModel &PROTEIN_MODEL) :
      storage::ObjectNDHashMap< 2, biol::AAData, storage::Pair< linal::Vector< double>, double> >(),
      m_Sequences()
    {

      // iterate over the chains in the model
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
          chain_itr( PROTEIN_MODEL.GetChains().Begin()), chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end; ++chain_itr
      )
      {
        // get the sequence
        const util::ShPtr< biol::AASequence> &sp_sequence( ( *chain_itr)->GetSequence());

        // make sure the sequence matches the size
        BCL_Assert
        (
          sp_sequence->GetSize() >= s_MinimumNumberResidues,
          "The given sequence needs to be at least " + util::Format()( s_MinimumNumberResidues) +
          " residues long, but it has only " + util::Format()( sp_sequence->GetSize()) + " residues"
        );

        // insert into sequences.
        m_Sequences.PushBack( sp_sequence);
      }

      // fill the map
      FillMap();
    }

    //! @brief virtual copy constructor
    PredictionMap *PredictionMap::Clone() const
    {
      return new PredictionMap( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PredictionMap::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief iterates over the stored sequences and finds the one with matching CHAIN_ID
    //! @param CHAIN_ID chain id of the sequence that is beings searched
    //! @return AASequence with the searched chain id
    const util::ShPtr< biol::AASequence> &PredictionMap::GetSequence( const char CHAIN_ID) const
    {
      // iterate over every sequence stored
      for
      (
        util::ShPtrVector< biol::AASequence>::const_iterator sequence_itr( m_Sequences.Begin()),
          sequence_itr_end( m_Sequences.End());
        sequence_itr != sequence_itr_end;
        ++sequence_itr
      )
      {
        // if the sequence id matches return it
        if( ( *sequence_itr)->GetChainID() == CHAIN_ID)
        {
          return *sequence_itr;
        }
      }

      // else Exit
      BCL_Exit( "No sequence with the provided chain id \'" + util::Format()( CHAIN_ID) + "\' is not stored!!", -1);
      return m_Sequences.FirstElement();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Returns predictions for the provided AA Data pointers
    //! @param AA_DATA_POINTERS pair of pointers to AAData of amino acids of interest
    //! @return the pair of vector of predictions and the merged prediction
    const storage::Pair< linal::Vector< double>, double> &PredictionMap::GetPredictions
    (
      const storage::VectorND< 2, util::SiPtr< const biol::AAData> > &AA_DATA_POINTERS
    ) const
    {
      // initialize undefined predictions vector
      static const storage::Pair< linal::Vector< double>, double> s_undefined_predictions
      (
        linal::Vector< double>( Types::s_NumberValidTypes, util::GetUndefined< double>()),
        util::GetUndefined< double>()
      );

      // search in the hash map
      storage::HashMap< size_t, storage::Pair< linal::Vector< double>, double> >::const_iterator itr
      (
        Find( AA_DATA_POINTERS)
      );

      // if not found return undefined
      if( itr == End())
      {
        return s_undefined_predictions;
      }

      // end
      return itr->second;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PredictionMap::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &PredictionMap::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

    //! @brief helper function to read predictions from a file
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PredictionMap::ReadPredictionMap
    (
      std::istream &ISTREAM
    )
    {
      // header that is going to be searched
      std::string header;

      // while reading header and not end of file yet
      while( ISTREAM >> header && !ISTREAM.eof())
      {
        // if header matches the identifier
        if( header == s_Identifier)
        {

          // read the two chain ids
          char chain_id_a, chain_id_b;
          ISTREAM >> chain_id_a >> chain_id_b;

          // check the chain_id_a exists
          BCL_Assert
          (
            GetSequence( chain_id_a).IsDefined(),
            "No sequence is stored in the map with given chain id: " + std::string( 1, chain_id_a)
          );

          // check the chain_id_b exists
          BCL_Assert
          (
            GetSequence( chain_id_b).IsDefined(),
            "No sequence is stored in the map with given chain id: " + std::string( 1, chain_id_b)
          );

          // read the predictions for two sequences with the provided chain ids
          ReadPredictions
          (
            ISTREAM,
            *GetSequence( chain_id_a),
            *GetSequence( chain_id_b)
          );
        }
      }

      // end
      return ISTREAM;
    }

    //! @brief helper function to write predictions to a file
    //! @param OSTREAM output stream to write to
    //! @param MERGED_THRESHOLD if merged predictions are below this threshold do not print them out
    //! @return output stream which was written to
    std::ostream &PredictionMap::WritePredictionMap
    (
      std::ostream &OSTREAM,
      const double MERGED_THRESHOLD
    ) const
    {

      // iterate over every sequence
      for
      (
        util::ShPtrVector< biol::AASequence>::const_iterator seq_itr_a( m_Sequences.Begin()),
          seq_itr_end( m_Sequences.End());
        seq_itr_a != seq_itr_end;
        ++seq_itr_a
      )
      {
        // versus every other sequence
        for
        (
          util::ShPtrVector< biol::AASequence>::const_iterator seq_itr_b( m_Sequences.Begin());
          seq_itr_b != seq_itr_end;
          ++seq_itr_b
        )
        {
          // output the header
          OSTREAM << s_Identifier << ' '
                  << util::Format()( (*seq_itr_a)->GetChainID()) << ' '
                  << util::Format()( (*seq_itr_b)->GetChainID()) << '\n';

          // write predictions
          WritePredictions( OSTREAM, **seq_itr_a, **seq_itr_b, MERGED_THRESHOLD);
        }
      }
      return OSTREAM;
    }

    //! @brief Reads predictions for the provided sequences
    //! @param ISTREAM input stream
    //! @param SEQUENCE_A first sequence
    //! @param SEQUENCE_B second sequence
    //! @return istream which was read from
    std::istream &PredictionMap::ReadPredictions
    (
      std::istream &ISTREAM,
      const biol::AASequence &SEQUENCE_A,
      const biol::AASequence &SEQUENCE_B
    )
    {

      // initialize sequence id and AAType pairs and the vector hold predictions
      int seq_id_a, seq_id_b;
      char type_a, type_b;
      linal::Vector< double> predictions( Types::s_NumberValidTypes, util::GetUndefined< double>());

      // while it is possible to read
      while( ISTREAM >> seq_id_a >> type_a >> seq_id_b >> type_b && !ISTREAM.eof())
      {
        // make the sure sequences match
        BCL_Assert
        (
          SEQUENCE_A.GetAA( seq_id_a - 1)->GetType()->GetOneLetterCode() == type_a &&
          SEQUENCE_B.GetAA( seq_id_b - 1)->GetType()->GetOneLetterCode() == type_b,
          " The provided contact pair does not match the provided sequences" +
          util::Format()( SEQUENCE_A.GetAA( seq_id_a - 1)->GetType()->GetOneLetterCode()) +
          " vs " + util::Format()( type_a) + " and " +
          util::Format()( SEQUENCE_B.GetAA( seq_id_b - 1)->GetType()->GetOneLetterCode()) +
          " vs " + util::Format()( type_b)
        );

        // read the predictions vector
        ISTREAM >> predictions( GetTypes().HELIX_HELIX)
                >> predictions( GetTypes().HELIX_SHEET)
                >> predictions( GetTypes().SHEET_HELIX)
                >> predictions( GetTypes().STRAND_STRAND)
                >> predictions( GetTypes().SHEET_SHEET);

        // read the merged prediction from the file
        double merged_prediction;
        ISTREAM >> merged_prediction;

        // set prediction for unhandled HELIX_STRAND and STRAND_HELIX to 0
        predictions( GetTypes().HELIX_STRAND) = 0.0;
        predictions( GetTypes().STRAND_HELIX) = 0.0;

        // insert the predictions into the predictionmap
        InsertPredictions( *SEQUENCE_A.GetAA( seq_id_a - 1), *SEQUENCE_B.GetAA( seq_id_b - 1), predictions);
      }

      // end
      return ISTREAM;
    }

    //! @brief Writes predictions for the provided sequences
    //! @param OSTREAM output stream
    //! @param SEQUENCE_A first sequence
    //! @param SEQUENCE_B second sequence
    //! @param MERGED_THRESHOLD if merged predictions are below this threshold do not print them out
    //! @return ostream which was written to
    std::ostream &PredictionMap::WritePredictions
    (
      std::ostream &OSTREAM,
      const biol::AASequence &SEQUENCE_A,
      const biol::AASequence &SEQUENCE_B,
      const double MERGED_THRESHOLD
    ) const
    {

      // iterate over first sequence
      for
      (
        biol::AASequence::const_iterator aa_itr_a( SEQUENCE_A.Begin()),
          aa_itr_a_end( SEQUENCE_A.End());
        aa_itr_a != aa_itr_a_end;
        ++aa_itr_a
      )
      {
        // iterate over second sequence
        for
        (
          biol::AASequence::const_iterator aa_itr_b( SEQUENCE_B.Begin()),
            aa_itr_b_end( SEQUENCE_B.End());
          aa_itr_b != aa_itr_b_end;
          ++aa_itr_b
        )
        {
          // get the predictions for the two amino acids
          const storage::Pair< linal::Vector< double>, double> predictions
          (
            GetPredictions
            (
              storage::VectorND< 2, util::SiPtr< const biol::AAData> >
              (
                ( *aa_itr_a)->GetData(),
                ( *aa_itr_b)->GetData()
              )
            )
          );

          // if predictions are not defined, meaning they were not calculated, skip this residue pair
          if( !predictions.First().IsDefined())
          {
            continue;
          }
          // if the merged prediction is below the threshold, skip this residue pair
          if( predictions.Second() < MERGED_THRESHOLD)
          {
            continue;
          }

          static const util::Format s_format( util::Format().W( 5).FFP( 3));
          // output values
          OSTREAM << util::Format().W( 4)( ( *aa_itr_a)->GetSeqID()) << ' '
                  << ( *aa_itr_a)->GetType()->GetOneLetterCode() << ' '
                  << util::Format().W( 4)( ( *aa_itr_b)->GetSeqID()) << ' '
                  << ( *aa_itr_b)->GetType()->GetOneLetterCode() << ' '
                  << s_format( predictions.First()( GetTypes().HELIX_HELIX)) << ' '
                  << s_format( predictions.First()( GetTypes().HELIX_SHEET)) << ' '
                  << s_format( predictions.First()( GetTypes().SHEET_HELIX)) << ' '
                  << s_format( predictions.First()( GetTypes().STRAND_STRAND)) << ' '
                  << s_format( predictions.First()( GetTypes().SHEET_SHEET)) << ' '
                  << s_format( predictions.Second()) << '\n';
        }
      }

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief This function generates the description for a single amino acid to be used for Contact ANNs
    //! @param AMINO_ACID amino acid of interest for which description will be generated
    //! @return the descriptor vector composed of 7(properties)+3(jufo)+20(blast)=30 values for a single amino acid
    linal::Vector< double> PredictionMap::GenerateDescription( const biol::AABase &AMINO_ACID)
    {
      // initialize an empty description vector
      linal::Vector< double> description_vector( s_NumberOfValuesPerResidue);

      //insert 7 amino acid properties
      description_vector.ReplaceElements( 0, AMINO_ACID.GetType()->GetPropertiesForANN());

      // insert jufo
      description_vector.ReplaceElements
      (
        biol::AATypeData::s_NumberPropertyTypesForANN,
        AMINO_ACID.GetSSPrediction( sspred::GetMethods().e_JUFO)->GetThreeStatePrediction()
      );

      // insert blast profile
      description_vector.ReplaceElements
      (
        biol::AATypeData::s_NumberPropertyTypesForANN + 3,
        AMINO_ACID.GetBlastProfile().GetProfile()
      );

      // end
      return description_vector;
    }

    //! @brief This function generates the descriptions for all amino acids in the provided SEQUENCE
    //! @brief to used for contact ANNS. It aims to prevent repetitive generation of descriptions
    //! @param SEQUENCE AASequence of interest for which descriptions will be generated
    //! @return the descriptor matrix which include descriptions for each individiual amino acid
    linal::Matrix< double> PredictionMap::GenerateDescription( const biol::AASequence &SEQUENCE)
    {
      // initialize empty matrix
      linal::Matrix< double> description_matrix( SEQUENCE.GetSize(), s_NumberOfValuesPerResidue);

      // initialize row number
      size_t row( 0), row_end( description_matrix.GetNumberRows());

      // iterate over every amino acid in the sequence
      for
      (
        biol::AASequence::const_iterator aa_itr( SEQUENCE.Begin()),
          aa_itr_end( SEQUENCE.End());
        aa_itr != aa_itr_end && row < row_end;
        ++aa_itr, ++row
      )
      {
        // insert description for amino acid behind aa_itr
        description_matrix.ReplaceRow( row, GenerateDescription( **aa_itr));
      }

      // end
      return description_matrix;
    }

    //! @brief Merges prediction from 5 contact types into one
    //! @param AMINO_ACID_A first amino acid
    //! @param AMINO_ACID_B second amino acid
    //! @param PREDICTIONS vector of predictions for 5 contact types
    //! @return merged value of 5 predictions
    double PredictionMap::MergePredictions
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B,
      const linal::Vector< double> &PREDICTIONS
    ) const
    {

      // store sspredictions for A and B
      linal::Vector3D sspredictions_a
      (
        AMINO_ACID_A.GetSSPrediction( sspred::GetMethods().e_JUFO)->GetThreeStatePrediction()
      );
      linal::Vector3D sspredictions_b
      (
        AMINO_ACID_B.GetSSPrediction( sspred::GetMethods().e_JUFO)->GetThreeStatePrediction()
      );

      // calculate the merged value
      const double merged_value
      (
        (
          PREDICTIONS( GetTypes().HELIX_HELIX) *
          sspredictions_a( biol::GetSSTypes().HELIX) *
          sspredictions_b( biol::GetSSTypes().HELIX)
        ) +
        (
          PREDICTIONS( GetTypes().HELIX_SHEET) *
          sspredictions_a( biol::GetSSTypes().HELIX) *
          sspredictions_b( biol::GetSSTypes().STRAND)
        ) +
        (
          PREDICTIONS( GetTypes().SHEET_HELIX) *
          sspredictions_a( biol::GetSSTypes().STRAND) *
          sspredictions_b( biol::GetSSTypes().HELIX)
        ) +
        (
          PREDICTIONS( GetTypes().STRAND_STRAND) *
          sspredictions_a( biol::GetSSTypes().STRAND) *
          sspredictions_b( biol::GetSSTypes().STRAND)
        ) +
        (
          PREDICTIONS( GetTypes().SHEET_SHEET) *
          sspredictions_a( biol::GetSSTypes().STRAND) *
          sspredictions_b( biol::GetSSTypes().STRAND)
        )
      );

      // return minimum of the merged value or 1.0
      return ( merged_value > 1.0 ? 1.0 : merged_value);
    }

    //! @brief for the given aa data pair, insert the predictions into the ObjectNDHashmap
    //! @param AMINO_ACID_A first amino acid
    //! @param AMINO_ACID_B second amino acid
    //! @param PREDICTIONS vector of predictions for 5 contact types
    void PredictionMap::InsertPredictions
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B,
      const linal::Vector< double> &PREDICTIONS
    )
    {
      linal::Vector< double> this_prediction( PREDICTIONS);

      // calculate the merged predictions from predictions
      const double merged_prediction( MergePredictions( AMINO_ACID_A, AMINO_ACID_B, PREDICTIONS));

      // insert h,i
      Insert
      (
        storage::VectorND< 2, util::SiPtr< const biol::AAData> >
        (
          util::SiPtr< const biol::AAData>( AMINO_ACID_A.GetData()),
          util::SiPtr< const biol::AAData>( AMINO_ACID_B.GetData())
        ),
        storage::Pair< linal::Vector< double>, double>( this_prediction, merged_prediction)
      );

    }

    //! @brief reverse the predictions so HELIX_SHEET and SHEET_HELIX values are switched
    //! @param PREDICTIONS vector of predictions which is going to be reversed
    //! @return vector where the predictions so HELIX_SHEET and SHEET_HELIX values are switched
    linal::Vector< double> PredictionMap::ReversePredictions( const linal::Vector< double> &PREDICTIONS)
    {
      // reverse the predictions
      linal::Vector< double> predictions_swapped( PREDICTIONS);
      std::swap( predictions_swapped( GetTypes().HELIX_SHEET), predictions_swapped( GetTypes().SHEET_HELIX));

      // end
      return predictions_swapped;
    }

    //! @brief Given a specific contact type and two sequence windows, this functions puts together the input
    //! @brief data and calls the corresponding ANN and returns the single prediction
    //! @param CONTACT_TYPE contact type two amino acids are in
    //! @param DESCRIPTION_A description for the first amino acid window
    //! @param DESCRIPTION_B description for the second amino acid window
    //! @param POSITIONS position descriptors for the pair of amino acids
    //! @return prediction from ANN specific to the provided contact type for the provided descriptions and positions
    double PredictionMap::PredictSingleContact
    (
      const Type &CONTACT_TYPE,
      const linal::Vector< double> &DESCRIPTION_A,
      const linal::Vector< double> &DESCRIPTION_B,
      const linal::Vector3D &POSITIONS
    )
    {
      //instantiate a vector for the input data
      linal::Vector< double> input_data( DESCRIPTION_A.GetSize() + DESCRIPTION_B.GetSize() + POSITIONS.GetSize());

      if( CONTACT_TYPE == GetTypes().SHEET_HELIX)
      {
        //fill the vector with the 3 inputdatas
        input_data.ReplaceElements( 0, DESCRIPTION_B);
        input_data.ReplaceElements( DESCRIPTION_B.GetSize(), DESCRIPTION_A);
        input_data.ReplaceElements( DESCRIPTION_A.GetSize() + DESCRIPTION_B.GetSize(), POSITIONS);
      }
      else
      {
        //fill the vector with the 3 inputdatas
        input_data.ReplaceElements( 0, DESCRIPTION_A);
        input_data.ReplaceElements( DESCRIPTION_A.GetSize(), DESCRIPTION_B);
        input_data.ReplaceElements( DESCRIPTION_A.GetSize() + DESCRIPTION_B.GetSize(), POSITIONS);
      }

      //switch over the contact types and call the according neural network for the prediction
      double value( util::GetUndefined< double>());

      if( CONTACT_TYPE == GetTypes().HELIX_HELIX)
      {
        value = ANN_CONTACT_HELIX_HELIX( input_data);
      }
      else if( CONTACT_TYPE == GetTypes().HELIX_SHEET)
      {
        value = ANN_CONTACT_HELIX_SHEET( input_data);
      }
      else if( CONTACT_TYPE == GetTypes().SHEET_HELIX)
      {
        value = ANN_CONTACT_HELIX_SHEET( input_data);
      }
      else if( CONTACT_TYPE == GetTypes().STRAND_STRAND)
      {
        value = ANN_CONTACT_STRAND_STRAND( input_data);
      }
      else if( CONTACT_TYPE == GetTypes().SHEET_SHEET)
      {
        value = ANN_CONTACT_SHEET_SHEET( input_data);
      }

      // no networks available for that currently
      else if( CONTACT_TYPE == GetTypes().HELIX_STRAND)
      {
        value = double( 0);
      }
      else if( CONTACT_TYPE == GetTypes().STRAND_HELIX)
      {
        value = double( 0);
      }

      else
      {
        return util::GetUndefined< double>();
      }

      if( util::IsDefined( value))
      {
        value = std::min( 1.0, std::max( 0.0, value));
      }
      return value;
    }

    //! @brief This function fills m_Map with prediction for each residue couple from each contact type ANN for
    //! @brief sequences stored in m_Sequences
    void PredictionMap::FillMap()
    {
      // check jufo and blast for profiles are stored for the provided sequences
      CheckBlastAndJufoForSequences();

      // iterate over every sequence stored
      for
      (
        util::ShPtrVector< biol::AASequence>::const_iterator seq_itr_a( m_Sequences.Begin()),
          seq_itr_end( m_Sequences.End());
        seq_itr_a != seq_itr_end;
        ++seq_itr_a
      )
      {
        FillMap( **seq_itr_a);

        // iterate over every other sequence
        for
        (
          util::ShPtrVector< biol::AASequence>::const_iterator seq_itr_b( seq_itr_a + 1);
          seq_itr_b != seq_itr_end;
          ++seq_itr_b
        )
        {
          FillMap( **seq_itr_a, **seq_itr_b);
        }
      }
    }

    //! @brief This function fills m_Map with prediction for each residue couple from each contacttype ANN for
    //! @brief two different sequences
    //! @param SEQUENCE AASequence to be used for filling the map with predictions
    void PredictionMap::FillMap( const biol::AASequence &SEQUENCE)
    {
      // fill up the matrix for data for every single residue
      linal::Matrix< double> data( GenerateDescription( SEQUENCE));

      // temporary vector for storing 5 predictions between 2 AA
      linal::Vector< double> predictions( Types::s_NumberValidTypes);

      // iterate over every residue h
      for
      (
        size_t aa_ctr_a( s_NumberSkippedResiduesInEnds);
        aa_ctr_a < SEQUENCE.GetSize() - s_NumberSkippedResiduesInEnds - s_MinSequenceDistance;
        ++aa_ctr_a
      )
      {
        // and every other following residue i
        for
        (
          size_t aa_ctr_b( aa_ctr_a + s_MinSequenceDistance);
          aa_ctr_b < SEQUENCE.GetSize() - s_NumberSkippedResiduesInEnds;
          ++aa_ctr_b
        )
        {
          // for every contact type form inputs and collect network results
          for( size_t type( GetTypes().HELIX_HELIX); type < Types::s_NumberValidTypes; ++type)
          {
            Type contact_type( type);

            predictions( contact_type) =
              PredictSingleContact
              (
                contact_type,
                linal::Vector< double>
                (
                  s_NumberOfValuesPerResidue * ( contact_type->GetWindowLengths().First()),
                  data[ aa_ctr_a - contact_type->GetWindowRadii().First()]
                ),
                linal::Vector< double>
                (
                  s_NumberOfValuesPerResidue * ( contact_type->GetWindowLengths().Second()),
                  data[ aa_ctr_b - contact_type->GetWindowRadii().Second()]
                ),
                linal::Vector3D( aa_ctr_a, aa_ctr_b - aa_ctr_a, SEQUENCE.GetSize() - aa_ctr_b)
              );
          }

          // insert the predictions for the selected residue pair
          InsertPredictions
          (
            *SEQUENCE.GetAA( aa_ctr_a), *SEQUENCE.GetAA( aa_ctr_b), predictions
          );

          // insert the predictions for the reversed residue pair
          InsertPredictions
          (
            *SEQUENCE.GetAA( aa_ctr_b), *SEQUENCE.GetAA( aa_ctr_a), ReversePredictions( predictions)
          );

        } // loop over i
      } // loop over h
    }

    //! @brief This function fills m_Map with prediction for each residue couple from each contacttype ANN for
    //! @brief two different sequences
    //! @param SEQUENCE_A First AASequence to be used for filling the map with predictions
    //! @param SEQUENCE_B Second AASequence to be used for filling the map with predictions
    void PredictionMap::FillMap
    (
      const biol::AASequence &SEQUENCE_A,
      const biol::AASequence &SEQUENCE_B
    )
    {
      // fill up the matrix for data for every single residue
      const linal::Matrix< double> data_a( GenerateDescription( SEQUENCE_A));
      const linal::Matrix< double> data_b( GenerateDescription( SEQUENCE_B));

      // temporary vector for storing 5 predictions between 2 AA
      linal::Vector< double> predictions( Types::s_NumberValidTypes);

      // iterate over every residue h
      for
      (
        size_t aa_ctr_a( s_NumberSkippedResiduesInEnds);
        aa_ctr_a < SEQUENCE_A.GetSize() - s_NumberSkippedResiduesInEnds;
        ++aa_ctr_a
      )
      {
        // and every other following residue i
        for
        (
          size_t aa_ctr_b( s_NumberSkippedResiduesInEnds);
          aa_ctr_b < SEQUENCE_B.GetSize() - s_NumberSkippedResiduesInEnds;
          ++aa_ctr_b
        )
        {
          // accumulate the predictions for each contacttype
          for( size_t type = GetTypes().HELIX_HELIX; type < Types::s_NumberValidTypes; ++type)
          {
            Type contact_type( type);

            predictions( contact_type) =
              PredictSingleContact
              (
                contact_type,
                linal::Vector< double>
                (
                  s_NumberOfValuesPerResidue * ( contact_type->GetWindowLengths().First()),
                  data_a[ aa_ctr_a - contact_type->GetWindowRadii().First()]
                ),
                linal::Vector< double>
                (
                  s_NumberOfValuesPerResidue * ( contact_type->GetWindowLengths().Second()),
                  data_b[ aa_ctr_b - contact_type->GetWindowRadii().Second()]
                ),
                linal::Vector3D( 40, 80, 40)
              );
          }

          // insert the predictions for the selected residue pair
          InsertPredictions
          (
            *SEQUENCE_A.GetAA( aa_ctr_a), *SEQUENCE_B.GetAA( aa_ctr_b), predictions
          );

          // insert the predictions for the reversed residue pair
          InsertPredictions
          (
            *SEQUENCE_B.GetAA( aa_ctr_b), *SEQUENCE_A.GetAA( aa_ctr_a), ReversePredictions( predictions)
          );

        } // loop over i
      } // loop over h
    }

    //! @brief this functions checks that jufo and blast profiles are stored for every sequence provided
    void PredictionMap::CheckBlastAndJufoForSequences() const
    {
      // iterate over every sequence stored
      for
      (
        util::ShPtrVector< biol::AASequence>::const_iterator seq_itr( m_Sequences.Begin()),
          seq_itr_end( m_Sequences.End());
        seq_itr != seq_itr_end;
        ++seq_itr
      )
      {
         //check that blast profile is stored
         BCL_Assert
         (
           ( *seq_itr)->GetFirstAA()->GetBlastProfilePtr().IsDefined(),
           "Blast Profile is not available for chain " + util::Format()( ( *seq_itr)->GetChainID())
         );

         //check that jufo is stored
         BCL_Assert
         (
           ( *seq_itr)->GetFirstAA()->GetSSPrediction( sspred::GetMethods().e_JUFO).IsDefined(),
           "Jufo Profile is not available for chain " + util::Format()( ( *seq_itr)->GetChainID())
         );
      }
    }

  } // namespace contact
} // namespace bcl
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
#include "contact/bcl_contact_recovery.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_list_container_generator_protein_model.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Recovery::s_Instance
    (
      GetObjectInstances().AddInstance( new Recovery())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a minimum sequence separation and cache boolean
    //! @param MIN_SEQUENCE_SEPARATION Minimum sequence separation
    //! @param CACHE boolean to decide whether a neighbor generator with cache should be used
    Recovery::Recovery
    (
      const size_t MIN_SEQUENCE_SEPARATION,
      const math::ContingencyMatrixMeasures::MeasureEnum &MEASURE
    ) :
      m_MinSequenceSeparation( MIN_SEQUENCE_SEPARATION),
      m_NeighborGenerator
      (
        assemble::AANeighborListContainerGeneratorProteinModel::AANeighborListGenerator
        (
          g_ContactCbDistanceCutoff,
          g_ContactMinSequenceSeparation,
          true,
          true
        )
      ),
      m_Measure( MEASURE),
      m_Normalization
      (
        util::EndsWith( math::ContingencyMatrixMeasures::GetMeasureName( m_Measure.GetMeasure()), "R")
        ? 100.0
        : 1.0
      )
      {
        BCL_Assert
        (
          MIN_SEQUENCE_SEPARATION >= g_ContactMinSequenceSeparation,
          "Given sequence separation " + util::Format()( m_MinSequenceSeparation) + " can't be smaller than " +
          util::Format()( g_ContactMinSequenceSeparation)
        );
      }

      //! @brief Clone function
      //! @return pointer to new Recovery
      Recovery *Recovery::Clone() const
      {
        return new Recovery( *this);
      }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Recovery::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate and return the percentage of recovered contacts
    //! @param TEMPLATE_MODEL Template model of interest
    //! @param PROTEIN_MODEL ProteinModel to be evaluated
    //! @return the percentage of recovered contacts
    double Recovery::operator()
    (
      const assemble::ProteinModel &TEMPLATE_MODEL,
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // generate the neighbor list containers for both models
      assemble::AANeighborListContainer neighbors_template( m_NeighborGenerator->operator()( TEMPLATE_MODEL));
      assemble::AANeighborListContainer neighbors_model( m_NeighborGenerator->operator()( PROTEIN_MODEL));

      // prune both by the min sequence separation
      if( m_MinSequenceSeparation > g_ContactMinSequenceSeparation)
      {
        neighbors_template.Prune( g_ContactCbDistanceCutoff, m_MinSequenceSeparation, true);
        neighbors_model.Prune( g_ContactCbDistanceCutoff, m_MinSequenceSeparation, true);
      }

      // calculate the contingency matrix
      const math::ContingencyMatrix contingency_matrix( CalculateContingencyMatrix( neighbors_template, neighbors_model));

      // if there were no positives found
      if( contingency_matrix.GetNumberActualPositives() == 0)
      {
        BCL_MessageStd( "No contacts found in the template model, returning 0");
        return double( 0.0);
      }

      return m_Normalization * m_Measure( contingency_matrix);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Recovery::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_MinSequenceSeparation, ISTREAM);
      io::Serialize::Read( m_Measure, ISTREAM);
      *this = Recovery( m_MinSequenceSeparation, m_Measure.GetMeasure());

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Recovery::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_MinSequenceSeparation, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Measure, OSTREAM, INDENT) << '\n';

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculate and return the contingency matrix for two neighbor list containers
    //! @param TEMPLATE_NEIGHBOR_CONTAINER AANeighborListContainer for a template model
    //! @param MODEL_NEIGHBOR_CONTAINER AANeighborListContainer for the model to be evaluated
    //! @return calculate and return the contingency matrix for two neighbor list containers
    math::ContingencyMatrix Recovery::CalculateContingencyMatrix
    (
      const assemble::AANeighborListContainer &TEMPLATE_NEIGHBOR_CONTAINER,
      const assemble::AANeighborListContainer &MODEL_NEIGHBOR_CONTAINER
    )
    {
      // store the number of contacts for template and model
      const size_t nr_contacts_template( TEMPLATE_NEIGHBOR_CONTAINER.GetNumberNeighbors());
      const size_t nr_contacts_model( MODEL_NEIGHBOR_CONTAINER.GetNumberNeighbors());
      const size_t nr_common_contacts( TEMPLATE_NEIGHBOR_CONTAINER.IntersectionSize( MODEL_NEIGHBOR_CONTAINER));
      const size_t min_seq_sep( TEMPLATE_NEIGHBOR_CONTAINER.GetMinimalSequenceSeparation());
      const size_t chain_size
      (
        std::max( std::max( MODEL_NEIGHBOR_CONTAINER.GetSize(), size_t( min_seq_sep + 1)), nr_contacts_model)
      );

      // initialize a contingency matrix from the intersect values
      const math::ContingencyMatrix matrix
      (
        nr_common_contacts,                        // true positives
        nr_contacts_model - nr_common_contacts,    // false positives
        nr_contacts_template - nr_common_contacts, // false negatives
        chain_size * ( chain_size - min_seq_sep - 1) - nr_contacts_model // true negatives
      );

      // end
      return matrix;
    }

  } // namespace contact
} // namespace bcl
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
#include "contact/bcl_contact_sse_prediction_map.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SSEPredictionMap::SSEPredictionMap() :
      storage::ObjectNDHashMap< 4, biol::AAData, double>(),
      m_PredictionMap( new PredictionMap())
    {
    }

    //! @brief constructor from a prediction map and ShPtrVector of SSEs
    //! @param PREDICTION_MAP PredictionMap that contains contact predictions for amino acids
    //! @param SSE_POOL ShPtrVector of SSEs of interest
    SSEPredictionMap::SSEPredictionMap
    (
      const util::ShPtr< PredictionMap> &PREDICTION_MAP,
      const util::ShPtrVector< assemble::SSE> &SSE_POOL
    ) :
      storage::ObjectNDHashMap< 4, biol::AAData, double>(),
      m_PredictionMap( PREDICTION_MAP)
    {
      FillMap( SSE_POOL);
    }

    //! @brief virtual copy constructor
    SSEPredictionMap *SSEPredictionMap::Clone() const
    {
      return new SSEPredictionMap( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &SSEPredictionMap::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns the prediction stored for given SSEs
    //! @param SSE_A first SSE of interest
    //! @param SSE_B second SSE of interest
    //! @return prediction stored for SSE_A and SSE_B
    double SSEPredictionMap::GetPrediction( const assemble::SSE &SSE_A, const assemble::SSE &SSE_B) const
    {
      storage::VectorND< 4, util::SiPtr< const biol::AAData> > keys
      (
        util::SiPtr< const biol::AAData>( SSE_A.GetFirstAA()->GetData()),
        util::SiPtr< const biol::AAData>( SSE_A.GetLastAA()->GetData()),
        util::SiPtr< const biol::AAData>( SSE_B.GetFirstAA()->GetData()),
        util::SiPtr< const biol::AAData>( SSE_B.GetLastAA()->GetData())
      );

      storage::HashMap< size_t, double>::const_iterator itr = Find( keys);
      if( itr == End())
      {
        return 0;
      }
      return itr->second;
    }

    //! @brief iterates over every SSE pair for a given ShPtrVector of SSEs and calculates the probabilities
    //! @param SSE_POOL ShPtrVector of SSEs of interest
    void SSEPredictionMap::FillMap( const util::ShPtrVector< assemble::SSE> &SSE_POOL)
    {
      // iterate over every sse
      for
      (
        util::ShPtrVector< assemble::SSE>::const_iterator sse_itr_a( SSE_POOL.Begin()), sse_itr_end( SSE_POOL.End());
        sse_itr_a != sse_itr_end;
        ++sse_itr_a
      )
      {
        // versus every other sse
        for
        (
          util::ShPtrVector< assemble::SSE>::const_iterator sse_itr_b( SSE_POOL.Begin());
          sse_itr_b != sse_itr_end;
          ++sse_itr_b
        )
        {
          // if they are not the same SSE or do not overlap
          if( sse_itr_a != sse_itr_b && !biol::DoOverlap( **sse_itr_a, **sse_itr_b))
          {
            // value to store the sum of all predictions of every possible residue pairing between t two sses
            double sum( CalculateContactProbability( **sse_itr_a, **sse_itr_b));

            storage::VectorND< 4, util::SiPtr< const biol::AAData> > keys
            (
              util::SiPtr< const biol::AAData>( ( *sse_itr_a)->GetFirstAA()->GetData()),
              util::SiPtr< const biol::AAData>( ( *sse_itr_a)->GetLastAA()->GetData()),
              util::SiPtr< const biol::AAData>( ( *sse_itr_b)->GetFirstAA()->GetData()),
              util::SiPtr< const biol::AAData>( ( *sse_itr_b)->GetLastAA()->GetData())
            );

            // insert this value with corresponding keys into the ObjectNDHashMap
            Insert( keys, sum);
          }
        }
      }
    }

    //! @brief calculates the contact probability given two SSEs
    //! @param SSE_A first SSE of interest
    //! @param SSE_B second SSE of interest
    //! @return contact probability for given pair of SSEs
    double SSEPredictionMap::CalculateContactProbability( const assemble::SSE &SSE_A, const assemble::SSE &SSE_B) const
    {
      // value to store the sum of all predictions of every possible residue pairing between t two sses
      double sum( 0.0);

      // value to store number of residue pairs for which no predictions existed
      size_t no_predictions( 0);

      // contact type of these two sses
      Type contact_type( GetTypes().TypeFromSSTypes( SSE_A, SSE_B));

      // iterate over every residue in SSE_A
      for
      (
        biol::AASequence::const_iterator aa_itr_a( SSE_A.GetData().Begin()),
          aa_itr_a_end( SSE_A.GetData().End());
        aa_itr_a != aa_itr_a_end;
        ++aa_itr_a
      )
      {
        // iterate over every residue in SSE_B
        for
        (
          biol::AASequence::const_iterator aa_itr_b( SSE_B.GetData().Begin()),
            aa_itr_b_end( SSE_B.GetData().End());
          aa_itr_b != aa_itr_b_end;
          ++aa_itr_b
        )
        {
          storage::VectorND< 2, util::SiPtr< const biol::AAData> > key
          (
            util::SiPtr< const biol::AAData>( ( *aa_itr_a)->GetData()),
            util::SiPtr< const biol::AAData>( ( *aa_itr_b)->GetData())
          );
          double prediction( ( m_PredictionMap->GetPredictions( key)).First()( contact_type));

          if( !util::IsDefined( prediction))
          {
            no_predictions++;
          }
          else
          {
            sum += prediction;
          }
        }
      }

      // return summed probability normalized with length of SSEs
      return sum / ( ( SSE_A.GetData().GetSize() * SSE_B.GetData().GetSize()) - no_predictions);
    }

  } // namespace contact
} // namespace bcl
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
#include "contact/bcl_contact_statistics.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container_generator_protein_model.h"
#include "biol/bcl_biol_aa_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {

    //! @brief StatisticType as string
    //! @param STATISTIC_TYPE the StatisticType
    //! @return the string for the StatisticType
    const std::string &Statistics::GetStatisticTypeDescriptor( const StatisticType &STATISTIC_TYPE)
    {
      static const std::string s_descriptors[] =
      {
        "NumberContacts",
        "NumberContactsShort",
        "NumberContactsMid",
        "NumberContactsLong",
        "RatioContactsShort",
        "RatioContactsMid",
        "RatioContactsLong",
        GetStaticClassName< StatisticType>()
      };

      return s_descriptors[ STATISTIC_TYPE];
    }

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> Statistics::s_Instance
    (
      GetObjectInstances().AddInstance( new Statistics())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Statistics::Statistics() :
      m_SequenceSeparationRange(),
      m_UseRatio(),
      m_NeighborGenerator()
    {
    }

    //! @brief constructor from a statistics type to be calculated and cache boolean
    //! @param STATISTIC_TYPE StatisticType to be calculated
    //! @param CACHE whether a neighbor generator with cache should be used
    Statistics::Statistics
    (
      const StatisticType &STATISTIC_TYPE,
      const bool CACHE
    ) :
      m_SequenceSeparationRange(),
      m_UseRatio(),
      m_NeighborGenerator
      (
        assemble::AANeighborListContainerGeneratorProteinModel::AANeighborListGenerator
        (
          g_ContactCbDistanceCutoff, g_ContactMinSequenceSeparation, true, CACHE
        )
      )
    {
      // switch over the statistic type given
      switch( STATISTIC_TYPE)
      {
        case e_NumberContacts:
        {
          m_SequenceSeparationRange = GetDefaultSequenceSeparationRange();
          m_UseRatio = false;
          break;
        }
        case e_NumberContactsShort:
        {
          m_SequenceSeparationRange = GetDefaultSequenceSeparationShortRange();
          m_UseRatio = false;
          break;
        }
        case e_NumberContactsMid:
        {
          m_SequenceSeparationRange = GetDefaultSequenceSeparationMidRange();
          m_UseRatio = false;
          break;
        }
        case e_NumberContactsLong:
        {
          m_SequenceSeparationRange = GetDefaultSequenceSeparationLongRange();
          m_UseRatio = false;
          break;
        }
        case e_RatioContactsShort:
        {
          m_SequenceSeparationRange = GetDefaultSequenceSeparationShortRange();
          m_UseRatio = true;
          break;
        }
        case e_RatioContactsMid:
        {
          m_SequenceSeparationRange = GetDefaultSequenceSeparationMidRange();
          m_UseRatio = true;
          break;
        }
        case e_RatioContactsLong:
        {
          m_SequenceSeparationRange = GetDefaultSequenceSeparationLongRange();
          m_UseRatio = true;
          break;
        }
        case s_NumberStatisticType:
        {
          break;
        }
      }
    }

    //! @brief constructor from a sequence separation range, whether to calculate ratio and whether to use cache
    //! @param SEQUENCE_SEPARATION_RANGE Sequence separation range
    //! @param CALCULATE_RATIO boolean indicating whether just the counts or the ratio wrt to total number of contacts should be calculated
    //! @param CACHE whether a neighbor generator with cache should be used
    Statistics::Statistics
    (
      const math::Range< size_t> &SEQUENCE_SEPARATION_RANGE,
      const bool CALCULATE_RATIO,
      const bool CACHE
    ) :
      m_SequenceSeparationRange( SEQUENCE_SEPARATION_RANGE),
      m_UseRatio( CALCULATE_RATIO),
      m_NeighborGenerator
      (
        assemble::AANeighborListContainerGeneratorProteinModel::AANeighborListGenerator
        (
          g_ContactCbDistanceCutoff, g_ContactMinSequenceSeparation, true, CACHE
        )
      )
      {
      }

    //! @brief Clone function
    //! @return pointer to new Statistics
    Statistics *Statistics::Clone() const
    {
      return new Statistics( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Statistics::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  /////////////////
  // operations  //
  /////////////////

    //! @brief calculate the contact statistics value for the given AANeighborListContainer
    //! @param NEIGHBOR_CONTAINER AANeighborListContainer of interest
    //! @return calculated statistics value for the given AANeighborListContainer
    double Statistics::Calculate( const assemble::AANeighborListContainer &NEIGHBOR_CONTAINER) const
    {
      // make a copy of the container
      assemble::AANeighborListContainer container
      (
        NEIGHBOR_CONTAINER,
        g_ContactCbDistanceCutoff,
        g_ContactMinSequenceSeparation,
        true
      );

      // if ratio is asked
      if( m_UseRatio)
      {
        return GetContactsRatio( m_SequenceSeparationRange, container);
      }
      // otherwise just return counts
      else
      {
        return double( GetNumberContacts( m_SequenceSeparationRange, container) / 2);
      }
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate the contact statistics value for the given ProteinModel
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @return calculuated statistics value for the given model
    double Statistics::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // calculate the neighbor list for the given ProteinModel and calculate and return the statistic value
      return Calculate( m_NeighborGenerator->operator()( PROTEIN_MODEL));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Statistics::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SequenceSeparationRange, ISTREAM);
      io::Serialize::Read( m_UseRatio, ISTREAM);
      io::Serialize::Read( m_NeighborGenerator, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Statistics::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SequenceSeparationRange, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_UseRatio, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NeighborGenerator, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculate the number of contacts within the given sequence separation range for the given NeighborListContainer
    //! @param SEQUENCE_SEPARATION_RANGE Sequence separation range to be used
    //! @param NEIGHBOR_CONTAINER AANeighborListContainer of interest
    //! @return the number of contacts in the given NeighborContainer
    size_t Statistics::GetNumberContacts
    (
      const math::Range< size_t> &SEQUENCE_SEPARATION_RANGE,
      const assemble::AANeighborListContainer &NEIGHBOR_CONTAINER
    )
    {
      // initialize count
      size_t count( 0);

      // iterate over the neighbor list container
      for
      (
        assemble::AANeighborListContainer::const_iterator
          itr( NEIGHBOR_CONTAINER.Begin()), itr_end( NEIGHBOR_CONTAINER.End());
        itr != itr_end; ++itr
      )
      {
        // create a reference on the center amino acid
        const biol::AABase &center_aa( *itr->second.GetCenterAminoAcid());

        // iterate over the neighbors listed
        for
        (
          assemble::AANeighborList::const_iterator neigh_itr( itr->second.Begin()), neigh_itr_end( itr->second.End());
          neigh_itr != neigh_itr_end; ++neigh_itr
        )
        {
          // if it's a different chain or the sequence separation is within the range then skip over
          if
          (
            center_aa.GetChainID() != neigh_itr->First()->GetChainID() ||
            SEQUENCE_SEPARATION_RANGE.IsWithin( biol::SequenceSeparation( center_aa, *neigh_itr->First()))
          )
          {
            // then increment count
            ++count;
          }
        }
      }

      // end
      return count;
    }

    //! @brief calculate and return the ratio of all contacts that are within given sequence separation range
    //! @param SEQUENCE_SEPARATION_RANGE sequence separation range to be used for calculating ratio
    //! @param NEIGHBOR_CONTAINER AANeighborListContainer of interest
    //! @return the ratio of all contacts that are within given sequence separation range
    double Statistics::GetContactsRatio
    (
      const math::Range< size_t> &SEQUENCE_SEPARATION_RANGE,
      const assemble::AANeighborListContainer &NEIGHBOR_CONTAINER
    )
    {
      // calculate the counts with the sequence separation
      const double counts_range( GetNumberContacts( SEQUENCE_SEPARATION_RANGE, NEIGHBOR_CONTAINER));

      // divide by total number of counts and return
      return 100.0 * counts_range / NEIGHBOR_CONTAINER.GetNumberNeighbors();
    }

  } // namespace contact
} // namespace bcl
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
#include "contact/bcl_contact_type_data.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> TypeData::s_Instance( GetObjectInstances().AddInstance( new TypeData()));

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct undefined contact type
    TypeData::TypeData() :
      m_WindowRadii( util::GetUndefined< size_t>(), util::GetUndefined< size_t>()),
      m_WindowLengths( util::GetUndefined< size_t>(), util::GetUndefined< size_t>()),
      m_IsValid( false),
      m_ResidueDistanceCutoff( util::GetUndefinedDouble()),
      m_MinimalSSEDistance( util::GetUndefinedDouble()),
      m_DistanceRange(),
      m_PreferredDistanceRange(),
      m_TiltAngleRange(),
      m_MinimalFragmentInterfaceLength( util::GetUndefinedDouble()),
      m_SSTypes( biol::GetSSTypes().e_Undefined)
    {
    }

    //! @brief construct contact type from provided data
    //! @param WINDOW_RADII radius of window used for representation
    //! @param WINDOW_LENGTH_PAIR window length used for representation
    //! @param IS_VALID whether a valid contact type
    //! @param RESIDUE_DISTANCE_CUTOFF cut off distance
    //! @param MINIMAL_SSE_DISTANCE minimal distance of two sses to not clash
    //! @param MINIMAL_FRAGMENT_INTERFACE_LENGTH minimal fragment interface length
    //! @param DISTANCE_RANGE distance range
    //! @param PREFERRED_DISTANCE_RANGE preferred distance range
    //! @param TILT_ANGLE_RANGE tilt angle range
    //! @param SS_TYPES SSTypes associated with this contact type
    TypeData::TypeData
    (
      const storage::Pair< size_t, size_t> &WINDOW_RADII,
      const storage::Pair< size_t, size_t> &WINDOW_LENGTH_PAIR,
      const bool IS_VALID,
      const double RESIDUE_DISTANCE_CUTOFF,
      const double MINIMAL_SSE_DISTANCE,
      const math::Range< double> &DISTANCE_RANGE,
      const math::Range< double> &PREFERRED_DISTANCE_RANGE,
      const math::Range< double> &TILT_ANGLE_RANGE,
      const double MINIMAL_FRAGMENT_INTERFACE_LENGTH,
      const storage::Set< biol::SSType> &SS_TYPES
    ) :
      m_WindowRadii( WINDOW_RADII),
      m_WindowLengths( WINDOW_LENGTH_PAIR),
      m_IsValid( IS_VALID),
      m_ResidueDistanceCutoff( RESIDUE_DISTANCE_CUTOFF),
      m_MinimalSSEDistance( MINIMAL_SSE_DISTANCE),
      m_DistanceRange( DISTANCE_RANGE),
      m_PreferredDistanceRange( PREFERRED_DISTANCE_RANGE),
      m_TiltAngleRange( TILT_ANGLE_RANGE),
      m_MinimalFragmentInterfaceLength( MINIMAL_FRAGMENT_INTERFACE_LENGTH),
      m_SSTypes( SS_TYPES)
    {
    }

    //! @brief virtual copy constructor
    TypeData *TypeData::Clone() const
    {
      return new TypeData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &TypeData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &TypeData::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_WindowRadii, ISTREAM);
      io::Serialize::Read( m_WindowLengths, ISTREAM);
      io::Serialize::Read( m_IsValid, ISTREAM);
      io::Serialize::Read( m_ResidueDistanceCutoff, ISTREAM);
      io::Serialize::Read( m_DistanceRange, ISTREAM);
      io::Serialize::Read( m_PreferredDistanceRange, ISTREAM);
      io::Serialize::Read( m_TiltAngleRange, ISTREAM);
      io::Serialize::Read( m_MinimalFragmentInterfaceLength, ISTREAM);
      io::Serialize::Read( m_SSTypes, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @return ostream which was written to
    std::ostream &TypeData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_WindowRadii, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_WindowLengths, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_IsValid, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ResidueDistanceCutoff, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DistanceRange, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PreferredDistanceRange, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_TiltAngleRange, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MinimalFragmentInterfaceLength, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SSTypes, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace contact
} // namespace bcl
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
#include "contact/bcl_contact_types.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {

    //! @brief default constructor that constructs all Types
                     //                                                                                                                                                                                                                                                                                                                                          is_valid  dist   min  distance                          preferred distance                   tilt angle                                                                      min frag.                  sstypes
    Types::Types() : //                                                      WINDOW_RADII                                                                                                                               WINDOW_LENGTH_PAIR                                                                                                                       cutoff ssedist        range                             range                                range                                                                           intlength
      HELIX_HELIX(             AddEnum( "HELIX_HELIX",             TypeData( storage::Pair< size_t, size_t>( biol::GetSSTypes().HELIX->GetContactWindowRadius(),  biol::GetSSTypes().HELIX->GetContactWindowRadius()),  storage::Pair< size_t, size_t>( biol::GetSSTypes().HELIX->GetContactWindowLength(),  biol::GetSSTypes().HELIX->GetContactWindowLength()),  true,   8.0,   4.0, math::Range< double>( 5.0, 16.0), math::Range< double>( 8.00,  12.00), math::Range< double>( math::Angle::Radian( -60.0), math::Angle::Radian( 45.0)), 4.0, storage::Set< biol::SSType>( biol::GetSSTypes().HELIX)                           ))),
      HELIX_SHEET(             AddEnum( "HELIX_SHEET",             TypeData( storage::Pair< size_t, size_t>( biol::GetSSTypes().HELIX->GetContactWindowRadius(),  biol::GetSSTypes().STRAND->GetContactWindowRadius()), storage::Pair< size_t, size_t>( biol::GetSSTypes().HELIX->GetContactWindowLength(),  biol::GetSSTypes().STRAND->GetContactWindowLength()), true,   8.0,   4.0, math::Range< double>( 5.0, 16.0), math::Range< double>( 8.00,  13.00), math::Range< double>( math::Angle::Radian( -60.0), math::Angle::Radian(  0.0)), 4.0, storage::Set< biol::SSType>( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND)))),
      SHEET_HELIX(             AddEnum( "SHEET_HELIX",             TypeData( storage::Pair< size_t, size_t>( biol::GetSSTypes().STRAND->GetContactWindowRadius(), biol::GetSSTypes().HELIX->GetContactWindowRadius()),  storage::Pair< size_t, size_t>( biol::GetSSTypes().STRAND->GetContactWindowLength(), biol::GetSSTypes().HELIX->GetContactWindowLength()),  true,   8.0,   4.0, math::Range< double>( 5.0, 16.0), math::Range< double>( 8.00,  13.00), math::Range< double>( math::Angle::Radian( -60.0), math::Angle::Radian(  0.0)), 4.0, storage::Set< biol::SSType>( biol::GetSSTypes().STRAND, biol::GetSSTypes().HELIX)))),
      HELIX_STRAND(            AddEnum( "HELIX_STRAND",            TypeData( storage::Pair< size_t, size_t>( biol::GetSSTypes().HELIX->GetContactWindowRadius(),  biol::GetSSTypes().STRAND->GetContactWindowRadius()), storage::Pair< size_t, size_t>( biol::GetSSTypes().HELIX->GetContactWindowLength(),  biol::GetSSTypes().STRAND->GetContactWindowLength()), true,   8.0,   4.0, math::Range< double>( 7.0, 18.0), math::Range< double>( 9.00,  16.00), math::Range< double>( math::Angle::Radian( -60.0), math::Angle::Radian(  0.0)), 4.0, storage::Set< biol::SSType>( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND)))),
      STRAND_HELIX(            AddEnum( "STRAND_HELIX",            TypeData( storage::Pair< size_t, size_t>( biol::GetSSTypes().STRAND->GetContactWindowRadius(), biol::GetSSTypes().HELIX->GetContactWindowRadius()),  storage::Pair< size_t, size_t>( biol::GetSSTypes().STRAND->GetContactWindowLength(), biol::GetSSTypes().HELIX->GetContactWindowLength()),  true,   8.0,   4.0, math::Range< double>( 7.0, 18.0), math::Range< double>( 9.00,  16.00), math::Range< double>( math::Angle::Radian( -60.0), math::Angle::Radian(  0.0)), 4.0, storage::Set< biol::SSType>( biol::GetSSTypes().STRAND, biol::GetSSTypes().HELIX)))),
      STRAND_STRAND(           AddEnum( "STRAND_STRAND",           TypeData( storage::Pair< size_t, size_t>( biol::GetSSTypes().STRAND->GetContactWindowRadius(), biol::GetSSTypes().STRAND->GetContactWindowRadius()), storage::Pair< size_t, size_t>( biol::GetSSTypes().STRAND->GetContactWindowLength(), biol::GetSSTypes().STRAND->GetContactWindowLength()), true,   8.0,   3.0, math::Range< double>( 3.5,  5.5), math::Range< double>( 4.00,   5.25), math::Range< double>( math::Angle::Radian( -30.0), math::Angle::Radian(  0.0)), 4.0, storage::Set< biol::SSType>( biol::GetSSTypes().STRAND)                          ))),
      SHEET_SHEET(             AddEnum( "SHEET_SHEET",             TypeData( storage::Pair< size_t, size_t>( biol::GetSSTypes().STRAND->GetContactWindowRadius(), biol::GetSSTypes().STRAND->GetContactWindowRadius()), storage::Pair< size_t, size_t>( biol::GetSSTypes().STRAND->GetContactWindowLength(), biol::GetSSTypes().STRAND->GetContactWindowLength()), true,   8.0,   3.0, math::Range< double>( 6.0, 14.0), math::Range< double>( 8.00,  12.00), math::Range< double>( math::Angle::Radian( -60.0), math::Angle::Radian(  0.0)), 4.0, storage::Set< biol::SSType>( biol::GetSSTypes().STRAND)                          ))),
      UNDEFINED_HELIX_STRAND(  AddEnum( "UNDEFINED_HELIX_STRAND",  TypeData( GetUndefinedLengthPair(),                                                                                                                  GetUndefinedLengthPair(),                                                                                       false, GetUnknownResidueDistanceCutoff(), 4.0, GetUnknownDistanceRange()       , GetUnknownPreferredDistanceRange() , GetUnknownTiltAngleRange()                                                    , 4.0, storage::Set< biol::SSType>( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND)))),
      UNDEFINED_STRAND_HELIX(  AddEnum( "UNDEFINED_STRAND_HELIX",  TypeData( GetUndefinedLengthPair(),                                                                                                                  GetUndefinedLengthPair(),                                                                                       false, GetUnknownResidueDistanceCutoff(), 4.0, GetUnknownDistanceRange()       , GetUnknownPreferredDistanceRange() , GetUnknownTiltAngleRange()                                                    , 4.0, storage::Set< biol::SSType>( biol::GetSSTypes().STRAND, biol::GetSSTypes().HELIX)))),
      UNDEFINED_STRAND_STRAND( AddEnum( "UNDEFINED_STRAND_STRAND", TypeData( GetUndefinedLengthPair(),                                                                                                                  GetUndefinedLengthPair(),                                                                                       false, GetUnknownResidueDistanceCutoff(), 3.0, GetUnknownDistanceRange()       , GetUnknownPreferredDistanceRange() , GetUnknownTiltAngleRange()                                                    , 4.0, storage::Set< biol::SSType>( biol::GetSSTypes().STRAND)                          )))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Types::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief cut off distance for unknown types;
    //! @return cut off distance for unknown types;
    const double &Types::GetUnknownResidueDistanceCutoff()
    {
      static const double s_distance_cutoff( 8.0);

      // end
      return s_distance_cutoff;
    }

    //! @brief return distance range for unknown types
    //! @return distance range for unknown types
    const math::Range< double> &Types::GetUnknownDistanceRange()
    {
      static const math::Range< double> s_distance_range( 5.0, 14.0);

      // end
      return s_distance_range;
    }

    //! @brief return preferred distance range for unknown types
    //! @return preferred distance range for unknown types
    const math::Range< double> &Types::GetUnknownPreferredDistanceRange()
    {
      static const math::Range< double> s_distance_range( 5.0, 14.0);

      // end
      return s_distance_range;
    }

    //! @brief return tilt angle range for unknown types
    //! @return tilt angle range for unknown types
    const math::Range< double> &Types::GetUnknownTiltAngleRange()
    {
      static const math::Range< double> s_tilt_angle_range( math::Angle::Radian( -30.0), math::Angle::Radian( 30.0));

      // end
      return s_tilt_angle_range;
    }

    //! @brief return undefined pair length
    //! @return undefined pair length
    const storage::Pair< size_t, size_t> &Types::GetUndefinedLengthPair()
    {
      static const storage::Pair< size_t, size_t> s_undefined_pair
      (
        util::GetUndefined< size_t>(), util::GetUndefined< size_t>()
      );

      // end
      return s_undefined_pair;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief returns the reverse pair for the contact type
    //! @param TYPE Type of interest
    //! @return the reverse pair for the contact type
    const Type &Types::Reverse( const Type &TYPE) const
    {
      if( TYPE == HELIX_SHEET)
      {
        return SHEET_HELIX;
      }
      else if( TYPE == SHEET_HELIX)
      {
        return HELIX_SHEET;
      }
      if( TYPE == HELIX_STRAND)
      {
        return STRAND_HELIX;
      }
      else if( TYPE == STRAND_HELIX)
      {
        return HELIX_STRAND;
      }
      else if( TYPE == UNDEFINED_HELIX_STRAND)
      {
        return UNDEFINED_STRAND_HELIX;
      }
      else if( TYPE == UNDEFINED_STRAND_HELIX)
      {
        return UNDEFINED_HELIX_STRAND;
      }

      return TYPE;
    }

    //! @brief Given two SSEGeometryInterfaces, merges their SSTypes to form a ContactType
    //! @param SSE_A first SSEGeometryInterface
    //! @param SSE_B second SSEGeometryInterface
    //! @return Contact type formed by SSE_A and SSE_B
    const Type &Types::TypeFromSSTypes
    (
      const assemble::SSEGeometryInterface &SSE_A,
      const assemble::SSEGeometryInterface &SSE_B
    ) const
    {
      // HELIX_HELIX
      if( SSE_A.GetType() == biol::GetSSTypes().HELIX && SSE_B.GetType() == biol::GetSSTypes().HELIX)
      {
        return HELIX_HELIX;
      }
      // HELIX_SHEET
      if( SSE_A.GetType() == biol::GetSSTypes().HELIX && SSE_B.GetType() == biol::GetSSTypes().STRAND)
      {
        return HELIX_SHEET;
      }
      // HELIX_SHEET
      if( SSE_A.GetType() == biol::GetSSTypes().STRAND && SSE_B.GetType() == biol::GetSSTypes().HELIX)
      {
        return SHEET_HELIX;
      }
      // STRAND_STRAND ( or SHEET_SHEET)
      if( SSE_A.GetType() == biol::GetSSTypes().STRAND && SSE_B.GetType() == biol::GetSSTypes().STRAND)
      {
        return STRAND_STRAND;
      }

      return e_Undefined;
    }

    //! @brief returns a map with distance ranges for valid contact types
    //! @return a map with distance ranges for valid contact types
    const storage::Map< Type, math::Range< double> > &Types::GetValidDistanceRanges() const
    {
      // initialize static storage for the distance cutoffs
      static const storage::Map< Type, math::Range< double> > s_distance_ranges( CollectValidDistanceRanges());

      // end
      return s_distance_ranges;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief generates a map with distance ranges for valid contact types
    //! @return a map with distance ranges for valid contact types
    storage::Map< Type, math::Range< double> > Types::CollectValidDistanceRanges() const
    {
      // initialize storage
      storage::Map< Type, math::Range< double> > distance_range_map;

      // iterator over standard amino acid types
      for
      (
        const_iterator type_itr( Begin()), type_itr_end( GetEnumIteratorFromIndex( s_NumberValidTypes));
        type_itr != type_itr_end;
        ++type_itr
      )
      {
        // store the cutoff distance thresholds
        distance_range_map[ *type_itr] = ( *type_itr)->GetDistanceRange();
      }
      // end
      return distance_range_map;
    }

    //! @brief retrieves Types enumerator
    //! @return Types enumerator
    const Types &GetTypes()
    {
      return Types::GetEnums();
    }

  } // namespace contact

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< contact::TypeData, contact::Types>;

  } // namespace util
} // namespace bcl
