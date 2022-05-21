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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "contact/bcl_contact_correlation_matrix.h"

// includes from bcl - sorted alphabetically
#include "align/bcl_align_aligner_progressive.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "contact/bcl_contact_calculate_correlations_sm.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_contact_correlation_matrix.cpp
  //! @brief TODO: add an detailed description for this example
  //!
  //! @author teixeipl
  //! @date Aug 1, 2011
  //! @remarks status incomplete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleContactCorrelationMatrix :
    public ExampleInterface
  {
  public:

    ExampleContactCorrelationMatrix *Clone() const
    {
      return new ExampleContactCorrelationMatrix( *this);
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

    int Run() const
    {
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // Create correlation matrix using alignments
      io::IFStream readstream;
      // read seq1 from fasta
      BCL_ExampleMustOpenInputFile( readstream, AddExampleInputPathToFilename( e_Biology, "1EPWtest1.fasta"));
      util::ShPtr< biol::AASequence> seq1
      (
        new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( readstream))
      );
      io::File::CloseClearFStream( readstream);
      // read seq2 from fasta
      BCL_ExampleMustOpenInputFile( readstream, AddExampleInputPathToFilename( e_Biology, "1EPWtest2.fasta"));
      util::ShPtr< biol::AASequence> seq2
      (
        new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( readstream))
      );
      io::File::CloseClearFStream( readstream);
      // read seq3 from fasta
      BCL_ExampleMustOpenInputFile( readstream, AddExampleInputPathToFilename( e_Biology, "1EPWtest3.fasta"));
      util::ShPtr< biol::AASequence> seq3
      (
        new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( readstream))
      );
      io::File::CloseClearFStream( readstream);

      // create alignment, alignment_list and engine
      util::ShPtr< align::AlignmentInterface< biol::AABase> > alignment1( new align::AlignmentLeaf< biol::AABase>( seq1));
      util::ShPtr< align::AlignmentInterface< biol::AABase> > alignment2( new align::AlignmentLeaf< biol::AABase>( seq2));
      util::ShPtr< align::AlignmentInterface< biol::AABase> > alignment3( new align::AlignmentLeaf< biol::AABase>( seq3));

      // create alignment list
      util::ShPtrList< align::AlignmentInterface< biol::AABase> > alignment_list;
      alignment_list.PushBack( alignment1);
      alignment_list.PushBack( alignment2);
      alignment_list.PushBack( alignment3);

      // construct a pairwise aligner
      util::ShPtr< align::PairwiseAlignerInterface< biol::AABase> >
        sp_pairwise_aligner( new align::AlignerDynamicProgramming< biol::AABase>());
      // create scorecalculate_correlations_sm
      score::AssignmentWithGap< biol::AABase> assign_score
      (
        util::CloneToShPtr( score::AAAssignmentBLOSUM( score::AAAssignmentBLOSUM::e_BLOSUM_45)),
        -1.0, -1.0, -1.0, -1.0
      );
      // test constructor taking pairwise aligner
      align::AlignerProgressive< biol::AABase> aligner( sp_pairwise_aligner, assign_score);
      // create calculate_correlations_sm using default constructor
      contact::CalculateCorrelationsSM calculate_correlations_sm;
      // align all three sequences and test size of alignment
      storage::Pair< align::AlignmentNode< biol::AABase>, double> alignment_pair( aligner.AlignMultiple( alignment_list));
      align::AlignmentNode< biol::AABase> result_alignment( alignment_pair.First());
      // Test out correlation matrix calculator on the short test sequences
      contact::CorrelationMatrix correlation_matrix( calculate_correlations_sm( result_alignment));

    /////////////////
    // data access //
    /////////////////

      // Test GetSize()
      size_t size_chk( 6);
      contact::CorrelationMatrix default_matrix( size_chk);
      BCL_ExampleIndirectCheck( default_matrix.GetSize(), size_chk, "Check Size");

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // Test normal value

      // Test NaN Conversion

    //////////////////////
    // input and output //
    //////////////////////

      // Test write and read
      WriteBCLObject( correlation_matrix);
      // Read in object
      contact::CorrelationMatrix read_object;
      ReadBCLObject( read_object);
      BCL_ExampleIndirectAssert( util::Format()( correlation_matrix), util::Format()( read_object), "Read and Write");

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleContactCorrelationMatrix

  const ExampleClass::EnumType ExampleContactCorrelationMatrix::s_Instance( GetExamples().AddEnum( ExampleContactCorrelationMatrix()));

} // namespace bcl
