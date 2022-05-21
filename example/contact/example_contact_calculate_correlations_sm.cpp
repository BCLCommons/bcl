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
#include "contact/bcl_contact_calculate_correlations_sm.h"

// includes from bcl - sorted alphabetically
#include "align/bcl_align_aligner_progressive.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_contact_calculate_correlations_sm.cpp
  //!
  //! @author teixeipl
  //! @date 08/01/11
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleContactCalculateCorrelationsSM :
    public ExampleInterface
  {
  public:

    ExampleContactCalculateCorrelationsSM *Clone() const
    {
      return new ExampleContactCalculateCorrelationsSM( *this);
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
      // declare AASequences and readstream which will be used to input the sequences which will be aligned and then
      // will be used by CalculateCorrelationsSM to create a CorrelationMatrix
      io::IFStream readstream;

      // read seq1 from fasta
      BCL_ExampleMustOpenInputFile( readstream, AddExampleInputPathToFilename( e_Biology, "1fms_.fasta"));
      util::ShPtr< biol::AASequence> long_seq1
      (
        new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( readstream))
      );
      io::File::CloseClearFStream( readstream);
      // read seq2 from fasta
      BCL_ExampleMustOpenInputFile( readstream, AddExampleInputPathToFilename( e_Biology, "1f5mA.fasta"));
      util::ShPtr< biol::AASequence> long_seq2
      (
        new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( readstream))
      );
      io::File::CloseClearFStream( readstream);
      // read seq3 from fasta
      BCL_ExampleMustOpenInputFile( readstream, AddExampleInputPathToFilename( e_Biology, "1IE9.fasta"));
      util::ShPtr< biol::AASequence> long_seq3
      (
        new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( readstream))
      );
      io::File::CloseClearFStream( readstream);

      // create alignment, alignment_list and engine
      util::ShPtr< align::AlignmentInterface< biol::AABase> > long_alignment1( new align::AlignmentLeaf< biol::AABase>( long_seq1));
      util::ShPtr< align::AlignmentInterface< biol::AABase> > long_alignment2( new align::AlignmentLeaf< biol::AABase>( long_seq2));
      util::ShPtr< align::AlignmentInterface< biol::AABase> > long_alignment3( new align::AlignmentLeaf< biol::AABase>( long_seq3));

      // create alignment list
      util::ShPtrList< align::AlignmentInterface< biol::AABase> > long_alignment_list;
      long_alignment_list.PushBack( long_alignment1);
      long_alignment_list.PushBack( long_alignment2);
      long_alignment_list.PushBack( long_alignment3);

      // create scorecalculate_correlations_sm
      score::AssignmentWithGap< biol::AABase> assign_score
      (
        util::CloneToShPtr( score::AAAssignmentBLOSUM( score::AAAssignmentBLOSUM::e_BLOSUM_45)),
        -0.8, -1.4, -1.7, -0.1
      );

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

      //  test sequences should be like these
      //   PVTI-N--*
      //   >1EPW:A|PDBID|CHAIN|SEQUENCE
      //   PVTIAN--*
      //   >1EPW:A|PDBID|CHAIN|SEQUENCE
      //   PRQIANAA*
      //   =std=bcl=> Calculated correlation matrix is :        bcl::linal::Matrix<double>

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct a pairwise aligner
      util::ShPtr< align::PairwiseAlignerInterface< biol::AABase> >
        sp_pairwise_aligner( new align::AlignerDynamicProgramming< biol::AABase>());

      // test constructor taking pairwise aligner
      align::AlignerProgressive< biol::AABase> aligner( sp_pairwise_aligner, assign_score);

      // create calculate_correlations_sm using default constructor
      contact::CalculateCorrelationsSM calculate_correlations_sm;

      // create pir_handler using default constructor
      align::HandlerPIR< biol::AABase> pir_handler;

    /////////////////
    // data access //
    /////////////////

    ////////////////
    // operations //
    ////////////////

      // align all three sequences and test size of alignment
      storage::Pair< align::AlignmentNode< biol::AABase>, double> alignment_pair( aligner.AlignMultiple( alignment_list));
      align::AlignmentNode< biol::AABase> result_alignment( alignment_pair.First());

      // Test out correlation matrix calculator on the short test sequences
      contact::CorrelationMatrix correlation_matrix( calculate_correlations_sm( result_alignment));
      contact::CorrelationMatrix original_correlation_matrix( calculate_correlations_sm( result_alignment));   // TODO: REMOVE THIS LINE

//      // Read the previously created short correlation matrix in
//      storage::SymmetricMatrix previous_short_matrix( ReadBCLObject( original_correlation_matrix));
//
//      // Compare the currently calculated short test sequence correlation matrix against the previous result
//      BCL_ExampleIndirectAssert( previous_short_matrix == correlation_matrix, true, "Compare to Previous Correlation Matrix");

      // Calculate and print out short correlation matrix using these alignments
      if( util::GetMessenger().GetCurrentMessageLevel() >= util::Message::e_Verbose)
      {
        // Output correlation matrix created
        // TODO: This is likely redundant, remove if and change BCL_Message
        // level to Verbose, right?
        BCL_MessageStd( +"Calculated correlation matrix is :" + util::Format().W( 8)( correlation_matrix));
      }

      // align all three sequences and test size of alignment
      storage::Pair< align::AlignmentNode< biol::AABase>, double> long_alignment_pair( aligner.AlignMultiple( long_alignment_list));
      align::AlignmentNode< biol::AABase> long_result_alignment( long_alignment_pair.First());

      // Test out correlation matrix calculator on the long test sequences
      storage::SymmetricMatrix< double> long_correlation_matrix( calculate_correlations_sm( long_result_alignment));
//      storage::SymmetricMatrix< double> original_long_correlation_matrix( calculate_correlations_sm( long_result_alignment)); // REMOVE THIS LINE

      // write the object for archiving the "correct" matrix, this line should remain commented out after first write
//      WriteBCLObject( original_long_correlation_matrix);  TODO: RE-ADD THIS
//
//      // Read the previously created long correlation matrix in
//      storage::SymmetricMatrix previous_long_matrix( ReadBCLObject( original_long_correlation_matrix));
//
//      // Compare the currently calculated long test sequence correlation matrix against the previous result
//      BCL_ExampleIndirectAssert( previous_long_matrix == long_correlation_matrix, true, "Compare to Previous Long Correlation Matrix");

      // Calculate and print out long correlation matrix using these alignments
      if( util::GetMessenger().GetCurrentMessageLevel() >= util::Message::e_Verbose)
      {
        // Output correlation matrix created
        BCL_MessageStd( +"Calculated correlation matrix is :" + util::Format().W( 8)( long_correlation_matrix));
      }

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

      // test Write() and Read()
      // Write out object
      WriteBCLObject( calculate_correlations_sm);
      // Read in object
      contact::CalculateCorrelationsSM read_object_sm;
      ReadBCLObject( read_object_sm);
      BCL_ExampleIndirectAssert
      (
        util::Format()( calculate_correlations_sm),
        util::Format()( read_object_sm),
        "Read and Write"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;
  };

  const ExampleClass::EnumType ExampleContactCalculateCorrelationsSM::s_Instance( GetExamples().AddEnum( ExampleContactCalculateCorrelationsSM()));

} // namespace bcl
