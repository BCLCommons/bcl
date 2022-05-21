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
#include "contact/bcl_contact_correlation_storage_file.h"

// includes from bcl - sorted alphabetically
#include "align/bcl_align_aligner_progressive.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "contact/bcl_contact_calculate_correlations_sm.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_contact_correlation_storage_file.cpp
  //! @brief TODO: add an detailed description for this example
  //!
  //! @author teixeipl
  //! @date Nov 6, 2011
  //! @remarks status incomplete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleContactCorrelationStorageFile :
    public ExampleInterface
  {
  public:
    //! typedef for better readability
    typedef storage::Pair< align::AlignmentNode< biol::AABase>, contact::CorrelationMatrix> AlignmentMatrixPair;

    ExampleContactCorrelationStorageFile *Clone() const
    {
      return new ExampleContactCorrelationStorageFile( *this);
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
      // Test Constructor
      contact::CorrelationStorageFile storage;

      // Test clone
      contact::CorrelationStorageFile *clone_ptr( storage.Clone());

      AlignmentMatrixPair am_pair;

      std::string initializer( AddExampleOutputPathToFilename( contact::GetNamespaceIdentifier(), "test_write_am_pair_storage_directory"));

    /////////////////
    // data access //
    /////////////////

      // Test GetDirectory, try to attach and if it doesn't exist yet then create it
      if( !storage.Initialize( initializer, contact::CorrelationStorageFile::e_Attach))
      {
        storage.Initialize( initializer, contact::CorrelationStorageFile::e_Create);
      }
      BCL_ExampleIndirectCheck( storage.GetDirectory()->GetPath(), initializer, "Testing GetDirectory");

    ////////////////
    // operations //
    ////////////////
      // Create alignment to create AMPair
      // read seq1 from fasta
      io::IFStream readstream;
      BCL_ExampleMustOpenInputFile( readstream, AddExampleInputPathToFilename( e_Biology, "1EPWtest1.fasta"));
      util::ShPtr< biol::AASequence> seq1
      (
        new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( readstream))
      );
      io::File::CloseClearFStream( readstream);
      seq1->SetFastaHeader( "1EPWtest1");
      // read seq2 from fasta
      BCL_ExampleMustOpenInputFile( readstream, AddExampleInputPathToFilename( e_Biology, "1EPWtest2.fasta"));
      util::ShPtr< biol::AASequence> seq2
      (
        new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( readstream))
      );
      io::File::CloseClearFStream( readstream);
      seq2->SetFastaHeader( "1EPWtest2");
      // read seq3 from fasta
      BCL_ExampleMustOpenInputFile( readstream, AddExampleInputPathToFilename( e_Biology, "1EPWtest3.fasta"));
      util::ShPtr< biol::AASequence> seq3
      (
        new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( readstream))
      );
      io::File::CloseClearFStream( readstream);
      seq3->SetFastaHeader( "1EPWtest3");
      // create alignment, alignment_list and engine
      util::ShPtr< align::AlignmentInterface< biol::AABase> > alignment1( new align::AlignmentLeaf< biol::AABase>( seq1));
      util::ShPtr< align::AlignmentInterface< biol::AABase> > alignment2( new align::AlignmentLeaf< biol::AABase>( seq2));
      util::ShPtr< align::AlignmentInterface< biol::AABase> > alignment3( new align::AlignmentLeaf< biol::AABase>( seq3));
      // create alignment list
      util::ShPtrList< align::AlignmentInterface< biol::AABase> > alignment_list;
      alignment_list.PushBack( alignment1);
      alignment_list.PushBack( alignment2);
      alignment_list.PushBack( alignment3);
      // create scorecalculate_correlations_sm
      score::AssignmentWithGap< biol::AABase> assign_score
      (
        util::CloneToShPtr( score::AAAssignmentBLOSUM( score::AAAssignmentBLOSUM::e_BLOSUM_45)),
        -1.0, -1.0, -1.0, -1.0
      );
      // construct a pairwise aligner
      util::ShPtr< align::PairwiseAlignerInterface< biol::AABase> >
        sp_pairwise_aligner( new align::AlignerDynamicProgramming< biol::AABase>());
      // Construct aligner
      align::AlignerProgressive< biol::AABase> aligner( sp_pairwise_aligner, assign_score);
      // create calculate_correlations_sm using default constructor
      contact::CalculateCorrelationsSM calculate_correlations_sm;
      // align all three sequences and test size of alignment
      storage::Pair< align::AlignmentNode< biol::AABase>, double> alignment_pair( aligner.AlignMultiple( alignment_list));
      align::AlignmentNode< biol::AABase> result_alignment( alignment_pair.First());
      // Calculate correlation matrix
      contact::CorrelationMatrix correlation_matrix( calculate_correlations_sm( result_alignment));

      // Create test AMPair
      AlignmentMatrixPair test_am_pair( result_alignment, correlation_matrix);

      // Test storage of one AMPair
      storage.Store( test_am_pair);

      // Test Storage of AMPair with given key
      storage.Store( test_am_pair, "test_key");

      // Test Storage of an ensemble of AMPairs

      // Test Retrieval of one AMPair

      // Test Retrieval of all AMPairs in source

      // Test Retrieval for a given list of keys

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

      // Test write and read
      BCL_ExampleCheck( TestBCLObjectIOForSymmetry( *clone_ptr, contact::CorrelationStorageFile()), true);
      BCL_ExampleCheck( clone_ptr->GetSize(), 0);

      // Todo: Fix counting subdirectories but none exist, check to see if it should
//      BCL_ExampleCheck( storage.GetSize(), 1);

    //////////////////////
    // helper functions //
    //////////////////////

      // Test filename creator helper function

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleContactCorrelationStorageFile

  const ExampleClass::EnumType ExampleContactCorrelationStorageFile::s_Instance
  (
    GetExamples().AddEnum( ExampleContactCorrelationStorageFile())
  );

} // namespace bcl
