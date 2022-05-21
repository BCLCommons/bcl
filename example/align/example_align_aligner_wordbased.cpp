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
#include "align/bcl_align_aligner_wordbased.h"

// includes from bcl - sorted alphabetically
#include "align/bcl_align_handler_pir.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "function/bcl_function_binary_sum.h"
#include "io/bcl_io_file.h"
#include "score/bcl_score_aa_assignment_pam.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_align_aligner_wordbased.cpp
  //!
  //! @author heinzes1
  //! @date 02/03/2011
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAlignAlignerWordbased :
    public ExampleInterface
  {
  public:

    ExampleAlignAlignerWordbased *Clone() const
    {
      return new ExampleAlignAlignerWordbased( *this);
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
      // prepare everything: read in a sequence to test
      io::IFStream read;

      // read seq_a from fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_test1.fasta"));
      util::ShPtr< biol::AASequence> seq_a( new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( read)));
      io::File::CloseClearFStream( read);

      // read seq from fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_test2.fasta"));
      util::ShPtr< biol::AASequence> seq_b( new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( read)));
      io::File::CloseClearFStream( read);

      // create alignment and alignment_list
      util::ShPtr< align::AlignmentInterface< biol::AABase> >
        alignment_a( new align::AlignmentLeaf< biol::AABase>( seq_a)),
        alignment_b( new align::AlignmentLeaf< biol::AABase>( seq_b));
      util::ShPtrList< align::AlignmentInterface< biol::AABase> > alignment_list( 1, alignment_b);

      // create assignment score
      score::AssignmentWithGap< biol::AABase> assign_score
      (
        util::CloneToShPtr( double( 1.0) * score::AAAssignmentPAM( score::AAAssignmentPAM::e_PAM_100)),
        -0.5, -0.5, -0.5, -0.5
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct with assign score
      align::AlignerWordbased< biol::AABase> aligner( assign_score);

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // align alignment_a to all alignments in the list
      storage::List< storage::Pair< align::AlignmentNode< biol::AABase>, double> >
        result_alignments( aligner.AlignPairwise( alignment_a, alignment_list));

      align::HandlerPIR< biol::AABase> handler_pir;
      for
      (
        storage::List< storage::Pair< align::AlignmentNode< biol::AABase>, double> >::const_iterator
          itr( result_alignments.Begin()),
          itr_end( result_alignments.End());
        itr != itr_end;
        ++itr
      )
      {
        util::GetLogger() << "Alignment: score=" << itr->Second() << '\n';
        handler_pir.WriteAlignment( util::GetLogger(), itr->First());
      }

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAlignAlignerWordbased

  const ExampleClass::EnumType ExampleAlignAlignerWordbased::s_Instance
  (
    GetExamples().AddEnum( ExampleAlignAlignerWordbased())
  );

} // namespace bcl
