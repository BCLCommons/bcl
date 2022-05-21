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
#include "align/bcl_align_handler_standard.h"

// includes from bcl - sorted alphabetically
#include "align/bcl_align_aligner_dynamic_programming.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "function/bcl_function_binary_sum.h"
#include "io/bcl_io_file.h"
#include "score/bcl_score_aa_assignment_blosum.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_align_handler_standard.cpp
  //!
  //! @author heinzes1, linders
  //! @date 02/19/2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAlignHandlerStandard :
    public ExampleInterface
  {
  public:

    ExampleAlignHandlerStandard *Clone() const
    {
      return new ExampleAlignHandlerStandard( *this);
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
      // declare AASequences and readstream which will be used to input the sequences
      io::IFStream readstream;

      // read seq1 from fasta
      BCL_ExampleMustOpenInputFile( readstream, AddExampleInputPathToFilename( e_Biology, "1fms_.fasta"));
      util::ShPtr< biol::AASequence> seq1
      (
        new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( readstream))
      );
      io::File::CloseClearFStream( readstream);
      // read seq2 from fasta
      BCL_ExampleMustOpenInputFile( readstream, AddExampleInputPathToFilename( e_Biology, "1f5mA.fasta"));
      util::ShPtr< biol::AASequence> seq2
      (
        new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( readstream))
      );
      io::File::CloseClearFStream( readstream);

      // create score
      score::AssignmentWithGap< biol::AABase> assign_score
      (
        util::CloneToShPtr( double( 1.0) * score::AAAssignmentBLOSUM( score::AAAssignmentBLOSUM::e_BLOSUM_45)),
        -1.0, -1.0, -1.0, -1.0
      );

      // create alignment, alignment_list and engine
      util::ShPtr< align::AlignmentInterface< biol::AABase> > alignment1( new align::AlignmentLeaf< biol::AABase>( seq1));
      util::ShPtr< align::AlignmentInterface< biol::AABase> > alignment2( new align::AlignmentLeaf< biol::AABase>( seq2));
      align::AlignerDynamicProgramming< biol::AABase> aligner( assign_score);
      // align sequences
      align::AlignmentNode< biol::AABase> result_alignment( aligner.AlignPair( alignment1, alignment2).First());

      // create empty alignment to test
      align::AlignmentNode< biol::AABase> empty_alignment( alignment1, alignment2);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test constructor taking assignment scoring function
      align::HandlerStandard< biol::AABase> standard_handler( util::CloneToShPtr( assign_score));

      // test Clone()
      util::ShPtr< align::HandlerStandard< biol::AABase> > standard_handler_cloned( standard_handler.Clone());

      standard_handler_cloned->WriteAlignment( util::GetLogger(), result_alignment);

    /////////////////
    // data access //
    /////////////////

      // test GetClassIdentifier()
      BCL_MessageStd( "1: GetClassIdentifier()")
      BCL_ExampleCheck
      (
        standard_handler.GetClassIdentifier(),
        "bcl::align::HandlerStandard<bcl::biol::AABase>"
      );

    ////////////////
    // operations //
    ////////////////

      // test writing a empty alignment
      BCL_MessageStd( "2: WriteAlignment() with empty alignment");
      standard_handler.WriteAlignment( util::GetLogger(), empty_alignment);

      // create output_filename and write_stream to output the alignment into a file
      const std::string output_filename( AddExampleOutputPathToFilename( standard_handler, "alignment.standard"));
      io::OFStream write_stream;
      // open the file, write the content, close the file
      BCL_ExampleMustOpenOutputFile( write_stream, output_filename);
      standard_handler.WriteAlignment( write_stream, result_alignment);
      io::File::CloseClearFStream( write_stream);

      // read in the file and compare
      BCL_ExampleMustOpenInputFile( readstream, output_filename);
      align::AlignmentNode< biol::AABase> read_alignment( *standard_handler.ReadAlignment( readstream, biol::AASequence()));
      io::File::CloseClearFStream( readstream);

      std::stringstream result_alignment_stream;
      standard_handler.WriteAlignment( result_alignment_stream, result_alignment);
      std::stringstream read_alignment_stream;
      standard_handler.WriteAlignment( read_alignment_stream, read_alignment);
      BCL_MessageStd( "3: WriteAlignment(), ReadAlignment()")
      BCL_ExampleIndirectCheck
      (
        result_alignment_stream.str(),
        read_alignment_stream.str(),
        "3: output of read in alignment is differs from written alignment"
      );

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

      // no test for Write() and Read(), class has no members

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAlignHandlerStandard

  const ExampleClass::EnumType ExampleAlignHandlerStandard::s_Instance
  (
    GetExamples().AddEnum( ExampleAlignHandlerStandard())
  );

} // namespace bcl
