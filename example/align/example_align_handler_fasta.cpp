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
#include "align/bcl_align_handler_fasta.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_align_handler_fasta.cpp
  //!
  //! @author heinzes1
  //! @date Apr 12, 2012
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAlignHandlerFasta :
    public ExampleInterface
  {
  public:

    ExampleAlignHandlerFasta *Clone() const
    {
      return new ExampleAlignHandlerFasta( *this);
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

      // create alignment, alignment_list and engine
      util::ShPtr< align::AlignmentInterface< biol::AABase> > alignment1( new align::AlignmentLeaf< biol::AABase>( seq1));
      util::ShPtr< align::AlignmentInterface< biol::AABase> > alignment2( new align::AlignmentLeaf< biol::AABase>( seq2));
      align::AlignerMerge< biol::AABase> aligner;
      // align sequences
      align::AlignmentNode< biol::AABase> result_alignment( aligner.AlignPair( alignment1, alignment2).First());

      // create empty alignment to test
      align::AlignmentNode< biol::AABase> empty_alignment( alignment1, alignment2);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create fasta_handler using default constructor
      align::HandlerFasta< biol::AABase> fasta_handler;

      // test Clone()
      util::ShPtr< align::HandlerFasta< biol::AABase> > fasta_handler_cloned( fasta_handler.Clone());

    /////////////////
    // data access //
    /////////////////

      // test GetClassIdentifier()
      BCL_MessageStd( "1: GetClassIdentifier()")
      const std::string class_identifier( "bcl::align::HandlerFasta<bcl::biol::AABase>");
      const std::string returned_class_identifier( fasta_handler.GetClassIdentifier());
      BCL_ExampleIndirectCheck
      (
        returned_class_identifier,
        class_identifier,
        "1: GetClassIdentifier()==" + util::Format()( returned_class_identifier) + "!=" + util::Format()( class_identifier) + "==class_identifier"
      );

    ////////////////
    // operations //
    ////////////////

      // test writing a empty alignment
      BCL_MessageStd( "2: WriteAlignment() with empty alignment");
      fasta_handler.WriteAlignment( util::GetLogger(), empty_alignment);

      // create output_filename and write_stream to output the alignment into a file
      const std::string output_filename( AddExampleOutputPathToFilename( fasta_handler, "alignment.fasta"));
      io::OFStream write_stream;
      // open the file, write the content, close the file
      BCL_ExampleMustOpenOutputFile( write_stream, output_filename);
      fasta_handler.WriteAlignment( write_stream, result_alignment);
      io::File::CloseClearFStream( write_stream);

      // read in the file and compare
      BCL_ExampleMustOpenInputFile( readstream, output_filename);
      align::AlignmentNode< biol::AABase>
        read_alignment( *fasta_handler.ReadAlignment( readstream, biol::AASequence()));
      io::File::CloseClearFStream( readstream);

      std::stringstream result_alignment_stream;
      fasta_handler.WriteAlignment( result_alignment_stream, result_alignment);
      std::stringstream read_alignment_stream;
      fasta_handler.WriteAlignment( read_alignment_stream, read_alignment);
      BCL_MessageStd( "5: WriteAlignment(), ReadAlignment()")
      BCL_ExampleIndirectCheck
      (
        result_alignment_stream.str(),
        read_alignment_stream.str(),
        "5: output of read in alignment is differs from written alignment"
      );

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAlignHandlerFasta

  const ExampleClass::EnumType ExampleAlignHandlerFasta::s_Instance
  (
    GetExamples().AddEnum( ExampleAlignHandlerFasta())
  );

} // namespace bcl
