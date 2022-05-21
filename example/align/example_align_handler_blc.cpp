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
#include "align/bcl_align_handler_blc.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_align_handler_blc.cpp
  //!
  //! @author heinzes1
  //! @date 2010/11/16
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAlignHandlerBLC :
    public ExampleInterface
  {
  public:

    ExampleAlignHandlerBLC *Clone() const
    {
      return new ExampleAlignHandlerBLC( *this);
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

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create blc_handler testing default constructor
      align::HandlerBLC< biol::AABase> blc_handler;

      // test Clone()
      util::ShPtr< align::HandlerBLC< biol::AABase> > blc_handler_cloned( blc_handler.Clone());

    /////////////////
    // data access //
    /////////////////

      // test GetClassIdentifier()
      BCL_MessageStd( "1: GetClassIdentifier()")
      const std::string class_identifier( "bcl::align::HandlerBLC<bcl::biol::AABase>");
      const std::string returned_class_identifier( blc_handler.GetClassIdentifier());
      BCL_ExampleIndirectCheck
      (
        returned_class_identifier,
        class_identifier,
        "1: GetClassIdentifier()==" + util::Format()( returned_class_identifier) + "!=" + util::Format()( class_identifier) + "==class_identifier"
      );

    ////////////////
    // operations //
    ////////////////

      // create output_filename and write_stream to output the alignment into a file
      const std::string output_filename( AddExampleOutputPathToFilename( blc_handler, "alignment.blc"));
      io::OFStream write_stream;
      // open the file, write the content, close the file
      BCL_ExampleMustOpenOutputFile( write_stream, output_filename);
      blc_handler.WriteAlignment( write_stream, result_alignment);
      io::File::CloseClearFStream( write_stream);

      // read in the file and compare
      BCL_ExampleMustOpenInputFile( readstream, output_filename);
      align::AlignmentNode< biol::AABase> read_alignment( *blc_handler.ReadAlignment( readstream, biol::AASequence()));
      io::File::CloseClearFStream( readstream);

      std::stringstream result_alignment_stream;
      blc_handler.WriteAlignment( result_alignment_stream, result_alignment);
      std::stringstream read_alignment_stream;
      blc_handler.WriteAlignment( read_alignment_stream, read_alignment);
      BCL_MessageStd( "2: WriteAlignment(), ReadAlignment()")
      BCL_ExampleIndirectCheck
      (
        result_alignment_stream.str(),
        read_alignment_stream.str(),
        "2: output of read in alignment is differs from written alignment"
      );

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

      // do not test Write() and Read(), they are empty

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAlignHandlerBLC

  const ExampleClass::EnumType ExampleAlignHandlerBLC::s_Instance( GetExamples().AddEnum( ExampleAlignHandlerBLC()));

} // namespace bcl
