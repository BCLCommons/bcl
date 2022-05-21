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
#include "align/bcl_align_handler_blocked.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_align_handler_blocked.cpp
  //!
  //! @author heinzes1
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAlignHandlerBlocked :
    public ExampleInterface
  {
  public:

    ExampleAlignHandlerBlocked *Clone() const
    {
      return new ExampleAlignHandlerBlocked( *this);
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

//      // declare to AASequences from AA and AACaCb
//      biol::AASequence seq1;
//      biol::AASequence seq2;
//      score::AssignmentWithGap< biol::AABase> assign
//      (
//        double( 0.6) * score::AAAssignmentBlastProfile()                          +
//        double( 0.2) * score::AAAssignmentSSPrediction( sspred::GetMethods().e_PSIPRED)            +
//        double( 0.1) * score::AAAssignmentPAM( score::AAAssignmentPAM::e_PAM_160) +
//        double( 0.1) * score::AAAssignmentBLOSUM( score::AAAssignmentBLOSUM::e_BLOSUM_45),
//        -0.5, -0.1, -0.2, -0.1
//      );
//
//      // create "dynamic_programming" which will perform the Alignment
//      align::EngineDynamicProgramming< biol::AABase> dynamic_programming( assign);
//
//      // create align::Factory which will handle the amino acid sequences
//      align::Factory< biol::AABase> align_factory;
//
//      // create "read" which will be used to input the sequences
//      io::IFStream read;
//
//      // read seq1 from fasta
//      BCL_ExampleMustOpenInputFile( e_Biology, "1fms_.fasta");
//      seq1.ReadFasta( read);
//      io::File::CloseClearFStream( read);
//
//      // read psipred ss prediction
//      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_.psipred_ss2"));
//      sspred::MethodHandler::ReadPredictionsForAASequence( read, seq1, sspred::GetMethods().e_PSIPRED);
//      io::File::CloseClearFStream( read);
//
//      // read blast profile
//      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_.ascii"));
//      biol::BlastProfileHandler::ReadProfileForAASequence( read, seq1);
//      io::File::CloseClearFStream( read);
//
//      // read seq2 from fasta
//      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1f5mA.fasta"));
//      seq2.ReadFasta( read);
//      io::File::CloseClearFStream( read);
//
//      // read psipred ss prediction
//      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1f5mA.psipred_ss2"));
//      sspred::MethodHandler::ReadPredictionsForAASequence( read, seq2, sspred::GetMethods().e_PSIPRED);
//      io::File::CloseClearFStream( read);
//
//      // read blast profile
//      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1f5mA.ascii"));
//      biol::BlastProfileHandler::ReadProfileForAASequence( read, seq2);
//      io::File::CloseClearFStream( read);
//
//      // align sequences
//      seq1.WriteFasta( util::GetLogger());
//      seq2.WriteFasta( util::GetLogger());
//
////      // create Pair "alignment_and_score" to hold the alignment of "seq1" and "seq2"
////      storage::Pair< align::AlignmentSimple< biol::AABase>, double> alignment_and_score
////      (
////        align_factory.Aligning( dynamic_programming, seq1.GetData(), seq2.GetData())
////      );
////
////      // create AlignmentHandlerBlocked "write_blocked" to output the Alignment of "alignment_and_score"
////      align::HandlerBlocked< biol::AABase> write_blocked( false, 50, assign);
////
////      // write the Alignment of "alignment_and_score" in blocked format with gaps
////      write_blocked.WriteAlignment( util::GetLogger(), alignment_and_score.First());
////
////      // change "write_blocked" to not print gaps
////      write_blocked.SetUnGapped( true);
////
////      // write the Alignment of "alignment_and_score" in blocked format without gaps
////      write_blocked.WriteAlignment( util::GetLogger(), alignment_and_score.First());
////
////      BCL_ExampleCheck
////      (
////        math::EqualWithinTolerance( 4.7275, alignment_and_score.Second()), true
////      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAlignHandlerBlocked

  const ExampleClass::EnumType ExampleAlignHandlerBlocked::s_Instance
  (
    GetExamples().AddEnum( ExampleAlignHandlerBlocked())
  );

} // namespace bcl
