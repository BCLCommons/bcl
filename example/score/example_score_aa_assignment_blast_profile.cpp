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
#include "score/bcl_score_aa_assignment_blast_profile.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "biol/bcl_biol_blast_profile_handler.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_aa_assignment_blast_profile.cpp
  //!
  //! @author karakam
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreAAAssignmentBlastProfile :
    public ExampleInterface
  {
  public:

    ExampleScoreAAAssignmentBlastProfile *Clone() const
    { return new ExampleScoreAAAssignmentBlastProfile( *this);}

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
      // declare to AASequences from AA and AACaCb
      util::ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> >
        score_assignment( new score::AAAssignmentBlastProfile());

      io::IFStream read;

      // read seq1 from fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_.fasta"));
      biol::AASequence seq1( biol::AASequenceFactory::BuildSequenceFromFASTA( read));
      io::File::CloseClearFStream( read);

      // read blast profile
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_.ascii"));
      biol::BlastProfileHandler::ReadProfileForAASequence( read, seq1);
      io::File::CloseClearFStream( read);

      // read seq2 from fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1f5mA.fasta"));
      biol::AASequence seq2( biol::AASequenceFactory::BuildSequenceFromFASTA( read));
      io::File::CloseClearFStream( read);

      // read blast profile
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1f5mA.ascii"));
      biol::BlastProfileHandler::ReadProfileForAASequence( read, seq2);
      io::File::CloseClearFStream( read);

      BCL_MessageStd( "blast profile of first aa of 1fms_ " + util::Format()( seq1.GetFirstAA()->GetBlastProfile()));
      BCL_MessageStd( "blast profile of first aa of 1f5mA " + util::Format()( seq2.GetFirstAA()->GetBlastProfile()));
      BCL_MessageStd( "scoring alignment of first amino acids of 1fms_ and 1f5mA");

      const double assign_score( score_assignment->operator()( *seq1.GetFirstAA(), *seq2.GetFirstAA()));
      const double expected_assign_score( 0.0357336);
      BCL_MessageStd( "assignment score is " + util::Format()( assign_score));
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_assign_score, assign_score), "Assignment score does not return expected score of " + util::Format()( expected_assign_score)
      );

    //////////////////////
    // input and output //
    //////////////////////

      // test Write
      WriteBCLObject( *score_assignment);
      // create instance of class "AAAssignmentBlastProfile" and read from file
      score::AAAssignmentBlastProfile read_score;
      ReadBCLObject( read_score);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreAAAssignmentBlastProfile

  const ExampleClass::EnumType ExampleScoreAAAssignmentBlastProfile::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreAAAssignmentBlastProfile())
  );

} // namespace bcl
