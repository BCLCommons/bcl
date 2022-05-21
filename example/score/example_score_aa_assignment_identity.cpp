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
#include "score/bcl_score_aa_assignment_identity.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_aa_assignment_identity.cpp
  //!
  //! @author riddeljs
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreAAAssignmentIdentity :
    public ExampleInterface
  {
  public:

    ExampleScoreAAAssignmentIdentity *Clone() const
    { return new ExampleScoreAAAssignmentIdentity( *this);}

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

      // default constructor
      score::AAAssignmentIdentity aaa_assign_id1;

      // clone
      util::ShPtr< score::AAAssignmentIdentity> sp_aaa_assign_id2( aaa_assign_id1.Clone());

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

      // only operator() needs to be tested, so set up the parameters

      // declare to AASequences from AA and AACaCb
      util::ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> >
        score_assignment_identity( new score::AAAssignmentIdentity());

      io::IFStream read;

      // read seq1 from fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_.fasta"));
      biol::AASequence seq1( biol::AASequenceFactory::BuildSequenceFromFASTA( read));
      io::File::CloseClearFStream( read);

      // read seq2 from fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1f5mA.fasta"));
      biol::AASequence seq2( biol::AASequenceFactory::BuildSequenceFromFASTA( read));
      io::File::CloseClearFStream( read);

      BCL_MessageStd
      (
        "aatype of first aa of 1fms_ and first aa of 1f5mA " +
        seq1.GetFirstAA()->GetType()->GetThreeLetterCode() + " " + seq2.GetFirstAA()->GetType()->GetThreeLetterCode()
      );

      BCL_MessageStd
      (
        "aatype of second aa of 1fms_ and first aa of 1f5mA " +
        seq1.GetData()( 1)->GetType()->GetThreeLetterCode() + " " + seq2.GetData()( 0)->GetType()->GetThreeLetterCode()
      );

      // now call operator()
      const double assign_score_11( score_assignment_identity->operator()( *seq1.GetFirstAA(), *seq2.GetFirstAA()));
      const double assign_score_12( score_assignment_identity->operator()( *seq1.GetData()( 1), *seq2.GetData()( 0)));
      const double expected_assign_score_11( 0.0);
      const double expected_assign_score_12( 1.0);
      BCL_MessageStd( "assignment score 11 is " + util::Format()( assign_score_11));
      BCL_MessageStd( "assignment score 12 is " + util::Format()( assign_score_12));
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_assign_score_11, assign_score_11),
        "Assignment score 11 does not return expected score of " + util::Format()( assign_score_11)
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_assign_score_12, assign_score_12),
        "Assignment score 12 does not return expected score of " + util::Format()( assign_score_12)
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // test Write
      WriteBCLObject( aaa_assign_id1);
      // create instance of class "AAAssignmentIdentity" and read from file
      score::AAAssignmentIdentity read_score;
      ReadBCLObject( read_score);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreAAAssignmentIdentity

  const ExampleClass::EnumType ExampleScoreAAAssignmentIdentity::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreAAAssignmentIdentity())
  );

} // namespace bcl
