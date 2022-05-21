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
#include "score/bcl_score_aa_assignment_blosum.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_aa_assignment_blosum.cpp
  //!
  //! @author riddeljs
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreAAAssignmentBlosum :
    public ExampleInterface
  {
  public:

    ExampleScoreAAAssignmentBlosum *Clone() const
    { return new ExampleScoreAAAssignmentBlosum( *this);}

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

      // default constructor from a blosum_table
      score::AAAssignmentBLOSUM AAA_assign_blosum1( score::AAAssignmentBLOSUM::e_BLOSUM_90);

      // clone
      util::ShPtr< score::AAAssignmentBLOSUM> sp_AAA_assign_blosum2
      (
        AAA_assign_blosum1.Clone()
      );

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ///////////////
    // operations //
    ///////////////

      // test Probability, but first we must set up the parameters.

      util::ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> >
        score_assignment_blosum_90( new score::AAAssignmentBLOSUM( score::AAAssignmentBLOSUM::e_BLOSUM_90));
      util::ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> >
        score_assignment_blosum_80( new score::AAAssignmentBLOSUM( score::AAAssignmentBLOSUM::e_BLOSUM_80));
      util::ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> >
        score_assignment_blosum_62( new score::AAAssignmentBLOSUM( score::AAAssignmentBLOSUM::e_BLOSUM_62));
      util::ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> >
        score_assignment_blosum_45( new score::AAAssignmentBLOSUM( score::AAAssignmentBLOSUM::e_BLOSUM_45));

      io::IFStream read;

      // read seq1 from fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_.fasta"));
      biol::AASequence seq1( biol::AASequenceFactory::BuildSequenceFromFASTA( read));
      io::File::CloseClearFStream( read);

      // read seq2 from fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1f5mA.fasta"));
      biol::AASequence seq2( biol::AASequenceFactory::BuildSequenceFromFASTA( read));
      io::File::CloseClearFStream( read);

      // now verify that probability
      BCL_MessageStd
      (
        "blosum 90 matrix entry of first aa of 1fms_ and first aa of 1f5mA " +
        util::Format()( score::AAAssignmentBLOSUM::Probability( score::AAAssignmentBLOSUM::e_BLOSUM_90, seq1.GetFirstAA()->GetType(), seq2.GetFirstAA()->GetType()))
      );
      BCL_MessageStd
      (
        "blosum 80 matrix entry of first aa of 1fms_ and first aa of 1f5mA " +
        util::Format()( score::AAAssignmentBLOSUM::Probability( score::AAAssignmentBLOSUM::e_BLOSUM_80, seq1.GetFirstAA()->GetType(), seq2.GetFirstAA()->GetType()))
      );
      BCL_MessageStd
      (
        "blosum 62 matrix entry of first aa of 1fms_ and first aa of 1f5mA " +
        util::Format()( score::AAAssignmentBLOSUM::Probability( score::AAAssignmentBLOSUM::e_BLOSUM_62, seq1.GetFirstAA()->GetType(), seq2.GetFirstAA()->GetType()))
      );
      BCL_MessageStd
      (
        "blosum 45 matrix entry of first aa of 1fms_ and first aa of 1f5mA " +
        util::Format()( score::AAAssignmentBLOSUM::Probability( score::AAAssignmentBLOSUM::e_BLOSUM_45, seq1.GetFirstAA()->GetType(), seq2.GetFirstAA()->GetType()))
      );

      // test operator()
      const double assign_score_90( score_assignment_blosum_90->operator()( *seq1.GetFirstAA(), *seq2.GetFirstAA()));
      const double assign_score_80( score_assignment_blosum_80->operator()( *seq1.GetFirstAA(), *seq2.GetFirstAA()));
      const double assign_score_62( score_assignment_blosum_62->operator()( *seq1.GetFirstAA(), *seq2.GetFirstAA()));
      const double assign_score_45( score_assignment_blosum_45->operator()( *seq1.GetFirstAA(), *seq2.GetFirstAA()));
      const double expected_assign_score_90( -0.1);
      const double expected_assign_score_80( -0.1);
      const double expected_assign_score_62( 0.0);
      const double expected_assign_score_45( 0.0);
      BCL_MessageStd( "assignment score 90 is " + util::Format()( assign_score_90));
      BCL_MessageStd( "assignment score 80 is " + util::Format()( assign_score_80));
      BCL_MessageStd( "assignment score 62 is " + util::Format()( assign_score_62));
      BCL_MessageStd( "assignment score 45 is " + util::Format()( assign_score_45));
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_assign_score_90, assign_score_90),
        "Assignment score 90 does not return expected score of " + util::Format()( assign_score_90)
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_assign_score_80, assign_score_80),
        "Assignment score 80 does not return expected score of " + util::Format()( assign_score_80)
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_assign_score_62, assign_score_62),
        "Assignment score 62 does not return expected score of " + util::Format()( assign_score_62)
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_assign_score_45, assign_score_45),
        "Assignment score 45 does not return expected score of " + util::Format()( assign_score_45)
      );

    //////////////////////
    // input and output //
    //////////////////////

      // test Write
      WriteBCLObject( AAA_assign_blosum1);
      // create instance of class "AAAssignmentBlastProfile" and read from file
      score::AAAssignmentBLOSUM read_score;
      ReadBCLObject( read_score);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreAAAssignmentBlosum

  const ExampleClass::EnumType ExampleScoreAAAssignmentBlosum::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreAAAssignmentBlosum())
  );

} // namespace bcl
