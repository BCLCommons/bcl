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
#include "score/bcl_score_epr_accessibility.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_epr_accessibility.cpp
  //!
  //! @author alexanns
  //! @date Nov 6, 2010
  //! @remarks status empty
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreEPRAccessibility :
    public ExampleInterface
  {
  public:

    ExampleScoreEPRAccessibility *Clone() const
    { return new ExampleScoreEPRAccessibility( *this);}

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
//    /////////////////
//    // preparation //
//    /////////////////
//
//      //instantiate pdb
//      io::IFStream read;
//      const std::string pdbfilename( AddExampleInputPathToFilename( e_Biology, "2LZM.pdb"));
//      //instantiate protein model of chains
//      const assemble::ProteinModel protein_model
//      (
//        pdb::Factory( biol::GetAAClasses().e_AABackBone).ProteinModelFromPDBFilename( pdbfilename)
//      );
//
//    //////////////////////////////////
//    // construction and destruction //
//    //////////////////////////////////
//
//      // construct from scheme and histogram filename
//      const score::EPRAccessibility score_accessibility
//      (
//        "epr_accessibility", score::EPRAccessibility::GetDefaultHistogramFilename()
//      );
//
//      BCL_Example_Check
//      (
//        ExampleClass::ExampleResult::e_Trivial,
//        score_accessibility.GetScheme() == "epr_accessibility",
//        "scheme was not initialize to \"epr_accessibility\" but is: " + score_accessibility.GetScheme()
//      );
//
//      // test copy constructor
//      const score::EPRAccessibility copy( score_accessibility);
//
//      // test clone
//      const util::ShPtr< score::EPRAccessibility> clone( copy.Clone());
//
//      BCL_Example_Check
//      (
//        ExampleClass::ExampleResult::e_Trivial,
//        score_accessibility.GetScheme() == clone->GetScheme()
//        && score_accessibility.GetScheme() == copy.GetScheme(),
//        "cloning or copying did not work since scheme does not match: "
//        " original: " + score_accessibility.GetScheme() +
//        " copy scheme: " + copy.GetScheme() +
//        " cloned " + clone->GetScheme()
//      );
//
//    /////////////////
//    // data access //
//    /////////////////
//
//      // output scheme
//      BCL_Message
//      (
//        util::Message::e_Standard, "scheme for constructed score::EPRAccessibility: " + score_accessibility.GetScheme()
//      );
//
//      // class identifier and static class name
//      BCL_Example_Check
//      (
//        ExampleClass::ExampleResult::e_Trivial,
//           GetStaticClassName< score::EPRAccessibility>() == "bcl::score::EPRAccessibility"
//        && clone->GetClassIdentifier() == GetStaticClassName< score::EPRAccessibility>(),
//        "class identifier returns wrong identifier: " + clone->GetClassIdentifier() + " != " +
//        GetStaticClassName< score::EPRAccessibility>()
//      );
//
//    ///////////////
//    // operators //
//    ///////////////
//
//      // create HandlerAtomDistanceAssigned "handler"
//      restraint::HandlerAccessibilityAA handler;
//
//      // create string "restraint_filename" which has path for example restraint file
//      const std::string restraint_filename
//      (
//        AddExampleInputPathToFilename( e_Biology, "2LZM_accessibility_aa_nc_seq_excl_3.cst")
//      );
//
//      // create stream to "restraint_filename"
//      BCL_ExampleMustOpenInputFile( read, restraint_filename);
//
//      // test CreateRestraints function
//      BCL_MessageStd( "Create Restraints");
//      // create ShPtrVector "restraint" and initialize with the restraints in "restraint_filename"
//      util::ShPtr
//      <
//        util::ShPtrVector
//        <
//          restraint::RestraintInterface
//          <
//            assemble::ProteinModel,
//            storage::Map< restraint::AccessibilityAA::EnvironmentEnum, double>,
//            double,
//            biol::AABase
//          >
//        >
//      > restraints( handler.CreateRestraints( read).Clone()); //< initialize with restraints creatd by "handler"
//
//      io::File::CloseClearFStream( read);
//      util::ShPtr
//      <
//        math::FunctionInterface
//        <
//          restraint::Assignment
//          <
//            storage::Map< restraint::AccessibilityAA::EnvironmentEnum, double>, double, biol::AABase
//          >,
//          double
//        >
//      > sh_ptr_score_accessibility( score_accessibility.Clone());
//      score::Restraint
//      <
//        assemble::ProteinModel, storage::Map< restraint::AccessibilityAA::EnvironmentEnum,
//        double>, double, biol::AABase
//      > restraint_score
//      (
//        restraints,
//        sh_ptr_score_accessibility
//      );
//
////      const double score( restraint_score( protein_model));
////
////      const double correct_score( -4.20641);
////      BCL_MessageStd( "Score is " + util::Format()( score));
////      BCL_Example_Check
////      (
////        ExampleClass::ExampleResult::e_Trivial,
////        math::EqualWithinTolerance( score, correct_score),
////        "correct score is " + util::Format()( correct_score) + " but computed score is " + util::Format()( score)
////      );
//
//      // histogram filename
//      BCL_Message
//      (
//        util::Message::e_Standard,
//        "histogram file used for generating energies: " + score_accessibility.GetHistogramFilename()
//      );
//
//      // GetEnergyFunction
//      for( double slnc_cbnc( -30.0); slnc_cbnc <= 30; slnc_cbnc += 0.1)
//      {
//        BCL_Message
//        (
//          util::Message::e_Debug, "NCsl - NCcb example value : " + util::Format()( slnc_cbnc) + " score : " +
//          util::Format()( score_accessibility.GetEnergyFunction().F( slnc_cbnc))
//        );
//      }
//
//    ////////////////
//    // operations //
//    ////////////////
//
//    //////////////////////
//    // input and output //
//    //////////////////////
//
//      // write score_accessibility to file
//      WriteBCLObject( score_accessibility);
//      // read score back
//      score::EPRAccessibility read_score_accessibility;
//      ReadBCLObject( read_score_accessibility);
//
//      BCL_Example_Check
//      (
//        ExampleClass::ExampleResult::e_Trivial,
//           score_accessibility.GetScheme() == read_score_accessibility.GetScheme()
//        && score_accessibility.GetHistogramFilename() == read_score_accessibility.GetHistogramFilename(),
//        "score::EPRAccessibility was not read properly!"
//      );
//
//    //////////////////////
//    // helper functions //
//    //////////////////////
//
//      // helper functions are private and are implicitly checked

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreEPRAccessibility

  const ExampleClass::EnumType ExampleScoreEPRAccessibility::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreEPRAccessibility())
  );

} // namespace bcl

