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
#include "score/bcl_score_accessibility.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_accessibility.cpp
  //!
  //! @author alexanns
  //! @date Apr 15, 2011
  //! @remarks status empty
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreAccessibility :
    public ExampleInterface
  {
  public:

    ExampleScoreAccessibility *Clone() const
    {
      return new ExampleScoreAccessibility( *this);
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
//      const std::string dir("t4l");
//      // const std::string filename( "/home/alexanns/accessibility/bcl_accessibility_restraint/crystallin/biomol/accessibility_epr.bcl_cst");
//      const std::string filename( "/home/alexanns/accessibility/bcl_accessibility_restraint/"+dir+"/accessibility_epr.bcl_cst");
//
//      io::IFStream read;
//
//      BCL_ExampleMustOpenInputFile( read, filename);
//
//      const restraint::AccessibilityProfile profile
//      (
//        restraint::HandlerAccessibilityAA
//        (
//          score::AANeighborVector::GetDefaultMinimalSequenceSeparation(),
//          score::AANeighborVector::GetDefaultThresholdLowHigh()
//        ).CreateRestraints( read)
//      );
//
//      util::ShPtr< restraint::AccessibilityProfile> access_data( profile.Clone());
//
//      storage::Map< biol::SSType, size_t> ss_types
//      (
//        storage::Map< biol::SSType, size_t>::Create
//        (
//          std::pair< biol::SSType, size_t>( biol::GetSSTypes().COIL, 999),
//          std::pair< biol::SSType, size_t>( biol::GetSSTypes().HELIX, 0),
//          std::pair< biol::SSType, size_t>( biol::GetSSTypes().STRAND, 0)
//        )
//      );
//
//      assemble::ProteinModel model
//      (
//        pdb::Factory( biol::GetAAClasses().e_AAComplete).ProteinModelFromPDBFilename
//          (
//            "/home/alexanns/accessibility/bcl_accessibility_restraint/"+dir+"/2LZM.pdb", ss_types
//          )
//      );
//
//      util::ShPtr< assemble::ProteinModelData> pmd( new assemble::ProteinModelData());
//
//      pmd->Insert( assemble::ProteinModelData::e_AccessibilityEPR, access_data);
//
//      model.SetProteinModelData( pmd);
//
//      const restraint::AccessibilityProfileAssignment assignment( profile.GenerateAssignment( model));
//      BCL_MessageDbg( "sse assignments size " + util::Format()( assignment.GetSSEAssignments().GetSize()));
//
//      score::AccessibilityDistribution score( restraint::AccessibilityAA::e_NiEDDA);
//
//      score::Accessibility model_score
//      (
//        util::ShPtr< score::AccessibilityDistribution>( score.Clone()),
//        assemble::ProteinModelData::e_AccessibilityEPR
//      );
//
//      model_score( model);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreAccessibility

  const ExampleClass::EnumType ExampleScoreAccessibility::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreAccessibility())
  );

} // namespace bcl
