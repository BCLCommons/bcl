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
#include "score/bcl_score_accessibility_hydrophobic_moment.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_vector.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "restraint/bcl_restraint_handler_accessibility_aa.h"
#include "util/bcl_util_colors.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_accessibility_hydrophobic_moment.cpp
  //!
  //! @author alexanns
  //! @date Apr 8, 2011
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreAccessibilityHydrophobicMoment :
    public ExampleInterface
  {
  public:

    ExampleScoreAccessibilityHydrophobicMoment *Clone() const
    {
      return new ExampleScoreAccessibilityHydrophobicMoment( *this);
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
      const std::string cst_filename( AddExampleInputPathToFilename( e_Biology, "2lzm_cst.access_bcl"));

      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, cst_filename);

      const util::ShPtr< assemble::AAExposureInterface> exposure_calculator( new assemble::AANeighborVector());
      restraint::HandlerAccessibilityAA handler( exposure_calculator);

      const restraint::AccessibilityProfile profile( handler.ReadRestraints( read));

      const std::string pdb_file( AddExampleInputPathToFilename( e_Biology, "2LZM_idealized.pdb"));

      pdb::Factory factory;

      const assemble::ProteinModel model( factory.ProteinModelFromPDBFilename( pdb_file));

      const restraint::AccessibilityProfileAssignment assignments( profile.GenerateAssignment( model));

      const restraint::AccessibilityProfileAssignment profile_assignment( profile.GenerateAssignment( model));

      // window size 1
      {
        const size_t window_size( 1);
        storage::Map< biol::SSType, size_t> window_sizes
        (
          storage::Map< biol::SSType, size_t>::Create
          (
            std::pair< biol::SSType, size_t>( biol::GetSSTypes().HELIX, window_size)
          )
        );
        const score::AccessibilityHydrophobicMoment scorer
        (
          restraint::AccessibilityAA::e_Oxygen, window_sizes
        );
        const double score( scorer( profile_assignment));
        const double correct_score( -1);
        BCL_ExampleCheck( score, correct_score);
        BCL_MessageDbg( "score is " + util::Format()( score));
        io::OFStream write;
        std::string windows_out_filename( AddExampleOutputPathToFilename( scorer, "ShowHydrophobicMomentWindows_a.py"));
        BCL_ExampleMustOpenOutputFile( write, windows_out_filename);
        BCL_MessageDbg( "printing to file " + windows_out_filename);
        const restraint::AccessibilityAA::EnvironmentType environment( restraint::AccessibilityAA::e_Oxygen);
        storage::List< score::AccessibilityHydrophobicMoment::Window> window_list
        (
          score::AccessibilityHydrophobicMoment::CalculateHydrophobicMomentWindows
          (
            profile_assignment.GetSSEAssignments().Begin()->second, window_size, environment
          )
        );
        score::AccessibilityHydrophobicMoment::ShowHydrophobicMomentWindows
        (
          window_list,
          write,
          *profile_assignment.GetSSEAssignments().Begin()->first,
          std::string( "sse_a"),
          util::GetColors().e_Cyan
        );
        io::File::CloseClearFStream( write);
        BCL_ExampleCheck( io::File::FilesMatch( windows_out_filename, windows_out_filename + ".correct"), true);
      }

      // window size 9
      {
        const size_t window_size( 9);
        storage::Map< biol::SSType, size_t> window_sizes
        (
          storage::Map< biol::SSType, size_t>::Create
          (
            std::pair< biol::SSType, size_t>( biol::GetSSTypes().HELIX, window_size)
          )
        );
        const score::AccessibilityHydrophobicMoment scorer
        (
          restraint::AccessibilityAA::e_Oxygen, window_sizes
        );
        const double score( scorer( profile_assignment));
        const double correct_score( -1);
        BCL_ExampleCheckWithinTolerance( score, correct_score, 0.0000001);
        BCL_MessageDbg( "score is " + util::Format()( score));
        io::OFStream write;
        std::string windows_out_filename( AddExampleOutputPathToFilename( scorer, "ShowHydrophobicMomentWindows_b.py"));
        BCL_ExampleMustOpenOutputFile( write, windows_out_filename);
        BCL_MessageDbg( "printing to file " + windows_out_filename);
        const restraint::AccessibilityAA::EnvironmentType environment( restraint::AccessibilityAA::e_Oxygen);
        storage::List< score::AccessibilityHydrophobicMoment::Window> window_list
        (
          score::AccessibilityHydrophobicMoment::CalculateHydrophobicMomentWindows
          (
            profile_assignment.GetSSEAssignments().Begin()->second, window_size, environment
          )
        );
        score::AccessibilityHydrophobicMoment::ShowHydrophobicMomentWindows
        (
          window_list,
          write,
          *profile_assignment.GetSSEAssignments().Begin()->first,
          std::string( "sse_a"),
          util::GetColors().e_Cyan
        );
        io::File::CloseClearFStream( write);
        BCL_ExampleCheck( io::File::FilesMatch( windows_out_filename, windows_out_filename + ".correct"), true);
      }

      // window size 4
      {
        const size_t window_size( 4);
        storage::Map< biol::SSType, size_t> window_sizes
        (
          storage::Map< biol::SSType, size_t>::Create
          (
            std::pair< biol::SSType, size_t>( biol::GetSSTypes().HELIX, window_size)
          )
        );
        const score::AccessibilityHydrophobicMoment scorer
        (
          restraint::AccessibilityAA::e_Oxygen, window_sizes
        );
        const double score( scorer( profile_assignment));
        const double correct_score( -1);
        BCL_ExampleCheckWithinTolerance( score, correct_score, 0.0000001);
        BCL_MessageDbg( "score is " + util::Format()( score));
        io::OFStream write;
        std::string windows_out_filename( AddExampleOutputPathToFilename( scorer, "ShowHydrophobicMomentWindows_c.py"));
        BCL_ExampleMustOpenOutputFile( write, windows_out_filename);
        BCL_MessageDbg( "printing to file " + windows_out_filename);
        const restraint::AccessibilityAA::EnvironmentType environment( restraint::AccessibilityAA::e_Oxygen);
        storage::List< score::AccessibilityHydrophobicMoment::Window> window_list
        (
          score::AccessibilityHydrophobicMoment::CalculateHydrophobicMomentWindows
          (
            profile_assignment.GetSSEAssignments().Begin()->second, window_size, environment
          )
        );
        score::AccessibilityHydrophobicMoment::ShowHydrophobicMomentWindows
        (
          window_list,
          write,
          *profile_assignment.GetSSEAssignments().Begin()->first,
          std::string( "sse_a"),
          util::GetColors().e_Cyan
        );
        io::File::CloseClearFStream( write);
        BCL_ExampleCheck( io::File::FilesMatch( windows_out_filename, windows_out_filename + ".correct"), true);
      }

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

  }; //end ExampleScoreAccessibilityHydrophobicMoment

  const ExampleClass::EnumType ExampleScoreAccessibilityHydrophobicMoment::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreAccessibilityHydrophobicMoment())
  );

} // namespace bcl
