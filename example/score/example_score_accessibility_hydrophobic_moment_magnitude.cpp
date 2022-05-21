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
#include "score/bcl_score_accessibility_hydrophobic_moment_magnitude.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_vector.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_ifstream.h"
#include "pdb/bcl_pdb_factory.h"
#include "restraint/bcl_restraint_handler_accessibility_aa.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_accessibility_hydrophobic_moment_magnitude.cpp
  //! @brief TODO: add an detailed description for this example
  //!
  //! @author alexanns
  //! @date Feb 7, 2012
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreAccessibilityHydrophobicMomentMagnitude :
    public ExampleInterface
  {
  public:

    ExampleScoreAccessibilityHydrophobicMomentMagnitude *Clone() const
    {
      return new ExampleScoreAccessibilityHydrophobicMomentMagnitude( *this);
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
        const score::AccessibilityHydrophobicMomentMagnitude scorer
        (
          restraint::AccessibilityAA::e_Oxygen, window_sizes
        );
        const double score( scorer( profile_assignment));
        const double correct_score( -1);
        BCL_ExampleCheckWithinTolerance( score, correct_score, 0.0000001);
        BCL_MessageDbg( "score is " + util::Format()( score));

        scorer.WriteDetailedSchemeAndValues( profile_assignment, util::GetLogger());
      }

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    /////////////////
    // data access //
    /////////////////

    ////////////////
    // operations //
    ////////////////

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

  }; //end ExampleScoreAccessibilityHydrophobicMomentMagnitude

  const ExampleClass::EnumType ExampleScoreAccessibilityHydrophobicMomentMagnitude::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreAccessibilityHydrophobicMomentMagnitude())
  );

} // namespace bcl
