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
#include "restraint/bcl_restraint_analyze_sas.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "io/bcl_io_file.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_analyze_saxs.cpp
  //! @brief Tests the functionality of restraint analyze saxs.
  //!
  //! @author putnamdk
  //! @date Sept 24, 2014
  //! @remarks status incomplete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintAnalyzeSaxs :
    public ExampleInterface
  {
  public:

      ExampleRestraintAnalyzeSaxs *Clone() const
    {
      return new ExampleRestraintAnalyzeSaxs( *this);
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

      // Read in the protein model
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ENH.pdb"), biol::GetAAClasses().e_AAComplete)
      );

      std::string exp_data( AddExampleInputPathToFilename( e_Biology, "1ENH00.saxs"));

      // Read in SASA Data
      std::string sasa_data( AddExampleInputPathToFilename( e_Biology, "1ENH.area"));

      // Create protein Ensemble
      assemble::ProteinEnsemble protein_ensemble;

      // Add protein_model to ensemble
      protein_ensemble.InsertElement( util::CloneToShPtr( protein_model));

      // Test chi function without error using transformations
      {
        std::string parameters
        (
          "AnalyzeSas( c1=1, c2=0, experimental_profile=" + exp_data +
          ", default_search_grid=false, scoring_function=chi,use_errors=0, approximate_side_chains=false, "
          "approximate_loops=false, use_sans=0, transformations( Scale, SetYMax, Log10, Derivative ), y_max=1)"
        );

        restraint::AnalyzeSas saxs;
        BCL_ExampleCheck( saxs.TryRead( parameters, util::GetLogger()), true);

        BCL_ExampleIndirectCheckWithinAbsTolerance
        (
          util::ConvertStringToNumericalValue< double>( saxs( protein_ensemble)),
          0.222487,
          0.001,
          "Analyze Saxs with " + parameters
        );
      }

      // Test parameter optimization  if default search grid is used, a chi score of 0.103277 is found
      // With the stepsizes altered, a less optimal score is found, but the function still works in principle
      // The search window is relaxed to save time in the example check procedure
      {
        std::string parameters
        (
          "AnalyzeSas( sasa_profile=" + sasa_data + " , experimental_profile="
          + exp_data + ", scoring_function=chi, optimize_hydration_parameters=true, c1=1, c2=0, "
          "default_search_grid=false, c1_min=0.8, c1_max=1.20, c2_min=0, c2_max=4, c1_stepsize=0.2, c2_stepsize=2, "
          "use_errors=0, approximate_side_chains=false, approximate_loops=false, "
          "use_sans=0, transformations( Scale, SetYMax, Log10, Derivative), y_max=1)"
        );

        restraint::AnalyzeSas saxs;
        BCL_ExampleCheck( saxs.TryRead( parameters, util::GetLogger()), true);

        BCL_ExampleIndirectCheckWithinAbsTolerance
        (
          util::ConvertStringToNumericalValue< double>( saxs( protein_ensemble)),
          0.282548,
          0.001,
          "Analyze Saxs with " + parameters
        );

      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintAnalyzeSaxs

  const ExampleClass::EnumType ExampleRestraintAnalyzeSaxs::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintAnalyzeSaxs())
  );

} // namespace bcl
