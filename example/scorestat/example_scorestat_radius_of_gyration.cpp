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
#include "scorestat/bcl_scorestat_radius_of_gyration.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "biol/bcl_biol_aa_classes.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_scorestat_radius_of_gyration.cpp
  //! @brief this example tests the implementation of RadiusOfGyration in the scorestat namespace for extracting
  //!        for computing radius of gyration statistics from an ensemble of protein models
  //!
  //! @author lib14
  //! @date February 02, 2015
  //! @remarks status complete
  //! @remarks reviewed by
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScorestatRadiusOfGyration :
    public ExampleInterface
  {

  public:

      //! single instance of this class
      static const ExampleClass::EnumType s_Instance;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

      //! @brief Clone function
      //! @return a pointer to a new ExampleScorestatRadiusOfGyration object
      ExampleScorestatRadiusOfGyration *Clone() const
      {
        return new ExampleScorestatRadiusOfGyration( *this);
      }

  /////////////////
  // data access //
  /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

  ////////////////
  // operations //
  ////////////////

      //! @brief run routine
      //! @detail performs the actual testing

      int Run() const
      {

        // protein ensemble filename
        const std::string protein_ensemble_filename
        (
          AddExampleInputPathToFilename( e_Scorestat, "test_pdb_chains.list")
        );

        // paths to example pdb files
        const std::string example_pdb_paths( AddExampleInputPathToFilename( e_Biology, ""));

        // create a protein ensemble
        const assemble::ProteinEnsemble protein_ensemble
        (
          protein_ensemble_filename, 0, biol::GetAAClasses().e_AAComplete, example_pdb_paths
        );

        // create an instance of RadiusOfGyration
        const scorestat::RadiusOfGyration radius_of_gyration;

        // do the analysis
        const std::string output_prefix( AddExampleOutputPathToFilename( "bcl::scorestat", "received_"));
        radius_of_gyration.WriteAnalysisFile( output_prefix, protein_ensemble);

        // check results
        const std::string written_file( output_prefix + radius_of_gyration.GetOutFilePostfix());
        BCL_MessageDbg( "The analysis file is being written to " + util::Format()( written_file));

        const std::string expected_file
        (
          AddExampleInputPathToFilename( e_Scorestat, "") + "expected_" + radius_of_gyration.GetOutFilePostfix()
        );
        BCL_MessageDbg( "The expected file is " + util::Format()( written_file));
        BCL_ExampleCheck
        (
          io::File::FilesMatchWithinAbsoluteTolerance( expected_file, written_file, 1.0e-5),
          true
        );

        return 0;
      }
  };

    //! single instance of this class
    const ExampleClass::EnumType ExampleScorestatRadiusOfGyration::s_Instance
    (
      GetExamples().AddEnum( ExampleScorestatRadiusOfGyration())
    );
} // namespace bcl
