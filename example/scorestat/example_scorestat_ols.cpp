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
#include "scorestat/bcl_scorestat_ols.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "biol/bcl_biol_aa_classes.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_scorestat_ols.cpp
  //! @brief this example tests the implementation of OLS for extracting neighbor count statistics from protein
  //!        ensemble
  //!
  //! @author lib14, mendenjl
  //! @date May 18, 2015
  //! @remarks status complete
  //! @remarks reviewed by
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScorestatOLS :
    public ExampleInterface
  {

  public:

      //! single instance of this class
      static const ExampleClass::EnumType s_Instance;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return a pointer to a new ExampleScorestatOLS object
    ExampleScorestatOLS *Clone() const
    {
      return new ExampleScorestatOLS( *this);
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

      // read in protein ensemble file
      const std::string protein_ensemble_filename
      (
        AddExampleInputPathToFilename( e_Scorestat, "test_pdb_chains.list")
      );

      // paths for example pdb files in test_pdb_chains.list
      const std::string example_pdb_paths( AddExampleInputPathToFilename( e_Biology, ""));

      // create an object of ProteinEnsemble
      const assemble::ProteinEnsemble protein_ensemble
      (
        protein_ensemble_filename, 0, biol::GetAAClasses().e_AAComplete, example_pdb_paths
      );

      // print the size of this protein ensemble
      BCL_MessageStd( "The size of this protein ensemble is " + util::Format()( protein_ensemble.GetSize()));

      // create a vector of OLS
      storage::Vector< scorestat::OLS> ols_vec;
      scorestat::OLS ols_default( false);
      ols_default.AssertRead( util::ObjectDataLabel( "(chain_option=One)"));
      ols_vec.PushBack( ols_default);
      ols_default.AssertRead( util::ObjectDataLabel( "(chain_option=All)"));
      ols_vec.PushBack( ols_default);
      scorestat::OLS ols_env( true);
      ols_env.AssertRead( util::ObjectDataLabel( "(chain_option=One)"));
      ols_vec.PushBack( ols_env);
      ols_env.AssertRead( util::ObjectDataLabel( "(chain_option=All)"));
      ols_vec.PushBack( ols_env);

      // run the analysis and write the results
      for
      (
        storage::Vector< scorestat::OLS>::const_iterator
          ols_itr( ols_vec.Begin()), ols_itr_end( ols_vec.End());
        ols_itr != ols_itr_end;
        ++ols_itr
      )
      {
        const std::string output_prefix( AddExampleOutputPathToFilename( "bcl::scorestat", "received_"));
        ols_itr->WriteAnalysisFile( output_prefix, protein_ensemble);

        // check results
        const std::string written_file( output_prefix + ols_itr->GetOutFilePostfix());

        const std::string expected_file
        (
          AddExampleInputPathToFilename( e_Scorestat, "") +
          "expected_" + ols_itr->GetOutFilePostfix()
        );

        BCL_ExampleCheck( io::File::FilesMatchWithinAbsoluteTolerance( written_file, expected_file, 15.0), true);
      }

      return 0;
    }
  }; // class ExampleScorestatOLS

    //! single instance of this class
    const ExampleClass::EnumType ExampleScorestatOLS::s_Instance
    (
      GetExamples().AddEnum( ExampleScorestatOLS())
    );

} // namespace bcl
