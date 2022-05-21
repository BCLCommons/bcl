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
#include "score/bcl_score_sas_type.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"
#include "restraint/bcl_restraint_sas_debye.h"
#include "restraint/bcl_restraint_sas_transformation.h"
#include "score/bcl_score_restraint_saxs.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_saxs_type.cpp
  //! @details Example file to test different SAXS scores comparing two SAXS profiles
  //!
  //! @author putnamdk
  //! @date Sept 15, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreSasType :
    public ExampleInterface
  {
  public:

    ExampleScoreSasType *Clone() const
    {
      return new ExampleScoreSasType( *this);
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
      // Get the Saxs Data
      // read in crysol generated data for ubiquitin

      util::ShPtr< restraint::SasScatteringData> sp_exp_data( new restraint::SasScatteringData());
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1ubi_saxs_crysol.gnom"));
      sp_exp_data->ReadFromDataFile( read);
      io::File::CloseClearFStream( read);

      // Get the protein model
      // read in the protein model
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"), biol::GetAAClasses().e_AAComplete)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // constructor taking parameters
      const std::string scheme( "test_saxs_score");

      std::string parameters
      (
        "SasDebye(consider loops=0, analytic=0, excluded volume=1.0, hydration shell=0.0)"
      );

      // use the non-opencl version of the debye formula
      util::Implementation< restraint::SasDebyeInterface> saxs( parameters);
      saxs->SetExperimentalData( sp_exp_data);

      // Compute the raw saxs profile
      restraint::SasExperimentalAndCalculatedData raw_data( saxs->operator()( protein_model));

      //Test the Chi Score

      // Define Data Transformations
      storage::Vector< restraint::SasTransformation::TransformationTypeEnum> transforms;

      transforms.PushBack( restraint::SasTransformation::e_ScaleCalculatedProfile);
      transforms.PushBack( restraint::SasTransformation::e_SetYScale);
      transforms.PushBack( restraint::SasTransformation::e_Log10Profiles);

      // transform the profile into the desired form
      restraint::SasExperimentalAndCalculatedData chi_transformed_data
      (
        restraint::SasTransformation( transforms, false, false, 1.0)( raw_data)
      );

      // calculate the score
      double chi( score::SasType( false, score::SasType::e_chi)( chi_transformed_data));

      //Test the derivative chi score
      transforms.PushBack( restraint::SasTransformation::e_DerivativeProfiles);
      restraint::SasExperimentalAndCalculatedData dchi_transformed_data
      (
        restraint::SasTransformation( transforms, false, false, 1.0)( raw_data)
      );

      double dchi( score::SasType( false, score::SasType::e_chi)( dchi_transformed_data));

      //Test the cumulative integral score
      storage::Vector< restraint::SasTransformation::TransformationTypeEnum> cumulative_integral_transforms;
      cumulative_integral_transforms.PushBack( restraint::SasTransformation::e_Log10Profiles);
      cumulative_integral_transforms.PushBack( restraint::SasTransformation::e_NormalizeData);

      // transform the profile into the desired form
      restraint::SasExperimentalAndCalculatedData cumulative_transformed_data
      (
        restraint::SasTransformation( cumulative_integral_transforms, false, false, 1.0)( raw_data)
      );

      double cumulative( score::SasType( false, score::SasType::e_cumulative)( cumulative_transformed_data));

      //Test the stovgaard score
      double stovgaard( score::SasType( false, score::SasType::e_stovgaard)( raw_data));

    ///////////////
    // operators //
    ///////////////

      // test different scores
      BCL_ExampleCheckWithinTolerance( chi,        0.0604384, 0.001);
      BCL_ExampleCheckWithinTolerance( dchi,       2.7469,   0.001);
      BCL_ExampleCheckWithinTolerance( cumulative, 0.246789,  0.001);
      BCL_ExampleCheckWithinTolerance( stovgaard,  12.4172,   0.001);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintSaxs

  const ExampleClass::EnumType ExampleScoreSasType::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreSasType())
  );

} // namespace bcl
