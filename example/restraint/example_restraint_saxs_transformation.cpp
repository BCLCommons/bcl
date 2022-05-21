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
#include "restraint/bcl_restraint_sas_transformation.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"
#include "restraint/bcl_restraint_sas_debye.h"
#include "score/bcl_score_sas_type.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_saxs_transformation.cpp
  //! @details Example file to test different SAXS transformations
  //!
  //! @author putnamdk
  //! @date Sept 17, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintSasTransformation :
    public ExampleInterface
  {
  public:

      ExampleRestraintSasTransformation *Clone() const
    {
      return new ExampleRestraintSasTransformation( *this);
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
      // Get the Experimental and Calculated Data from a file
      util::ShPtr< restraint::SasExperimentalAndCalculatedData>
        sp_raw_data( new restraint::SasExperimentalAndCalculatedData());

      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1ubi_exp_and_cal.data"));
      sp_raw_data->ReadFromDataFile( read);
      io::File::CloseClearFStream( read);

      //Test the scale calculated profile transformation
      storage::Vector< restraint::SasTransformation::TransformationTypeEnum> scaled_transform;
      scaled_transform.PushBack( restraint::SasTransformation::e_ScaleCalculatedProfile);

      // transform the raw data into the desired form
      restraint::SasExperimentalAndCalculatedData scaled_data
      (
        restraint::SasTransformation( scaled_transform, false, false, 1.0)( *sp_raw_data)
      );

      double experimental_scaled_point( scaled_data.GetExperimentalData().GetScatteringLocation( 5).GetIntensity());
      double experimental_error_point( scaled_data.GetExperimentalData().GetScatteringLocation( 5).GetError());
      double calculated_scaled_point( scaled_data.GetCalculatedData().GetScatteringLocation( 5).GetIntensity());

      BCL_ExampleCheckWithinTolerance( experimental_scaled_point, 1.877e+07,   0.001);
      BCL_ExampleCheckWithinTolerance( experimental_error_point,  722.2,       0.001);
      BCL_ExampleCheckWithinTolerance( calculated_scaled_point,   1.88738e+07, 0.001);

      // transform the raw data into the desired form using errors
      restraint::SasExperimentalAndCalculatedData scaled_data_with_error
      (
        restraint::SasTransformation( scaled_transform, false, true, 1.0)( *sp_raw_data)
      );

      double experimental_scaled_point_with_error
      (
        scaled_data_with_error.GetExperimentalData().GetScatteringLocation( 5).GetIntensity()
      );

      double experimental_error_point_with_error
      (
        scaled_data_with_error.GetExperimentalData().GetScatteringLocation( 5).GetError()
      );

      double calculated_scaled_point_with_error
      (
        scaled_data_with_error.GetCalculatedData().GetScatteringLocation( 5).GetIntensity()
      );

      BCL_ExampleCheckWithinTolerance( experimental_scaled_point_with_error, 1.877e+07,   0.001);
      BCL_ExampleCheckWithinTolerance( experimental_error_point_with_error,  722.2,       0.001);
      BCL_ExampleCheckWithinTolerance( calculated_scaled_point_with_error,   1.87983e+07, 0.001);

      //Test the Normalize Data transformation
      storage::Vector< restraint::SasTransformation::TransformationTypeEnum> normalize_data_transform;
      normalize_data_transform.PushBack( restraint::SasTransformation::e_NormalizeData);

      // transform the raw data into the desired form
      restraint::SasExperimentalAndCalculatedData normalize_data
      (
        restraint::SasTransformation( normalize_data_transform, false, false, 1.0)( *sp_raw_data)
      );

      double experimental_normalized_point( normalize_data.GetExperimentalData().GetScatteringLocation( 5).GetIntensity());
      double experimental_normalized_error_point( normalize_data.GetExperimentalData().GetScatteringLocation( 5).GetError());
      double calculated_normalized_point( normalize_data.GetCalculatedData().GetScatteringLocation( 5).GetIntensity());

      BCL_ExampleCheckWithinTolerance( experimental_normalized_point,        0.287177, 0.001);
      BCL_ExampleCheckWithinTolerance( experimental_normalized_error_point,  1.10495e-05,    0.001);
      BCL_ExampleCheckWithinTolerance( calculated_normalized_point,          0.288816, 0.001);

      // transform the raw data into the desired form using errors
      restraint::SasExperimentalAndCalculatedData normalized_data_with_error
      (
        restraint::SasTransformation( normalize_data_transform, false, true, 1.0)( *sp_raw_data)
      );

      double experimental_normalized_point_with_error
      (
        normalized_data_with_error.GetExperimentalData().GetScatteringLocation( 5).GetIntensity()
      );

      double experimental_normalized_error_point_with_error
      (
        normalized_data_with_error.GetExperimentalData().GetScatteringLocation( 5).GetError()
      );

      double calculated_normalized_point_with_error
      (
        normalized_data_with_error.GetCalculatedData().GetScatteringLocation( 5).GetIntensity()
      );

      BCL_ExampleCheckWithinTolerance( experimental_normalized_point_with_error,       0.287177, 0.001);
      BCL_ExampleCheckWithinTolerance( experimental_normalized_error_point_with_error, 1.10495e-05,    0.001);
      BCL_ExampleCheckWithinTolerance( calculated_normalized_point_with_error,         0.288816, 0.001);

      //Test the Set Y Scale Transformation
      storage::Vector< restraint::SasTransformation::TransformationTypeEnum> normalize_dataset_transform;
      normalize_dataset_transform.PushBack( restraint::SasTransformation::e_SetYScale);

      // transform the raw data into the desired form
      restraint::SasExperimentalAndCalculatedData normalize_dataset
      (
        restraint::SasTransformation( normalize_dataset_transform, false, false, 1.0)( *sp_raw_data)
      );

      double experimental_normalized_set_point( normalize_dataset.GetExperimentalData().GetScatteringLocation( 5).GetIntensity());
      double experimental_normalized_set_error_point( normalize_dataset.GetExperimentalData().GetScatteringLocation( 5).GetError());
      double calculated_normalized_set_point( normalize_dataset.GetCalculatedData().GetScatteringLocation( 5).GetIntensity());

      BCL_ExampleCheckWithinTolerance( experimental_normalized_set_point,         0.89084, 0.001);
      BCL_ExampleCheckWithinTolerance( experimental_normalized_set_error_point,   3.42762e-05, 0.001);
      BCL_ExampleCheckWithinTolerance( calculated_normalized_set_point,           0.039545, 0.001);

      // transform the raw data into the desired form using errors
      restraint::SasExperimentalAndCalculatedData normalized_set_data_with_error
      (
        restraint::SasTransformation( normalize_dataset_transform, false, true, 1.0)( *sp_raw_data)
      );

      double experimental_normalized_set_point_with_error
      (
        normalized_set_data_with_error.GetExperimentalData().GetScatteringLocation( 5).GetIntensity()
      );

      double experimental_normalized_set_error_point_with_error
      (
        normalized_set_data_with_error.GetExperimentalData().GetScatteringLocation( 5).GetError()
      );

      double calculated_normalized_set_point_with_error
      (
        normalized_set_data_with_error.GetCalculatedData().GetScatteringLocation( 5).GetIntensity()
      );

      BCL_ExampleCheckWithinTolerance( experimental_normalized_set_point_with_error,       0.89084,     0.001);
      BCL_ExampleCheckWithinTolerance( experimental_normalized_set_error_point_with_error, 3.42762e-05, 0.001);
      BCL_ExampleCheckWithinTolerance( calculated_normalized_set_point_with_error,         0.039545,    0.001);

      //Test the Log10 transformation
      storage::Vector< restraint::SasTransformation::TransformationTypeEnum> log10_transform;
      log10_transform.PushBack( restraint::SasTransformation::e_Log10Profiles);

      // transform the raw data into the desired form
      restraint::SasExperimentalAndCalculatedData log10_data
      (
        restraint::SasTransformation( log10_transform, false, false, 1.0)( *sp_raw_data)
      );

      double experimental_log10_point( log10_data.GetExperimentalData().GetScatteringLocation( 5).GetIntensity());
      double experimental_log_10_error_point( log10_data.GetExperimentalData().GetScatteringLocation( 5).GetError());
      double calculated_log10_point( log10_data.GetCalculatedData().GetScatteringLocation( 5).GetIntensity());

      BCL_ExampleCheckWithinTolerance( experimental_log10_point,         7.27346, 0.001);
      BCL_ExampleCheckWithinTolerance( experimental_log_10_error_point,  1.671e-05, 0.001);
      BCL_ExampleCheckWithinTolerance( calculated_log10_point,           5.92076, 0.001);

      // transform the raw data into the desired form using errors
      restraint::SasExperimentalAndCalculatedData log10_data_with_error
      (
        restraint::SasTransformation( log10_transform, false, true, 1.0)( *sp_raw_data)
      );

      double experimental_log10_point_with_error
      (
        log10_data_with_error.GetExperimentalData().GetScatteringLocation( 5).GetIntensity()
      );

      double experimental_log10_error_point_with_error
      (
        log10_data_with_error.GetExperimentalData().GetScatteringLocation( 5).GetError()
      );

      double calculated_log10_point_with_error
      (
        log10_data_with_error.GetCalculatedData().GetScatteringLocation( 5).GetIntensity()
      );

      BCL_ExampleCheckWithinTolerance( experimental_log10_point_with_error,       7.27346, 0.001);
      BCL_ExampleCheckWithinTolerance( experimental_log10_error_point_with_error, 1.671e-05, 0.001);
      BCL_ExampleCheckWithinTolerance( calculated_log10_point_with_error,         5.92076, 0.001);

      //Test the Derivative transformation
      storage::Vector< restraint::SasTransformation::TransformationTypeEnum> derivative_transform;
      derivative_transform.PushBack( restraint::SasTransformation::e_DerivativeProfiles);

      // transform the raw data into the desired form
      restraint::SasExperimentalAndCalculatedData derivative_data
      (
        restraint::SasTransformation( derivative_transform, false, false, 1.0)( *sp_raw_data)
      );

      double experimental_derivative_point( derivative_data.GetExperimentalData().GetScatteringLocation( 5).GetIntensity());
      double experimental_derivative_error_point( derivative_data.GetExperimentalData().GetScatteringLocation( 5).GetError());
      double calculated_derivative_point( derivative_data.GetCalculatedData().GetScatteringLocation( 5).GetIntensity());

      BCL_ExampleCheckWithinTolerance( experimental_derivative_point,        -8.7e+07, 0.001);
      BCL_ExampleCheckWithinTolerance( experimental_derivative_error_point,  722.2, 0.001);
      BCL_ExampleCheckWithinTolerance( calculated_derivative_point,          -3.98425e+06, 0.001);

      // transform the raw data into the desired form using errors
      restraint::SasExperimentalAndCalculatedData derivative_data_with_error
      (
        restraint::SasTransformation( derivative_transform, false, true, 1.0)( *sp_raw_data)
      );

      double experimental_derivative_point_with_error
      (
        derivative_data_with_error.GetExperimentalData().GetScatteringLocation( 5).GetIntensity()
      );

      double experimental_derivative_error_point_with_error
      (
        derivative_data_with_error.GetExperimentalData().GetScatteringLocation( 5).GetError()
      );

      double calculated_derivative_point_with_error
      (
        derivative_data_with_error.GetCalculatedData().GetScatteringLocation( 5).GetIntensity()
      );

      BCL_ExampleCheckWithinTolerance( experimental_derivative_point_with_error,       -8.7e+07, 0.001);
      BCL_ExampleCheckWithinTolerance( experimental_derivative_error_point_with_error, 722.2, 0.001);
      BCL_ExampleCheckWithinTolerance( calculated_derivative_point_with_error,         -3.98425e+06, 0.001);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintSaxs

  const ExampleClass::EnumType ExampleRestraintSasTransformation::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintSasTransformation())
  );

} // namespace bcl
