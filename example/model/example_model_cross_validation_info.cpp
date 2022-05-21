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
#include "model/bcl_model_cross_validation_info.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_cross_validation_info.cpp
  //!
  //! @author mendenjl
  //! @date May 16, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelCrossValidationInfo :
    public ExampleInterface
  {
  public:

    ExampleModelCrossValidationInfo *Clone() const
    {
      return new ExampleModelCrossValidationInfo( *this);
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

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      model::CrossValidationInfo def_cross_validation_info;

      // file containing a data set
      const std::string data_set_filename( AddExampleInputPathToFilename( e_Model, "example_data_set.bcl"));

      util::ObjectDataLabel monitoring_label
      (
        "File(filename=" + data_set_filename + ", number chunks=4, chunks=\"[0]+[2]\")"
      );
      util::ObjectDataLabel independent_label
      (
        "File(filename=" + data_set_filename + ", number chunks=4, chunks=\"[1]\")"
      );
      util::ObjectDataLabel training_label
      (
        "File(filename=" + data_set_filename + ", number chunks=4, chunks=\"[3]\")"
      );
      util::ObjectDataLabel objective_label( "RMSD");
      util::ObjectDataLabel iterate_label
      (
        "LinearRegression(objective function=" + objective_label.ToString() + ")"
      );
      const float result( 0.235);
      const opti::ImprovementTypeEnum improvement_type
      (
        opti::e_SmallerEqualIsBetter
      );
      const std::string independent_predictions
      (
        io::File::MakeAbsolutePath( AddExampleInputPathToFilename( e_Model, "prediction_merge_1.txt"))
      );

      // constructor from members
      model::CrossValidationInfo cross_validation_info
      (
        result,
        improvement_type,
        independent_label,
        monitoring_label,
        training_label,
        objective_label,
        iterate_label,
        independent_predictions,
        model::FeatureLabelSet(),
        1
      );

    /////////////////
    // data access //
    /////////////////

      // Result should be undefined for defaultly constructed info
      BCL_ExampleCheck( util::IsDefined( def_cross_validation_info.GetResult()), false);

      // Check each member and access function
      BCL_ExampleCheck( cross_validation_info.GetResult(), result);
      BCL_ExampleCheck( cross_validation_info.GetMonitoringDatasetRetriever(), monitoring_label);
      BCL_ExampleCheck( cross_validation_info.GetIndependentDatasetRetriever(), independent_label);
      BCL_ExampleCheck( cross_validation_info.GetTrainingDatasetRetriever(), training_label);
      BCL_ExampleCheck( cross_validation_info.GetImprovementType(), improvement_type);
      BCL_ExampleCheck( cross_validation_info.GetObjective(), objective_label);
      BCL_ExampleCheck( cross_validation_info.GetIterate(), iterate_label);
      BCL_ExampleCheck( cross_validation_info.GetIndependentPredictionsFilename(), independent_predictions);

      // check that the default constructed version can be initialized with the label from the fully-initialized verison
      BCL_MessageStd( "Cross validation info label: " + util::Format()( cross_validation_info.GetLabel()));
      def_cross_validation_info.TryRead( cross_validation_info.GetLabel(), util::GetLogger());
      monitoring_label.SetName( "monitoring");
      independent_label.SetName( "independent");
      training_label.SetName( "training");
      objective_label.SetName( "objective");
      iterate_label.SetName( "iterate");
      BCL_ExampleIndirectCheckWithinAbsTolerance( def_cross_validation_info.GetResult(), result, 0.001, "GetLabel/TryRead");
      BCL_ExampleIndirectCheck( def_cross_validation_info.GetMonitoringDatasetRetriever(), monitoring_label, "GetLabel/TryRead");
      BCL_ExampleIndirectCheck( def_cross_validation_info.GetIndependentDatasetRetriever(), independent_label, "GetLabel/TryRead");
      BCL_ExampleIndirectCheck( def_cross_validation_info.GetTrainingDatasetRetriever(), training_label, "GetLabel/TryRead");
      BCL_ExampleIndirectCheck( def_cross_validation_info.GetImprovementType(), improvement_type, "GetLabel/TryRead");
      BCL_ExampleIndirectCheck( def_cross_validation_info.GetObjective(), objective_label, "GetLabel/TryRead");
      BCL_ExampleIndirectCheck( def_cross_validation_info.GetIterate(), iterate_label, "GetLabel/TryRead");
      BCL_ExampleIndirectCheck( def_cross_validation_info.GetIndependentPredictionsFilename(), independent_predictions, "GetLabel/TryRead");

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelCrossValidationInfo

  const ExampleClass::EnumType ExampleModelCrossValidationInfo::s_Instance
  (
    GetExamples().AddEnum( ExampleModelCrossValidationInfo())
  );

} // namespace bcl
