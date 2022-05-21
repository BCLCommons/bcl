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
#include "model/bcl_model_data_set_statistics.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_data_set_statistics.cpp
  //!
  //! @author mendenjl
  //! @date Aug 17, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelDataSetStatistics :
    public ExampleInterface
  {
  public:

    ExampleModelDataSetStatistics *Clone() const
    {
      return new ExampleModelDataSetStatistics( *this);
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

      // load in the example score dataset
      const std::string model_directory( AddExampleInputPathToFilename( e_Model, ""));
      util::Implementation< model::RetrieveDataSetBase>
        dataset_retriever( "File(filename=" + model_directory + "example_data_set_score.bcl)");

      model::DataSetStatistics stats;
      stats.AddDataSet( dataset_retriever);

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( stats.GetMinMaxFeatures().GetMin(), linal::MakeVector< float>( 0.0, 0.0));
      BCL_ExampleCheck( stats.GetMinMaxFeatures().GetMax(), linal::MakeVector< float>( 1.0, 2.0));
      BCL_ExampleCheck( stats.GetMinMaxResults().GetMin(), linal::MakeVector< float>( -6.0));
      BCL_ExampleCheck( stats.GetMinMaxResults().GetMax(), linal::MakeVector< float>( 6.0));
      BCL_ExampleCheckWithinAbsTolerance
      (
        stats.GetAveStdFeatures().GetAverage(),
        linal::MakeVector< float>( 0.5, 0.75),
        0.0001
      );
      BCL_ExampleCheckWithinAbsTolerance
      (
        stats.GetAveStdResults().GetAverage(),
        linal::MakeVector< float>( 0.0),
        0.0001
      );
      BCL_ExampleCheckWithinAbsTolerance
      (
        stats.GetAveStdFeatures().GetStandardDeviation(),
        linal::MakeVector< float>( math::Sqrt( 8.4 / 72.0), math::Sqrt( 25.5 / 72.0)),
        0.0001
      );
      BCL_ExampleCheckWithinAbsTolerance
      (
        stats.GetAveStdResults().GetStandardDeviation(),
        linal::MakeVector< float>( math::Sqrt( 271.4333333 / 72.0)),
        0.0001
      );

      BCL_MessageStd( "Stats: " + util::Format()( stats));

      // end
      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelDataSetStatistics

  const ExampleClass::EnumType ExampleModelDataSetStatistics::s_Instance
  (
    GetExamples().AddEnum( ExampleModelDataSetStatistics())
  );

} // namespace bcl
