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
#include "opencl/bcl_opencl_dataset_min_max.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opencl_dataset_min_max.cpp
  //!
  //! @author woetzen
  //! @date Nov 16, 2010
  //! @remarks status empty
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOpenclDatasetMinMax :
    public ExampleInterface
  {
  public:

    ExampleOpenclDatasetMinMax *Clone() const
    {
      return new ExampleOpenclDatasetMinMax( *this);
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
      // do not try to run opencl commands if no queue was found
      if( !opencl::GetTools().HasCommandQueues())
      {
        return 1;
      }

      static const float data[] =
      {
         5,  6,  7,  8,
         1,  2,  3,  4,
        13, 14, 15, 16,
         9, 10, 11, 12
      };

      const linal::Matrix< float> matrix( 4, 4, data);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      opencl::DataSetMinMax< float> min_max_function( opencl::GetTools().GetFirstCommandQueue());

      const linal::Vector< float> max( min_max_function( matrix));
      BCL_MessageStd( "result max: " + util::Format()( max));

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
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleOpenclDatasetMinMax

  const ExampleClass::EnumType ExampleOpenclDatasetMinMax::s_Instance
  (
    GetExamples().AddEnum( ExampleOpenclDatasetMinMax())
  );

} // namespace bcl
