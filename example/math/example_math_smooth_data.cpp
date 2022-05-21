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
#include "math/bcl_math_smooth_data.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_smooth_data.cpp
  //!
  //! @author butkiem1
  //! @date
  //! @remarks status empty
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathSmoothData :
    public ExampleInterface
  {
  public:

    ExampleMathSmoothData *Clone() const
    { return new ExampleMathSmoothData( *this);}

    static int SmoothVector()
    {
      const util::Format format_vector;
      //construct a test dataset
      const double data_array_vector[] = { 0.0, 1.1, 1.9, 3.1, 3.9, 5.1, 5.9, 6.1, 7.0};
      const linal::Vector< double> test_data( 9, data_array_vector);

      const linal::Vector< double> smoothed_data( math::SmoothData::SmoothVector( test_data, 0.5));
      const linal::Vector< double> smoothed_data_borders( math::SmoothData::SmoothVector( test_data, 0.5, true));

      BCL_MessageStd( "this is the smoothed dataset linear mit 0.5 weight:\n" +
                             format_vector( smoothed_data));
      BCL_MessageStd( "this is the smoothed dataset linear mit 0.5 weight and smoothed borders:\n" +
                             format_vector( smoothed_data_borders));

      return 0;
    }

    static int SmoothMatrix()
    {
      const util::Format format_matrix;

      //construct a test dataset
      const double data_array_matrix[] =
                          { 0.0, 1.0, 2.0, 3.0,
                            1.0, 2.0, 3.0, 4.0,
                            2.0, 3.0, 4.0, 5.0,
                            3.0, 4.0, 5.0, 6.0};

      const linal::Matrix< double> test_data( 4, 4, data_array_matrix);

      const linal::Matrix< double> smoothed_data( math::SmoothData::SmoothMatrix( test_data, 0.5));
      const linal::Matrix< double> smoothed_data_borders( math::SmoothData::SmoothMatrix( test_data, 0.5, true));

      BCL_MessageStd( "this is the smoothed dataset linear mit 0.5 weight:\n" +
                             format_matrix( smoothed_data));
      BCL_MessageStd( "this is the smoothed dataset linear mit 0.5 weight and smoothed borders:\n" +
                             format_matrix( smoothed_data_borders));

      return 0;
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
      return
        SmoothVector() ||
        SmoothMatrix();
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathSmoothData

  const ExampleClass::EnumType ExampleMathSmoothData::s_Instance
  (
    GetExamples().AddEnum( ExampleMathSmoothData())
  );

} // namespace bcl
