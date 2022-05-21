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
#include "model/bcl_model_support_vector_kernel_rbf.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_dataset.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_support_vector_kernel_rbf.cpp
  //!
  //! @author butkiem1
  //! @date Mar 7, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelSupportVectorKernelRbf :
    public ExampleInterface
  {
  public:

    ExampleModelSupportVectorKernelRbf *Clone() const
    {
      return new ExampleModelSupportVectorKernelRbf( *this);
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
      // monitoring dataset
      storage::Vector< storage::VectorND< 2, linal::Vector< float> > > monitoring_data;

      // generate monitoring data for dataset
      for( size_t counter( 0); counter < 10; ++counter)
      {
        monitoring_data.PushBack
        (
          storage::VectorND< 2, linal::Vector< float> >
          (
            linal::MakeVector< float>( counter / 10.0, ( counter + 1) / 10.0, ( counter + 2) / 10.0),
            linal::MakeVector< float>( float( counter % 2) / 10.0)
          )
        );
      }

      descriptor::Dataset monitor_frds( monitoring_data);

      model::FeatureDataSet< float> features( monitor_frds.GetFeaturesPtr()->GetMatrix());

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      model::SupportVectorKernelRBF kernel_default;

      // constructor with parameters
      util::Implementation< model::SupportVectorKernelBase> kernel( std::string( "RBF(gamma=0.33)"));

      // check clone
      util::ShPtr< model::SupportVectorKernelRBF> kernel_clone( kernel_default.Clone());

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

      // check operator
      BCL_Example_Check
      (
        math::EqualWithinTolerance
        (
          kernel->operator ()( features( 0), features( 5)),
          float( 0.78075),
          0.01
        ),
        "1) rbf kernel value is incorrect! should be 0.78075 but is "
        + util::Format()( kernel->operator ()( features( 0), features( 5)))
      );

      // check operator
      BCL_Example_Check
      (
        math::EqualWithinTolerance
        (
          kernel->operator ()( features( 3), features( 7)),
          float( 0.853508),
          0.01
        ),
        "2) rbf kernel value is incorrect! should be 0.853508 but is "
        + util::Format()( kernel->operator ()( features( 3), features( 7)))
      );

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
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelSupportVectorKernelRbf

  const ExampleClass::EnumType ExampleModelSupportVectorKernelRbf::s_Instance
  (
    GetExamples().AddEnum( ExampleModelSupportVectorKernelRbf())
  );

} // namespace bcl
