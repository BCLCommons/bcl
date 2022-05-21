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
#include "math/bcl_math_tensor.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_statistics.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_tensor.cpp
  //!
  //! @author butkiem1
  //! @date
  //! @remarks status empty
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathTensor :
    public ExampleInterface
  {
  public:

    ExampleMathTensor *Clone() const
    { return new ExampleMathTensor( *this);}

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
      math::Tensor< double> tensor1( 4, 4, 4, double( 1.0));
      math::Tensor< double> tensor2( 4, 4, 4);

      math::Statistics::SetRand( tensor1.Begin(), tensor1.End(), 0.0, 10.0);
      BCL_MessageStd( "this is a 4 x 4 x 4 tensor with random values");
      BCL_MessageStd( util::Format()( tensor1));

      BCL_MessageStd
      (
        "this is element ( 1, 2, 3) which means ( layer, row, col): " + util::Format()( tensor1( 1, 2, 3))
      );

    //////////////////////
    // input and output //
    //////////////////////

      // write to file
      WriteBCLObject( tensor1);
      ReadBCLObject( tensor2);
      BCL_MessageStd( util::Format()( tensor2));

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathTensor

  const ExampleClass::EnumType ExampleMathTensor::s_Instance
  (
    GetExamples().AddEnum( ExampleMathTensor())
  );

} // namespace bcl
