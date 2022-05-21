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
#include "linal/bcl_linal_matrix_inversion_cholesky.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_operations.h"
#include "linal/bcl_linal_vector_reference.h"
#include "math/bcl_math_statistics.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_linal_matrix_inversion_cholesky.cpp
  //!
  //! @author mendenjl
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleLinalMatrixInversionCholesky :
    public ExampleInterface
  {
  public:

    ExampleLinalMatrixInversionCholesky *Clone() const
    { return new ExampleLinalMatrixInversionCholesky( *this);}

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
      const util::Format format_matrix( util::Format().W( 7).FFP( 3));

      // instantiate two matrices - 4x4 and 2x4 with 2.5 as value
      linal::Matrix< double> ma( 8, 8);
      // set all elements of matrix a to random elements
      math::Statistics::SetRand( ma.Begin(), ma.End(), 0.0, 10.0);
      ma = ma * ma.Transposed();
      // output of ma
      BCL_MessageStd( " this is a random 8x8 symmetric matrix ma\n" + format_matrix( ma));

      // invert and output t1
      BCL_MessageStd( "mai is the inverse of matrix ma");
      linal::MatrixInversionCholesky< double> mai( ma);
      BCL_ExampleCheck( linal::MatrixInversionCholesky< double>( ma).IsDefined(), true);
      BCL_ExampleCheckWithinAbsTolerance
      (
        ( mai.ComputeInverse() * ma).AsVector().SquareNorm(),
        double( ma.GetNumberRows()),
        0.00001
      );

      BCL_MessageStd( format_matrix( mai.ComputeInverse()));

      // proof that t2 is indeed t1
      BCL_MessageStd( "as we can easily prove: ma * mai = unity matrix");
      BCL_MessageStd( format_matrix( ma * mai.ComputeInverse()));

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // end ExampleLinalMatrixInversionCholesky

  const ExampleClass::EnumType ExampleLinalMatrixInversionCholesky::s_Instance
  (
    GetExamples().AddEnum( ExampleLinalMatrixInversionCholesky())
  );

} // namespace bcl
