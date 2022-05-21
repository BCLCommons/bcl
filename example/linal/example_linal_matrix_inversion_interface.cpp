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
#include "linal/bcl_linal_matrix_inversion_interface.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_operations.h"
#include "linal/bcl_linal_vector_reference.h"
#include "math/bcl_math_range.h"
#include "math/bcl_math_statistics.h"
#include "random/bcl_random_distribution_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_linal_matrix_inversion_interface.cpp
  //!
  //! @author mendenjl
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleLinalMatrixInversionInterface :
    public ExampleInterface
  {
  public:

    ExampleLinalMatrixInversionInterface *Clone() const
    { return new ExampleLinalMatrixInversionInterface( *this);}

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
      const int dim( 8);
      linal::Matrix< double> ma( dim, dim);
      for( int i( 0); i < dim; ++i)
      {
        ma( i, i) = random::GetGlobalRandom().Double( math::Range< double>( 0.0, 10.0));
        if( i)
        {
          ma( i, i - 1) = random::GetGlobalRandom().Double( math::Range< double>( 0.0, 10.0));
        }
        if( i + 1 < dim)
        {
          ma( i, i + 1) = random::GetGlobalRandom().Double( math::Range< double>( 0.0, 10.0));
        }
      }
      BCL_ExampleCheck( ma.IsTriDiagonal(), true);
      // output of ma
      BCL_MessageStd( " this is a random 8x8 tridiagonal matrix ma\n" + format_matrix( ma));

      // create storage for the inverse
      linal::Matrix< double> inv_storage( dim, dim);

      // invert and output t1
      BCL_MessageStd( "mai is the inverse of matrix ma");
      BCL_ExampleCheck
      (
        linal::MatrixInversionInterface< double>::TryInvertTridiagonalMatrix( inv_storage, ma),
        true
      );
      BCL_ExampleCheckWithinAbsTolerance
      (
        ( ma * inv_storage).AsVector().SquareNorm(),
        double( ma.GetNumberRows()),
        0.00001
      );

      BCL_MessageStd( format_matrix( inv_storage));

      // proof that t2 is indeed t1
      BCL_MessageStd( "as we can easily prove: ma * mai = unity matrix");
      BCL_MessageStd( format_matrix( ma * inv_storage));
      BCL_ExampleCheckWithinAbsTolerance
      (
        linal::Matrix< double>( ma).Inverse(),
        inv_storage,
        0.00001
      );
      linal::Vector< double> random_vector( 8);
      random_vector.SetRand( 0.0, 10.0);
      BCL_ExampleCheckWithinAbsTolerance
      (
        ma.Inverse() * random_vector,
        linal::MatrixInversionInterface< double>::SolveTridiagonalMatrix( ma, random_vector),
        0.00001
      );
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // end ExampleLinalMatrixInversionInterface

  const ExampleClass::EnumType ExampleLinalMatrixInversionInterface::s_Instance
  (
    GetExamples().AddEnum( ExampleLinalMatrixInversionInterface())
  );

} // namespace bcl
