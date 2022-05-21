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
#include "linal/bcl_linal_vector_operations.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_linal_vector_operations.cpp
  //!
  //! @author woetzen
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleLinalVectorOperations :
    public ExampleInterface
  {
  public:

    ExampleLinalVectorOperations *Clone() const
    {
      return new ExampleLinalVectorOperations( *this);
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
      // generate two Vector< double>
      double c[] = { double( 4), double( -3), double( 7)};

    //////////////////
    // construction //
    //////////////////

      // default constructed - empty vector
      linal::Vector< double> v0;

      // constructed from length - initialize all values to double (0.0);
      linal::Vector< double> v1( 3);

      // construct from length and pointer to data
      linal::Vector< double> v2( 3, c);

      // construct from length and pointer to data
      linal::Vector< double> v3( 3, c);

    ///////////////////////
    // binary operations //
    ///////////////////////

      // add v2 to v2
      v2 += v3;
      BCL_MessageStd( "this is v2 += v3: " + util::Format()( v2));

      // multiply v3 times 2
      v3 *= double( 2);
      BCL_MessageStd( "this is v3 *= 2: " + util::Format()( v3));

      BCL_Example_Check
      (
        v2 == v3, "adding v2 to v2 and 2*v3 should be identical"
      );

      // divide v2 by 2
      v2 /= double( 2);
      BCL_MessageStd( "this is v2 /= 2: " + util::Format()( v2));

      // subtract v2 from v3
      v3 -= v2;
      BCL_MessageStd( "this is v3 -= v2: " + util::Format()( v3));

      BCL_Example_Check
      (
        v2 == v3, "dividing v2 by two and subtracting v2 from v3 should be identical"
      );

      // add scalar to v2
      v2 += double( 2);
      BCL_MessageStd( "this is v2 += 2: " + util::Format()( v2));

      // add and subtract sclalr to v3
      v3 -= double( 2);
      v3 += double( 4);
      BCL_MessageStd( "this is v3 -= 2 and v3 += 4: " + util::Format()( v3));

      BCL_Example_Check
      (
        v2 == v3, "adding v2 +=2 and subtracting 2 and adding 4 to v3 should be identical"
      );

      // divide by vector
      v3 /= v3;
      BCL_MessageStd( "v3 /= v3: " + util::Format()( v3));
      BCL_Example_Check
      (
        v3 == double( 1), "vector divided by itself should only have 1.0's"
      );

    //////////////////////////////
    // binary logical operators //
    //////////////////////////////

      const bool comp1( v0 == v3);
      BCL_MessageStd( "compare v0 == v3 " + util::Format()( comp1));
      const bool comp2( v0 != v3);
      BCL_MessageStd( "compare v0 != v3 " + util::Format()( comp2));
      const bool comp3( v1 == v3);
      BCL_MessageStd( "compare v1 == v3 " + util::Format()( comp3));
      const bool comp4( v1 != v3);

      v3 = v2;
      BCL_MessageStd( "compare v1 != v3 " + util::Format()( comp4));
      const bool comp5( v2 == v3);
      BCL_MessageStd( "compare v2 == v3 " + util::Format()( comp5));
      const bool comp6( v2 != v3);
      BCL_MessageStd( "compare v2 != v3 " + util::Format()( comp6));
      const bool comp7( v0 == double( 0));
      BCL_MessageStd( "compare v0 == 0.0 " + util::Format()( comp7));
      const bool comp8( v0 != double( 0));
      BCL_MessageStd( "compare v0 != 0.0 " + util::Format()( comp8));
      const bool comp9( v1 == double( 0));
      BCL_MessageStd( "compare v1 == 0.0 " + util::Format()( comp9));
      const bool comp10( v1 != double( 0));
      BCL_MessageStd( "compare v1 != 0.0 " + util::Format()( comp10));
      const bool comp11( v2 == double( 0));
      BCL_MessageStd( "compare v2 == 0.0 " + util::Format()( comp11));
      const bool comp12( v2 != double( 0));
      BCL_MessageStd( "compare v2 != 0.0 " + util::Format()( comp12));
      const bool comp13( double( 0) == v0);
      BCL_MessageStd( "compare 0.0 == v0 " + util::Format()( comp13));
      const bool comp14( double( 0) != v0);
      BCL_MessageStd( "compare 0.0 != v0 " + util::Format()( comp14));
      const bool comp15( double( 0) == v1);
      BCL_MessageStd( "compare 0.0 == v1 " + util::Format()( comp15));
      const bool comp16( double( 0) != v1);
      BCL_MessageStd( "compare 0.0 != v1 " + util::Format()( comp16));
      const bool comp17( double( 0) == v2);
      BCL_MessageStd( "compare 0.0 == v2 " + util::Format()( comp17));
      const bool comp18( double( 0) != v2);
      BCL_MessageStd( "compare 0.0 != v2 " + util::Format()( comp18));

      BCL_Example_Check
      (
        !comp1 && comp2 &&
        !comp3 && comp4 &&
        comp5 && !comp6 &&
        comp7 && !comp8 &&
        comp9 && !comp10 &&
        !comp11 && comp12 &&
        comp13 && !comp14 &&
        comp15 && !comp16 &&
        !comp17 && comp18,
        "one of the binary comparisons did not return the correct result"
      );

    //////////////////////
    // binary operators //
    //////////////////////

      // plus
      linal::Vector< double> v_plus( v2 + v2);
      BCL_MessageStd( "v2 + v2 = " + util::Format()( v_plus));
      BCL_Example_Check
      (
        v_plus == double( 2) * v2, "v2 + v2 should be equal 2 * v2"
      );

      // minus
      linal::Vector< double> v_minus( v2 - v2);
      BCL_MessageStd( "v2 - v2 = " + util::Format()( v_minus));
      BCL_Example_Check
      (
        v_minus == linal::Vector< double>( v2.GetSize(), double( 0)), "v2 + v2 should be equal 2 * v2"
      );

      // divide
      linal::Vector< double> v_divide( v2 / v2);
      BCL_MessageStd( "v2 / v2 = " + util::Format()( v_divide));
      BCL_Example_Check
      (
        v_divide == double( 1), "v2 / v2 should be equal 1.0"
      );

      // multiply (scalar product)
      const double scalar_product_v2( v2 * v2);
      BCL_MessageStd( "scalar product v2 * v2 = " + util::Format()( scalar_product_v2));
      BCL_Example_Check
      (
        math::EqualWithinTolerance( scalar_product_v2, v2.SquareNorm()),
        "scalar product v2 * v2 should be equal to norm of v2"
      );

      // plus value
      linal::Vector< double> v_plus_value_r( v1 + double( 2));
      linal::Vector< double> v_plus_value_l( double( 2) + v1);
      BCL_MessageStd( "v1 + 2.0 = " + util::Format()( v_plus_value_r));
      BCL_Example_Check
      (
        v_plus_value_r == double( 2) && v_plus_value_l == double( 2),
        "adding 2 to zero vector should result vector with twos"
      );

      // minus value
      linal::Vector< double> v_minus_value_r( v1 - double( 2));
      linal::Vector< double> v_minus_value_l( double( 2) - v1);
      BCL_MessageStd( "v1 - 2.0 = " + util::Format()( v_minus_value_r));
      BCL_MessageStd( "2.0 - v1 = " + util::Format()( v_minus_value_l));
      BCL_Example_Check
      (
        v_minus_value_r == double( -2) && v_minus_value_l == double( 2),
        "subtracting 2.0 from zero vector should result vector with -twos or twos when subtracting zero vector from 2.0"
      );

      // divide by scalar
      v1 = double( 4);
      BCL_MessageStd( "new vi is assigned to be 4: " + util::Format()( v1));
      linal::Vector< double> v_divide_scalar_r( v1 / double( 2));
      linal::Vector< double> v_divide_scalar_l( double( 2) / v1);
      BCL_MessageStd( "v1 / 2.0 = " + util::Format()( v_divide_scalar_r));
      BCL_MessageStd( "2.0 / v1 = " + util::Format()( v_divide_scalar_l));
      BCL_Example_Check
      (
        v_divide_scalar_r == double( 2) && v_divide_scalar_l == double( 0.5),
        "dividing 2.0 by vector should result vector with two's, dividing 2.0 with vector of 4's should result vector with 0.5's"
      );

      // multiply with scalar
      linal::Vector< double> v_multiply_scalar_r( v1 * double( 1.5));
      linal::Vector< double> v_multiply_scalar_l( double( 1.5) * v1);
      BCL_MessageStd( "v1 * 1.5 = " + util::Format()( v_multiply_scalar_r));
      BCL_MessageStd( "1.5 * v1 = " + util::Format()( v_multiply_scalar_l));
      BCL_Example_Check
      (
        v_multiply_scalar_r == double( 6) && v_multiply_scalar_l == double( 6),
        "multiply vector with scalar 1.5 should result vector with 6.0's"
      );

      // value to the power of a vector returning a vector
      linal::Vector< double> v_power( double( 2) ^ v1);
      BCL_MessageStd( "2 ^ v1 = " + util::Format()( v_power));
      BCL_Example_Check
      (
        v_power == double( 16) && v_power.GetSize() == 3,
        "2 power 4 results 16 and this in a vector of size 3"
      );

    ///////////////////////
    // general functions //
    ///////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleLinalVectorOperations

  const ExampleClass::EnumType ExampleLinalVectorOperations::s_Instance
  (
    GetExamples().AddEnum( ExampleLinalVectorOperations())
  );

} // namespace bcl
