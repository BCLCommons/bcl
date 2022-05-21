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
#include "math/bcl_math.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math.cpp
  //!
  //! @author mueller, putnamdk
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by putnamdk on  Aug 17, 2013
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMath :
    public ExampleInterface
  {
  public:

    ExampleMath *Clone() const
    {
      return new ExampleMath( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // check absolute for different data types
      BCL_ExampleCheck( math::Absolute( int( -1)), 1);
      BCL_ExampleCheck( math::Absolute( int( 1)), 1);
      BCL_ExampleCheck( math::Absolute( size_t( 1)), size_t( 1));
      BCL_ExampleCheck( math::Absolute( float( -1.5)), float( 1.5));
      BCL_ExampleCheck( math::Absolute( float( 1.5)), float( 1.5));
      BCL_ExampleCheck( math::Absolute( double( -1.5)), double( 1.5));
      BCL_ExampleCheck( math::Absolute( double( 1.5)), double( 1.5));
      BCL_ExampleCheck
      (
        math::EqualWithinTolerance
        (
          math::Absolute( std::complex< double>( double( 1.5), double( 1.5))).real(),
          double( math::Sqrt( double( 4.5)))
        ),
        true
      );

      const double values[] =
      {
        -180.00, -172.00, -157.00, -151.00, -145.00, -130.00, -125.00, -114.00, -102.00,
         -90.00,  -81.00,  -76.00,  -60.00,  -53.00,  -47.00,  -30.00,  -27.00,  -15.00,
           0.00,    7.00,   19.00,   25.00,   40.00,   50.00,   59.00,   68.00,   81.00,
          86.00,   95.00,  110.00,  119.00,  121.00,  139.00,  150.00,  160.00,  170.00
      };

      const double abs_values[] =
      {
         180.00,  172.00,  157.00,  151.00,  145.00,  130.00,  125.00,  114.00,  102.00,
          90.00,   81.00,   76.00,   60.00,   53.00,   47.00,   30.00,   27.00,   15.00,
           0.00,    7.00,   19.00,   25.00,   40.00,   50.00,   59.00,   68.00,   81.00,
          86.00,   95.00,  110.00,  119.00,  121.00,  139.00,  150.00,  160.00,  170.00
      };

      storage::Vector< double> test_values( 36, values);
      storage::Vector< double> correct_values( 36, abs_values);

      math::Absolute( test_values);

      BCL_ExampleIndirectCheck
      (
        test_values,
        correct_values,
        "test container version of absolute values"
      );

      std::pair< double, double> pair( -0.5234, -0.6326);
      std::pair< double, double> abs_pair( 0.5234, 0.6326);

      BCL_ExampleCheck( math::Absolute( pair).first, abs_pair.first);
      BCL_ExampleCheck( math::Absolute( pair).second, abs_pair.second);

      int int_number( 5.0);
      float nonint_number( 4.2);

      // check Square Function for integer type
      BCL_ExampleCheckWithinTolerance( math::Sqr( int_number), 25.0, 0.000001);

      // check Square Function for non-integer type
      BCL_ExampleCheckWithinTolerance( math::Sqr( nonint_number), 17.64, 0.000001);

      // check Square Root Function for integer type
      BCL_ExampleCheckWithinTolerance( math::Sqrt( int_number), 2.236068, 0.000001);

      // check Square Root Function for non-integer type
      BCL_ExampleCheckWithinTolerance( math::Sqrt( nonint_number), 2.049390, 0.000001);

      // check Power Function
      BCL_ExampleCheckWithinTolerance( math::Pow( nonint_number, nonint_number), 414.616918, 0.000001);

      // check Size_t implementation of Pow

      size_t value( 2);
      size_t exponent( 7);

      BCL_ExampleCheckWithinTolerance( math::Pow( value, exponent), 128, 0.0001);

      // Check Sign Functions
      int b( -4.0);
      int a( 3.0);

      BCL_ExampleCheck( math::Sign( a, b), -3.0);
      BCL_ExampleCheck( math::Sign( a, a),  3.0);
      BCL_ExampleCheck( math::Sign( b, a),  4.0);
      BCL_ExampleCheck( math::Sign( b, b), -4.0);

      // Check pythagoreans' theorem
      BCL_ExampleCheckWithinTolerance( math::Pythag( a, b), 5.0, 0.0001);

      // Check relational operators for complex numbers

      std::complex< double> number_one( double( 1.5), double( 1.5));
      std::complex< double> number_two( double( 4.7), double( 3.1));

      BCL_ExampleCheck( math::operator<  ( number_one, number_two), true);
      BCL_ExampleCheck( math::operator<  ( number_two, number_one), false);
      BCL_ExampleCheck( math::operator<= ( number_one, number_one), true);
      BCL_ExampleCheck( math::operator>= ( number_two, number_two), true);
      BCL_ExampleCheck( math::operator>  ( number_two, number_one), true);
      BCL_ExampleCheck( math::operator>  ( number_one, number_two), false);

      // check if number is power of two
      BCL_ExampleCheck( math::IsPowerOfTwo(   0), false);
      BCL_ExampleCheck( math::IsPowerOfTwo(   1), true);
      BCL_ExampleCheck( math::IsPowerOfTwo(   2), true);
      BCL_ExampleCheck( math::IsPowerOfTwo(  36), false);
      BCL_ExampleCheck( math::IsPowerOfTwo( 512), true);

      // test of ARGUEMENT_A is equal within absolute Tolerance, in this case the error margin is epsilon
      BCL_ExampleCheck( math::EqualWithinAbsoluteTolerance( double( 2.615), double( 2.62)), false);

      // test of the same values from above with the error margin set to 0.1, in this case it should evaluate to true
      BCL_ExampleCheck( math::EqualWithinAbsoluteTolerance( double( 2.615), double( 2.62), 0.1), true);

      BCL_ExampleCheck( math::EqualWithinTolerance( double( 2.615), double( 2.62)), false);

      const double tolerance_values1[] = { -180.00, -172.00, -157.00, -151.00, -145.00, -130.00, -125.00, -114.00 };
      const double tolerance_values2[] = {  180.00,  172.00,  157.00,  151.00,  145.00,  130.00,  125.00,  114.00 };

      storage::Vector< double> tol1( 8, tolerance_values1);
      storage::Vector< double> tol2( 8, tolerance_values2);

      bool result( math::EqualWithinTolerance( tol1, tol2));

      BCL_ExampleIndirectCheck( result, false, "test EqualWithinTolerance");

      // test WeightBetweenZeroandPi_ThreeSections
      BCL_ExampleCheckWithinTolerance( math::WeightBetweenZeroAndPi_ThreeSections(6.256), 0.0, 0.000001);
      BCL_ExampleCheckWithinTolerance( math::WeightBetweenZeroAndPi_ThreeSections(0.002), 1.0, 0.000001);
      BCL_ExampleCheckWithinTolerance( math::WeightBetweenZeroAndPi_ThreeSections( 1.562), 0.513193, 0.000001);

      // test WeightBetweenZeroandPi
      BCL_ExampleCheckWithinTolerance( math::WeightBetweenZeroAndPi( 2.562), 0.0816571, 0.0000001);

      // test FilterValuesSmallerEqualLimit
      BCL_ExampleCheckWithinTolerance( math::FilterValuesSmallerEqualLimit( double( 2.562), double( 3.0)), 0.0, 0.001);
      BCL_ExampleCheckWithinTolerance( math::FilterValuesSmallerEqualLimit( double( 2.562), double( 2.0)), 2.562, 0.001);

      // test ConvertBooleanToSign
      BCL_ExampleCheckWithinTolerance( math::ConvertBooleanToSign( true), 1, 0.001);
      BCL_ExampleCheckWithinTolerance( math::ConvertBooleanToSign( false), -1, 0.001);

      // demonstrate that math::Absolute and std::abs always return the same values
      // the reason for using math::Absolute instead is that it ensures that we have all overloads of std::abs, otherwise
      // the compiler may choose to use the std::abs version for integers only, which would be disasterous if we had
      // a floating point type
      BCL_ExampleCheck( math::Absolute( int( -1)), std::abs( int( -1)));
      BCL_ExampleCheck( math::Absolute( int( 1)), std::abs( int( 1)));
      BCL_ExampleCheck( math::Absolute( size_t( 1)), size_t( 1));
      BCL_ExampleCheck( math::Absolute( float( -1.5)), std::abs( float( -1.5)));
      BCL_ExampleCheck( math::Absolute( float( 1.5)), std::abs( float( 1.5)));
      BCL_ExampleCheck( math::Absolute( double( -1.5)), std::abs( double( -1.5)));
      BCL_ExampleCheck( math::Absolute( double( 1.5)), std::abs( double( 1.5)));
      BCL_ExampleCheck( math::Absolute( std::complex< double>( 1.5, 1.5)), std::abs( std::complex< double>( 1.5, 1.5)));

      // test binomial coefficient
      BCL_ExampleCheck( math::BinomialCoefficient( 10, 4), math::BinomialCoefficient( 10, 6)); // test symmetry
      BCL_ExampleCheck( math::BinomialCoefficient( 10, 4), 210); // test an actual value
      BCL_ExampleCheck( math::Factorial( 4), UINT64_C( 24)); // test factorial

      // test Error Function
      BCL_ExampleCheckWithinTolerance( math::Erf( double( .1723)), 0.192513, 0.000001);

      // test complimentary error function
      BCL_ExampleCheckWithinTolerance( math::Erfc( double( .1723)), 0.807487, 0.000001);

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMath

  const ExampleClass::EnumType ExampleMath::s_Instance
  (
    GetExamples().AddEnum( ExampleMath())
  );

} // namespace bcl

