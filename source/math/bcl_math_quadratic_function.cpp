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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "math/bcl_math_quadratic_function.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> QuadraticFunction::s_Instance
    (
      GetObjectInstances().AddInstance( new QuadraticFunction())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default consructor
    QuadraticFunction::QuadraticFunction() :
      m_X( util::GetUndefinedDouble()),
      m_Y( util::GetUndefinedDouble()),
      m_A( util::GetUndefinedDouble())
    {
    }

    //! @brief construct from a, b, c variables according to the general form
    //! @param A_VARIABLE_GEN_FORM double which is the a constant of the general form
    //! @param B_VARIABLE_GEN_FORM double which is the b constant of the general form
    //! @param C_VARIABLE_GEN_FORM double which is the c constant of the general form
    QuadraticFunction::QuadraticFunction
    (
      double A_VARIABLE_GEN_FORM, double B_VARIABLE_GEN_FORM, double C_VARIABLE_GEN_FORM
    ) :
      m_X(),
      m_Y(),
      m_A()
    {
      // create VectorND "xya_values" and initialize the results of the CompletingTheSquare Function
      storage::VectorND< 3, double> xya_values
      (
        CompletingTheSquare( A_VARIABLE_GEN_FORM, B_VARIABLE_GEN_FORM, C_VARIABLE_GEN_FORM)
      );

      // assign member variables to the appropriate values
      m_X = xya_values.First();
      m_Y = xya_values.Second();
      m_A = xya_values.Third();
    }

    //! @brief construct from x and y coordinates of the vertex and an A variable according to the standard form
    //! @param XY_COORD VectorND which contains the x and y coordinates of the vertex from the standard form
    //! @param A_VARIABLE double which is the a variable from the standard form
    QuadraticFunction::QuadraticFunction
    (
      const storage::VectorND< 2, double> &XY_COORD,
      const double A_VARIABLE
    ) :
      m_X( XY_COORD.First()),
      m_Y( XY_COORD.Second()),
      m_A( A_VARIABLE)
    {
    }

    //! @brief virtual copy constructor
    QuadraticFunction *QuadraticFunction::Clone() const
    {
      return new QuadraticFunction( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &QuadraticFunction::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief GetX return m_X
    //! @return returns double which is m_X
    double QuadraticFunction::GetX() const
    {
      return m_X;
    }

    //! @brief GetY return m_Y
    //! @return returns double which is m_Y
    double QuadraticFunction::GetY() const
    {
      return m_Y;
    }

    //! @brief GetA return m_A
    //! @return returns double which is m_A
    double QuadraticFunction::GetA() const
    {
      return m_A;
    }

    //! @brief pretty-print the scheme of the linear function
    std::string QuadraticFunction::AsString() const
    {
      return util::Format()( m_A) + " ( x - " + util::Format()( m_X) + ")^2 + " + util::Format()( m_Y);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculates the derivative of the quadratic function
    //! @param QUADRATIC the function which will have its derivative taken
    //! @return the derivative of the quadratic function
    util::ShPtr< FunctionInterfaceSerializable< double, double> > QuadraticFunction::GetDerivative() const
    {
      // put the variables into the general form
      const storage::VectorND< 3, double> general_form( Distribution( m_X, m_Y, m_A));

      // use the constants from the general form to find the derivative
      const double a_linear_variable( general_form.First());
      const double b_linear_variable( general_form.Second());
      const double slope( a_linear_variable * 2);
      const double y_intercept( b_linear_variable);

      // create the derivative
      return util::ShPtr< FunctionInterfaceSerializable< double, double> >( new LinearFunction( slope, y_intercept));
    }

    //! @brief calculates the roots of the given quadratic function
    //! @param QUADRATIC the function whose roots will be determined
    //! @return the pair of roots of the quadratic function
    storage::VectorND< 2, double> QuadraticFunction::GetRoot() const
    {
      // put the variables into the general form
      const storage::VectorND< 3, double> general_form( Distribution( m_X, m_Y, m_A));
      const double a( general_form.First());
      const double b( general_form.Second());
      const double c( general_form.Third());

      // calculate the roots
      const double denominator( 2 * a);
      const double sqrt_result( Sqrt( b * b - 4.0 * a * c));
      const double numerator_pos( -b + sqrt_result);
      const double numerator_neg( -b - sqrt_result);
      const double root_a( numerator_pos / denominator);
      const double root_b( numerator_neg / denominator);
      return storage::VectorND< 2, double>( root_a, root_b);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator() taking an x-value and returning the y-value based on "m_X", "m_Y" and "m_A"
    //! @param ARGUMENT double which is the x value from which the y-value will be calculated
    //! @return returns a double which is the y-value based on ARGUMENT, "m_X", "m_Y", and "m_A"
     double QuadraticFunction::operator()( const double &ARGUMENT) const
     {
       return ( m_A * ( Sqr( ARGUMENT - m_X))) + m_Y;
     }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &QuadraticFunction::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_X, ISTREAM);
      io::Serialize::Read( m_Y, ISTREAM);
      io::Serialize::Read( m_A, ISTREAM);

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &QuadraticFunction::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_X, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Y, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_A, OSTREAM, INDENT);

      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief converts general form variables to standard form variables stored as x, y, and a in that order
    //! @param A_VALUE double that is the a variable from the general form
    //! @param B_VALUE double that is the b variable from the general form
    //! @param C_VALUE double that is the c variable from the general form
    //! @return VectorND which contains the standard form variables x, y, and a in that order
    storage::VectorND< 3, double> QuadraticFunction::CompletingTheSquare
    (
      const double A_VALUE, const double B_VALUE, const double C_VALUE
    )
    {
      // convert the general form parameters into the standard form parameters
      // check to make sure the A_VALUE is not equal to zero
      const double a_variable( A_VALUE);

      // if A is equal to zero return undefined doubles
      if( a_variable == 0)
      {
        return storage::VectorND< 3, double>
        (
          util::GetUndefinedDouble(), util::GetUndefinedDouble(), util::GetUndefinedDouble()
        );
      }

      // if A is not equal to zero complete the square
      else
      {
        const double x_coord( -B_VALUE / ( 2.0 * a_variable));
        const double y_coord( C_VALUE - Sqr( x_coord) * a_variable);

        // return the VectorND with the standard form "x_coord", "y_coord" and "a_variable" in that order
        return storage::VectorND< 3, double>( x_coord, y_coord, a_variable);
      }
    }

    //! @brief converts standard form variables to general form variables stored as a, b, and c in that order
    //! @param X_COORD double that is the x variable from the standard form
    //! @param Y_COORD double that is the y variable from the standard form
    //! @param A_VARIABLE double that is the a variable from the standard form
    //! @return VectorND which contains the general form variables a, b, and c in that order
    storage::VectorND< 3, double> QuadraticFunction::Distribution
    (
      const double X_COORD, const double Y_COORD, const double A_VARIABLE
    )
    {
      // convert the standard form parameters into the general form parameters
      const double a_value( A_VARIABLE);
      const double b_value( -2 * X_COORD * A_VARIABLE);
      const double c_value( Sqr( X_COORD) * A_VARIABLE + Y_COORD);

      // return the VectorND with the general form "a_value", "b_value" and "c_value" in that order
      return storage::VectorND< 3, double>( a_value, b_value, c_value);
    }

  } // namespace math
} // namespace bcl
