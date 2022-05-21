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

#ifndef BCL_MATH_QUADRATIC_FUNCTION_H_
#define BCL_MATH_QUADRATIC_FUNCTION_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_math_function_interface_serializable.h"
#include "bcl_math_linear_function.h"
#include "bcl_math_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class QuadraticFunction
    //! @brief represents a Quadratic Function
    //! @details  a function of the standard form f( x) = a ( x - h) ^ 2 + k
    //! This class also has the ability to convert between standard and general form f( x) = ax ^ 2 + bx + c
    //! http://en.wikipedia.org/wiki/Quadratic_function
    //!
    //! @see @link example_math_quadratic_function.cpp @endlink
    //! @author akinlr
    //! @date Jun 24, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API QuadraticFunction :
      public FunctionInterfaceSerializable< double, double>
    {

    private:

    //////////
    // data //
    //////////

      double m_X; //!< h variable and x coordinate of the vertex
      double m_Y; //!< k variable and y coordinate of the vertex
      double m_A; //!< the a variable from the standard form

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default consructor
      QuadraticFunction();

      //! @brief construct from a, b, c variables according to the general form
      //! @param A_VARIABLE_GEN_FORM double which is the a constant of the general form
      //! @param B_VARIABLE_GEN_FORM double which is the b constant of the general form
      //! @param C_VARIABLE_GEN_FORM double which is the c constant of the general form
      QuadraticFunction
      (
        const double A_VARIABLE_GEN_FORM,
        const double B_VARIABLE_GEN_FORM,
        const double C_VARIABLE_GEN_FORM
      );

      //! @brief construct from x and y coordinates of the vertex and an A variable according to the standard form
      //! @param XY_COORD VectorND which contains the x and y coordinates of the vertex from the standard form
      //! @param A_VARIABLE double which is the a variable from the standard form
      QuadraticFunction
      (
        const storage::VectorND< 2, double> &XY_COORD,
        const double A_VARIABLE
      );

      //! @brief Clone function
      //! @return pointer to new QuadraticFunction
      QuadraticFunction *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief GetX return m_X
      //! @return returns double which is m_X
      double GetX() const;

      //! @brief GetY return m_Y
      //! @return returns double which is m_Y
      double GetY() const;

      //! @brief GetA return m_A
      //! @return returns double which is m_A
      double GetA() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief calculates the derivative of the quadratic function
      //! @param QUADRATIC the function which will have its derivative taken
      //! @return the derivative of the quadratic function
      util::ShPtr< FunctionInterfaceSerializable< double, double> > GetDerivative() const;

      //! @brief calculates the roots of the given quadratic function
      //! @param QUADRATIC the function whose roots will be determined
      //! @return the pair of roots of the quadratic function
      storage::VectorND< 2, double> GetRoot() const;

      //! @brief pretty-print the function
      std::string AsString() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator() taking an x-value and returning the y-value based on "m_X", "m_Y" and "m_A"
      //! @param ARGUMENT double which is the x value from which the y-value will be calculated
      //! @return returns a double which is the y-value based on ARGUMENT, "m_X", "m_Y" and "m_A"
      double operator()( const double &ARGUMENT) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief converts general form variables to standard form variables stored as x, y, and a in that order
      //! @param A_VALUE double that is the a variable from the general form
      //! @param B_VALUE double that is the b variable from the general form
      //! @param C_VALUE double that is the c variable from the general form
      //! @return VectorND which contains the standard form variables x, y, and a in that order
      static storage::VectorND< 3, double> CompletingTheSquare
      (
        const double A_VALUE, const double B_VALUE, const double C_VALUE
      );

      //! @brief converts standard form variables to general form variables stored as a, b, and c in that order
      //! @param X_COORD double that is the x variable from the standard form
      //! @param Y_COORD double that is the y variable from the standard form
      //! @param A_VARIABLE double that is the a variable from the standard form
      //! @return VectorND which contains the general form variables a, b, and c in that order
      static storage::VectorND< 3, double> Distribution
      (
        const double X_COORD, const double Y_COORD, const double A_VARIABLE
      );

    }; // class QuadraticFunction

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_QUADRATIC_FUNCTION_H_
