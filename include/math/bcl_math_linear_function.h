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

#ifndef BCL_MATH_LINEAR_FUNCTION_H_
#define BCL_MATH_LINEAR_FUNCTION_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_function_interface_serializable.h"

// external includes - sorted alphabetically
#include <numeric>

namespace bcl
{
  namespace math
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LinearFunction
    //! @brief LinearFunction is a class for representing a linear function of the form y = mx + b. It can be used to
    //!        calculate the linear regression of a set of y and x values upon construction.
    //!
    //! @author alexanns
    //! @see @link example_math_linear_function.cpp @endlink
    //! @date 02/1/2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LinearFunction :
      public FunctionInterfaceSerializable< double, double>
    {

    private:

    //////////
    // data //
    //////////

      double m_Slope;             //!< the slope of the line (i.e. the "m" in y = mx + b)
      double m_OrdinateIntercept; //!< the point at which the line crosses the y-axis (i.e. the "b" in y = mx + b)

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LinearFunction();

      //! @brief construct from a slope and a y-intercept
      //! @param SLOPE double which is the slope
      //! @param ORDINATE_INTERCEPT double which is the y-intercept of the line
      LinearFunction( const double &SLOPE, const double &ORDINATE_INTERCEPT);

      //! @brief construct from a two iterator ranges and calculate the linear regression of the data to set members
      //!        The data ranges must be the same size, otherwise the linear regression will be calculated only for the
      //!        data which falls in the shorter of the two ranges. It is assumed that the corresponding x and y values
      //!        are at the same positions in the iterator ranges (i.e. first y value corresponds to first x value...)
      //! @param X_ITR_BEGIN iterator denoting the start of the x values
      //! @param X_ITR_END iterator denoting the end of the x values
      //! @param Y_ITR_BEGIN iterator denoting the start of the corresponding y values
      //! @param Y_ITR_END iterator denoting the end of the corresponding y values
      //! @param DATA_SIZE double which is the size of the data range (i.e. the number of datapoints)
      template< typename t_X_IteratorType, typename t_Y_IteratorType>
      LinearFunction
      (
        const t_X_IteratorType X_ITR_BEGIN, const t_X_IteratorType X_ITR_END,
        const t_Y_IteratorType Y_ITR_BEGIN, const t_Y_IteratorType Y_ITR_END,
        const size_t DATA_SIZE
      ) :
        m_Slope(),
        m_OrdinateIntercept()
      {
        // calculate "m_Slope" and "m_OrdinateIntercept" based on the provided data ranges
        CalculateLinearRegression( X_ITR_BEGIN, X_ITR_END, Y_ITR_BEGIN, Y_ITR_END, DATA_SIZE);
      }

      //! @brief create a linear function capable of rescaling between ranges
      //! @param FROM_RANGE the input range
      //! @param TO_RANGE the range that will be rescaled to
      template< typename t_DataType>
      LinearFunction( const Range< t_DataType> &FROM_RANGE, const Range< t_DataType> &TO_RANGE) :
        m_Slope( FROM_RANGE.GetWidth() ? double( TO_RANGE.GetWidth()) / double( FROM_RANGE.GetWidth()) : 0.0),
        m_OrdinateIntercept( TO_RANGE.GetMin() - m_Slope * FROM_RANGE.GetMin())
      {
      }

      //! @brief create a linear function capable of rescaling between ranges
      //! @param FROM_RANGE the input range
      //! @param TO_RANGE the range that will be rescaled to
      LinearFunction( const double &X1, const double &Y1, const double &X2, const double &Y2) :
        m_Slope( X1 != X2 ? ( Y2 - Y1) / ( X2 - X1) : 0.0),
        m_OrdinateIntercept( Y1 - m_Slope * X1)
      {
      }

      //! @brief virtual copy constructor
      LinearFunction *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief GetSlope return m_Slope
      //! @return returns double which is m_Slope
      const double &GetSlope() const;

      //! @brief GetY-Intercept return m_Y-Intercept
      //! @return returns double which is m_Y-Intercept
      const double &GetOrdinateIntercept() const;

      //! @brief pretty-print the scheme of the linear function
      std::string AsString() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator() taking an x-value and returning the y-value based on "m_Slope" and "m_OrdinateIntercept"
      //! @param X_VALUE double which is the x values from which the y-value will be calculated
      //! @return returns a double which is the y-value based on X_VALUE, "m_Slope" and "m_OrdinateIntercept"
      double operator()( const double &X_VALUE) const;

      //! @brief test equality
      //! @param OTHER the other linear function
      //! @return true if the linear functions have the same slope and intercept
      bool operator ==( const LinearFunction &FUNCTION) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read distance from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write distance to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief CalculateLinearRegression calculates the linear regression for a set of data denoted by two iterator
      //!        ranges. The data ranges must be the same size, otherwise the linear regression will be calculated
      //!        only for the data which falls in the x data range. It is assumed that the
      //!        corresponding x and y values are at the same positions in the iterator ranges (i.e. first y
      //!        value corresponds to first x value...)
      //! See Mortimer, R. "Mathematics for Physical Chemistry" 3rd Edition. Elsevier Academic Press: 2005. pp340 -341.
      //! @param X_ITR_BEGIN iterator denoting the start of the x values
      //! @param X_ITR_END iterator denoting the end of the x values
      //! @param Y_ITR_BEGIN iterator denoting the start of the corresponding y values
      //! @param Y_ITR_END iterator denoting the end of the corresponding y values
      //! @param DATA_SIZE double which is the size of the data range (i.e. the number of datapoints)
      template< typename t_X_IteratorType, typename t_Y_IteratorType>
      void CalculateLinearRegression
      (
        const t_X_IteratorType X_ITR_BEGIN, const t_X_IteratorType X_ITR_END,
        const t_Y_IteratorType Y_ITR_BEGIN, const t_Y_IteratorType Y_ITR_END,
        const size_t DATA_SIZE
      )
      {
        // create const double "x_sum" and initialize with the sum of the x data
        const double x_sum( std::accumulate( X_ITR_BEGIN, X_ITR_END, 0.0));

        // create const double "y_sum" and initialize with the sum of the y values
        const double y_sum( std::accumulate( Y_ITR_BEGIN, Y_ITR_END, 0.0));

        // create const double "xy_sum" and initialize with the sum of the product of the x and y values
        const double xy_sum( std::inner_product( X_ITR_BEGIN, X_ITR_END, Y_ITR_BEGIN, 0.0));

        // create const double "x_squared_sum" and initialize with the sum of the x values squared
        const double x_squared_sum( std::inner_product( X_ITR_BEGIN, X_ITR_END, X_ITR_BEGIN, 0.0));

        // create const double "d"
        const double d( DATA_SIZE * x_squared_sum - x_sum * x_sum);

        // create const double "slope" and initialize to the slope of the linear regression fit to the data
        const double slope( ( DATA_SIZE * xy_sum - x_sum * y_sum) / d);

        // create const double "ordinate_intercept" and initialize to the ordinate intercept of the linear regression
        // fit to the data
        const double ordinate_intercept( ( x_squared_sum * y_sum - x_sum * xy_sum) / d);

        // set "m_Slope" to "slope"
        m_Slope = slope;

        // set "m_OrdinateIntercept" to "ordinate_intercept"
        m_OrdinateIntercept = ordinate_intercept;
      }

    }; // class LinearFunction

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_LINEAR_FUNCTION_H_
