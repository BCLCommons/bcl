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

#ifndef BCL_MATH_TRIGONOMETRIC_TRANSITION_H_
#define BCL_MATH_TRIGONOMETRIC_TRANSITION_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class TrigonometricTransition
    //! @brief This class will allow the user to use and modify trigonometric functions as a function interface
    //!   it defines a cos function where (x0, y0) is one end of the function and (x1, y1) is the other end.  This
    //!   class places the cos function between these two points.  Function only works if X argument is between X1 and
    //!   X0
    //!
    //! @see @link example_math_trigonometric_transition.cpp @endlink
    //! @author akinlr, alexanns
    //! @date 06/30/2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API TrigonometricTransition :
      public FunctionInterfaceSerializable< double, double>
    {

    private:

    //////////
    // data //
    //////////

      double m_X0;  //! used to modify function as point ( x0, y0)
      double m_X1;  //! used to modify function as point ( x1, y1)
      double m_Y0;  //! used to modify function as point ( x0, y0)
      double m_Y1;  //! used to modify function as point ( x1, y1)

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      TrigonometricTransition();

      //! @brief constructor four doubles used to define a cosine function
      //! @param X0 used to modify function as point ( x0, y0)
      //! @param X1 used to modify function as point ( x1, y1)
      //! @param Y0 used to modify function as point ( x0, y0)
      //! @param Y1 used to modify function as point ( x1, y1)
      TrigonometricTransition
      (
        const double X0, const double X1, const double Y0, const double Y1
      );

      //! @brief Clone is the virtual copy constructor
      TrigonometricTransition *Clone() const
      {
        return new TrigonometricTransition( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief pretty-print the function
      std::string AsString() const;

      //! @brief get the beginning of the transition region
      //! @return the beginning of the transition region
      const double &GetTransitionBeginXAxis() const
      {
        return m_X0;
      }

      //! @brief get the beginning of the transition region
      //! @return the beginning of the transition region
      const double &GetTransitionEndXAxis() const
      {
        return m_X1;
      }

      //! @brief get the beginning of the transition region
      //! @return the beginning of the transition region
      const double &GetTransitionBeginYAxis() const
      {
        return m_Y0;
      }

      //! @brief get the beginning of the transition region
      //! @return the beginning of the transition region
      const double &GetTransitionEndYAxis() const
      {
        return m_Y1;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator() taking an x-value and returning the y-value based on "m_X0", "m_X1", "m_Y0" and "m_Y1"
      //! @param ARGUMENT double which is the x value from which the y-value will be calculated
      //! @return returns a double which is the y-value based on ARGUMENT, "m_X0", "m_X1", "m_Y0" and "m_Y1"
      double operator()( const double &X) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from istream
      //! @param ISTREAM is the input stream
      //! @return returns the input stream
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to ostream
      //! @param OSTREAM is the output stream
      //! @return returns the output stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class TrigonometricTransition

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_TRIGONOMETRIC_TRANSITION_H_
