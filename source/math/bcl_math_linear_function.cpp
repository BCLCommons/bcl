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
#include "math/bcl_math_linear_function.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> LinearFunction::s_Instance
    (
      GetObjectInstances().AddInstance( new LinearFunction())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LinearFunction::LinearFunction() :
      m_Slope( util::GetUndefinedDouble()),
      m_OrdinateIntercept( util::GetUndefinedDouble())
    {
    }

    //! @brief construct from a slope and a y-intercept
    //! @param SLOPE double which is the slope
    //! @param ORDINATE_INTERCEPT double which is the y-intercept of the line
    LinearFunction::LinearFunction( const double &SLOPE, const double &ORDINATE_INTERCEPT) :
      m_Slope( SLOPE),
      m_OrdinateIntercept( ORDINATE_INTERCEPT)
    {
    }

    //! @brief virtual copy constructor
    LinearFunction *LinearFunction::Clone() const
    {
      return new LinearFunction( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LinearFunction::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief GetSlope return m_Slope
    //! @return returns double which is m_Slope
    const double &LinearFunction::GetSlope() const
    {
      return m_Slope;
    }

    //! @brief GetOrdinateIntercept return m_OrdinateIntercept
    //! @return returns double which is m_OrdinateIntercept
    const double &LinearFunction::GetOrdinateIntercept() const
    {
      return m_OrdinateIntercept;
    }

    //! @brief pretty-print the scheme of the linear function
    std::string LinearFunction::AsString() const
    {
      return util::Format()( m_Slope) + " x + " + util::Format()( m_OrdinateIntercept);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator() taking an x-value and returning the y-value based on "m_Slope" and "m_OrdinateIntercept"
    //! @param X_VALUE double which is the x values from which the y-value will be calculated
    //! @return returns a double which is the y-value based on X_VALUE, "m_Slope" and "m_OrdinateIntercept"
    double LinearFunction::operator()( const double &X_VALUE) const
    {
      return m_Slope * X_VALUE + m_OrdinateIntercept;
    }

    //! @brief test equality
    //! @param OTHER the other linear function
    //! @return true if the linear functions have the same slope and intercept
    bool LinearFunction::operator ==( const LinearFunction &FUNCTION) const
    {
      return m_Slope == FUNCTION.m_Slope && m_OrdinateIntercept == FUNCTION.m_OrdinateIntercept;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read distance from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LinearFunction::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Slope, ISTREAM);
      io::Serialize::Read( m_OrdinateIntercept, ISTREAM);

      return ISTREAM;
    }

    //! @brief write distance to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &LinearFunction::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Slope, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_OrdinateIntercept, OSTREAM, 0);

      return OSTREAM;
    }

  } // namespace math
} // namespace bcl
