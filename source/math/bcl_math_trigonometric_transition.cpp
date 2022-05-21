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
#include "math/bcl_math_trigonometric_transition.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_message.h"
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
    const util::SiPtr< const util::ObjectInterface> TrigonometricTransition::s_Instance
    (
      GetObjectInstances().AddInstance( new TrigonometricTransition())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    TrigonometricTransition::TrigonometricTransition() :
      m_X0( util::GetUndefinedDouble()),
      m_X1( util::GetUndefinedDouble()),
      m_Y0( util::GetUndefinedDouble()),
      m_Y1( util::GetUndefinedDouble())
    {
    }

    //! @brief constructor four doubles used to define a cosine function
    //! @param X0 used to modify function as point ( x0, y0)
    //! @param X1 used to modify function as point ( x1, y1)
    //! @param Y0 used to modify function as point ( x0, y0)
    //! @param Y1 used to modify function as point ( x1, y1)
    TrigonometricTransition::TrigonometricTransition
    (
      const double X0, const double X1, const double Y0, const double Y1
    ) :
      m_X0( X0),
      m_X1( X1),
      m_Y0( Y0),
      m_Y1( Y1)
    {
      BCL_Assert( X1 >= X0, "Reversed x-range (was given as " + util::Format()( m_X1) + " < " + util::Format()( m_X1) + ")");
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &TrigonometricTransition::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief pretty-print the function
    std::string TrigonometricTransition::AsString() const
    {
      // print out the exact function that is used
      const double amplitude( m_Y1 - m_Y0);
      const double midpoint( 0.5 * ( m_Y1 + m_Y0));
      const double width( m_X1 - m_X0);
      const std::string x_begin( util::Format()( m_X0));
      const std::string x_end( util::Format()( m_X1));
      return
        "{ " + util::Format()( m_Y0) + " for x <= " + x_begin
        + " },{ "   + util::Format()( m_Y1) + " for x >= " + x_end
        + " },{ "   + util::Format()( midpoint) + " + " + util::Format()( amplitude)
        + " cos( pi ( x - " + x_begin + " ) / " + util::Format()( width)
        + ") for " + x_begin + " < x < " + x_end + " }";
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator() taking an x-value and returning the y-value based on "m_X0", "m_X1", "m_Y0" and "m_Y1"
    //! @param ARGUMENT double which is the x value from which the y-value will be calculated
    //! @return returns a double which is the y-value based on ARGUMENT, "m_X0", "m_X1", "m_Y0" and "m_Y1"
     double TrigonometricTransition::operator()( const double &X) const
     {
       if( X <= m_X0)
       {
         return m_Y0;
       }

       // calculate where x lies within the given cos function
       return X >= m_X1 ? m_Y1 : 0.5 * ( m_Y0 + m_Y1 - ( m_Y1 - m_Y0) * cos( g_Pi * ( ( X - m_X0) / ( m_X1 - m_X0))));
     }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &TrigonometricTransition::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_X0, ISTREAM);
      io::Serialize::Read( m_X1, ISTREAM);
      io::Serialize::Read( m_Y0, ISTREAM);
      io::Serialize::Read( m_Y1, ISTREAM);

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &TrigonometricTransition::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_X0, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_X1, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Y0, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Y1, OSTREAM, INDENT);

      return OSTREAM;
    }

  } // namespace math
} // namespace bcl
