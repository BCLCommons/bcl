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
#include "math/bcl_math_angle.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

  //////////
  // data //
  //////////

    //! string used to represent the degree character in gnuplot scripts
    const std::string Angle::s_DegreeSymbolGnuplot = "\\260";

    //! @brief conversion to a string from a Unit
    //! @param UNIT the Unit to get a string for
    //! @return a string representing that Unit
    const std::string &Angle::GetUnitName( const Unit &UNIT)
    {
      static const std::string s_descriptors[] =
      {
        "radian",
        "degree",
        "undefined",
        "Angle"
      };
      return s_descriptors[ size_t( UNIT)];
    }

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_assign.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Assign< double>;
    template class BCL_API Assign< float>;

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_assignment_by_comparison.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    // register the class with the enumerated class instance
    template< typename t_ArgumentType, template< typename> class t_ComparisonType>
    const util::SiPtr< const AssignmentOperationInterface< t_ArgumentType> >
    AssignmentByComparison< t_ArgumentType, t_ComparisonType>::s_Instance
    (
      util::Enumerated< AssignmentOperationInterface< t_ArgumentType> >::AddInstance
      (
        new AssignmentByComparison< t_ArgumentType, t_ComparisonType>()
      )
    );

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API AssignmentByComparison< float, std::less>;
    template class BCL_API AssignmentByComparison< float, std::less_equal>;
    template class BCL_API AssignmentByComparison< float, std::greater>;
    template class BCL_API AssignmentByComparison< float, std::greater_equal>;
    template class BCL_API AssignmentByComparison< float, std::equal_to>;
    template class BCL_API AssignmentByComparison< float, std::not_equal_to>;

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_assignments.hpp"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Assignments< double>;
    template class BCL_API Assignments< float>;

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_assignment_unary_interface.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_assignment_unary_standard.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    //! @brief get all the different unary assignments
    //! @return vector of all different unary assignment types
    const storage::Vector< util::Implementation< AssignmentUnaryInterface> >
      &AssignmentUnaryInterface::GetUnaryAssignments()
    {
      static const storage::Vector< util::Implementation< AssignmentUnaryInterface> > s_impls
      (
        storage::Vector< util::Implementation< AssignmentUnaryInterface> >::Create
        (
          util::Implementation< AssignmentUnaryInterface>( new AssignmentUnaryStandard( cos)),
          util::Implementation< AssignmentUnaryInterface>( new AssignmentUnaryStandard( sin)),
          util::Implementation< AssignmentUnaryInterface>( new AssignmentUnaryStandard( log)),
          util::Implementation< AssignmentUnaryInterface>( new AssignmentUnaryStandard( log10)),
          util::Implementation< AssignmentUnaryInterface>( new AssignmentUnaryStandard( exp)),
          util::Implementation< AssignmentUnaryInterface>( new AssignmentUnaryStandard( std::abs)),
          util::Implementation< AssignmentUnaryInterface>( new AssignmentUnaryStandard( std::sqrt)),
          util::Implementation< AssignmentUnaryInterface>( new AssignmentUnaryStandard( AssignmentUnaryStandard::e_Square)),
          util::Implementation< AssignmentUnaryInterface>( new AssignmentUnaryStandard( AssignmentUnaryStandard::e_Not)),
          util::Implementation< AssignmentUnaryInterface>( new AssignmentUnaryStandard( AssignmentUnaryStandard::e_Negative))
        )
      );
      return s_impls;
    }

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_assignment_unary_standard.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    namespace
    {
      //! @brief square without const ref
      //! @param VALUE the value to square
      double Square( double VALUE)
      {
        return VALUE * VALUE;
      }

      //! @brief square without const ref
      //! @param VALUE the value to not
      double Not( double VALUE)
      {
        return !VALUE;
      }

      //! @brief square without const ref
      //! @param VALUE the value to negate
      double Negative( double VALUE)
      {
        return -VALUE;
      }
    }

    //! @brief function to add all the instances
    util::SiPtr< const util::ObjectInterface> AssignmentUnaryStandard::AddInstances()
    {
      util::Enumerated< AssignmentUnaryInterface>::AddInstance( new AssignmentUnaryStandard( cos));
      util::Enumerated< AssignmentUnaryInterface>::AddInstance( new AssignmentUnaryStandard( sin));
      util::Enumerated< AssignmentUnaryInterface>::AddInstance( new AssignmentUnaryStandard( log));
      util::Enumerated< AssignmentUnaryInterface>::AddInstance( new AssignmentUnaryStandard( log10));
      util::Enumerated< AssignmentUnaryInterface>::AddInstance( new AssignmentUnaryStandard( exp));
      util::Enumerated< AssignmentUnaryInterface>::AddInstance( new AssignmentUnaryStandard( std::abs));
      util::Enumerated< AssignmentUnaryInterface>::AddInstance( new AssignmentUnaryStandard( std::sqrt));
      util::Enumerated< AssignmentUnaryInterface>::AddInstance( new AssignmentUnaryStandard( Square));
      util::Enumerated< AssignmentUnaryInterface>::AddInstance( new AssignmentUnaryStandard( Not));
      util::Enumerated< AssignmentUnaryInterface>::AddInstance( new AssignmentUnaryStandard( Negative));
      return util::SiPtr< const util::ObjectInterface>();
    }

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    // register the class with the enumerated class instance
    const util::SiPtr< const util::ObjectInterface> AssignmentUnaryStandard::s_Instance( AddInstances());

    //! @brief helper function, get alias and description for a function
    //! @param FUNCTION the function of interest
    std::pair< std::string, std::string> AssignmentUnaryStandard::GetUnaryAssignmentAliasDescription
    (
      double ( *FUNCTION)( double)
    )
    {
      typedef double ( *t_DoubleFunction)( double);
      if( FUNCTION == static_cast< double( *)( double)>( &cos))
      {
        return std::pair< std::string, std::string>( "Cos", "Takes the cosine");
      }
      else if( FUNCTION == static_cast< double( *)( double)>( &sin))
      {
        return std::pair< std::string, std::string>( "Sin", "Takes the sine");
      }
      else if( FUNCTION == static_cast< double( *)( double)>( &log))
      {
        return std::pair< std::string, std::string>( "Ln", "Takes the natural log");
      }
      else if( FUNCTION == static_cast< double( *)( double)>( &log10))
      {
        return std::pair< std::string, std::string>( "Log", "Takes the base-10 logarithm");
      }
      else if( FUNCTION == static_cast< double( *)( double)>( &exp))
      {
        return std::pair< std::string, std::string>( "Exp", "Takes the exponential");
      }
      else if( FUNCTION == t_DoubleFunction( &std::abs))
      {
        return std::pair< std::string, std::string>( "Abs", "Takes the absolute value");
      }
      else if( FUNCTION == t_DoubleFunction( &std::sqrt))
      {
        return std::pair< std::string, std::string>( "Sqrt", "Takes the square root");
      }
      else if( FUNCTION == &Square)
      {
        return std::pair< std::string, std::string>( "Sqr", "squares the argument");
      }
      else if( FUNCTION == &Not)
      {
        return std::pair< std::string, std::string>( "Not", "1 if the argument is exactly 0, otherwise returns 0");
      }
      else if( FUNCTION == &Negative)
      {
        return std::pair< std::string, std::string>( "Negative", "gives the negative of the given argument");
      }
      return std::pair< std::string, std::string>( "", "undefined");
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from alias
    //! @param FUNCTION the function to use
    AssignmentUnaryStandard::AssignmentUnaryStandard( double ( *FUNCTION)( double)) :
      m_Function( FUNCTION),
      m_Alias( GetUnaryAssignmentAliasDescription( m_Function).first)
    {
    }

    //! @brief constructor from local enum
    AssignmentUnaryStandard::AssignmentUnaryStandard( FunctionType FUNC) :
      m_Function( FUNC == e_Not ? &Not : FUNC == e_Negative ? &Negative : &Square),
      m_Alias( GetUnaryAssignmentAliasDescription( m_Function).first)
    {
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AssignmentUnaryStandard::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( GetUnaryAssignmentAliasDescription( m_Function).second);
      return parameters;
    }

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_bicubic_spline.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_cubic_spline.h"
#include "math/bcl_math_smooth_data.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> BicubicSpline::s_Instance
    (
      GetObjectInstances().AddInstance( new BicubicSpline())
    );

  //////////////
  // operator //
  //////////////

    //! @brief operator
    //! @param X argument x to spline
    //! @param Y argument y to spline
    //! @return the interpolate value at (x,y)
    double BicubicSpline::operator()( const double &X, const double &Y) const
    {
      return F( linal::MakeVector< double>( X, Y));
    }

  ////////////////
  // operations //
  ////////////////

    //! train BicubicSpline
    void BicubicSpline::Train( const SplineBorderType BORDER[ 2], const double START[ 2], const double DELTA[ 2], const linal::Matrix< double> &RESULTS, const bool LINCONT[ 2], const storage::Pair< double, double> FIRSTBE[ 2])
    {
      //check, if the points are given in positive direction
      BCL_Assert( DELTA[ coord::GetAxes().e_X] > 0, "deltax <= 0 not supported");
      BCL_Assert( DELTA[ coord::GetAxes().e_Y] > 0, "deltay <= 0 not supported");

      //determine values for all dimensions
      const size_t dimx( RESULTS.GetNumberRows());
      const size_t dimy( RESULTS.GetNumberCols());

      //assigning values
      m_Border[ coord::GetAxes().e_X] = BORDER[ coord::GetAxes().e_X];
      m_Border[ coord::GetAxes().e_Y] = BORDER[ coord::GetAxes().e_Y];
      m_Start[ coord::GetAxes().e_X]  = START[ coord::GetAxes().e_X];
      m_Start[ coord::GetAxes().e_Y]  = START[ coord::GetAxes().e_Y];
      m_Delta[ coord::GetAxes().e_X]  = DELTA[ coord::GetAxes().e_X];
      m_Delta[ coord::GetAxes().e_Y]  = DELTA[ coord::GetAxes().e_Y];
      m_DeltaNorm[ coord::GetAxes().e_X]  = Sqr( DELTA[ coord::GetAxes().e_X]) / 6.0; // cache for speed
      m_DeltaNorm[ coord::GetAxes().e_Y]  = Sqr( DELTA[ coord::GetAxes().e_Y]) / 6.0; // cache for speed

      BCL_Assert( m_Border[ 0] != e_NotAKnot && m_Border[ 1] != e_NotAKnot, "The specified boundary condition is not yet supported");

      m_Values  = RESULTS;
      m_Dsecox  = RESULTS;
      m_Dsecoy  = RESULTS;
      m_Dsecoxy = RESULTS;

      m_LinCont[ coord::GetAxes().e_X] = LINCONT[ coord::GetAxes().e_X];
      m_LinCont[ coord::GetAxes().e_Y] = LINCONT[ coord::GetAxes().e_Y];
      m_FirstBe[ coord::GetAxes().e_X] = FIRSTBE[ coord::GetAxes().e_X];
      m_FirstBe[ coord::GetAxes().e_Y] = FIRSTBE[ coord::GetAxes().e_Y];

      //train three times for fxx, fyy, fxxyy
      //reduction to Spline by training only rows/columns at the same time
      for( size_t row( 0); row < dimx; ++row)
      {
        CubicSpline cs;
        cs.Train( BORDER[ coord::GetAxes().e_Y], START[ coord::GetAxes().e_Y], DELTA[ coord::GetAxes().e_Y], RESULTS.GetRow( row), FIRSTBE[ coord::GetAxes().e_Y]);
        m_Dsecoy.ReplaceRow( row, cs.GetDsecox());
      }

      for( size_t col( 0); col < dimy; ++col)
      {
        CubicSpline cs;
        cs.Train( BORDER[ coord::GetAxes().e_X], START[ coord::GetAxes().e_X], DELTA[ coord::GetAxes().e_X], RESULTS.GetCol( col), FIRSTBE[ coord::GetAxes().e_X]);
        m_Dsecox.ReplaceCol( col, cs.GetDsecox());
      }

      for( size_t row( 0); row < dimx; ++row)
      {
        CubicSpline cs;
        cs.Train( BORDER[ coord::GetAxes().e_Y], START[ coord::GetAxes().e_Y], DELTA[ coord::GetAxes().e_Y], m_Dsecox.GetRow( row), FIRSTBE[ coord::GetAxes().e_Y]);
        m_Dsecoxy.ReplaceRow( row, cs.GetDsecox());
      }
      return;
    }

    void BicubicSpline::TrainWithPreprocessing( const SplineBorderType BORDER[ 2], const double START[ 2], const double DELTA[ 2], const linal::Matrix< double> &RESULTS, const bool LINCONT[ 2], const storage::Pair< double, double> FIRSTBE[ 2], const double PREPROC)
    {
      // passing of smoothed data set
      Train( BORDER, START, DELTA, SmoothData::SmoothMatrix( RESULTS, PREPROC, true), LINCONT, FIRSTBE);
    }

    //! return value at certain (x, y)
    double BicubicSpline::F( const linal::Vector< double> &ARGUMENTS) const
    {
      // check that there are two argument values given
      BCL_Assert( ARGUMENTS.GetSize() == 2, "number of arguments doesn't match");

      const double x( ARGUMENTS( coord::GetAxes().e_X));
      const double y( ARGUMENTS( coord::GetAxes().e_Y));

      const int dimx( m_Values.GetNumberRows());
      const int dimy( m_Values.GetNumberCols());

      // check if argument x is in range for non-periodic splines
      if( x < m_Start[ coord::GetAxes().e_X])
      {
        switch( m_Border[ coord::GetAxes().e_X])
        {
          case e_FirstDer: BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");
                           BCL_MessageDbg
                           (
                             "argument out of range (X < start : "
                             + util::Format()( x)
                             + " < "
                             + util::Format()( m_Start[ coord::GetAxes().e_X])
                             + "), using linear continuation"
                           );
                           return F( linal::MakeVector( m_Start[ coord::GetAxes().e_X], y)) + ( x - m_Start[ coord::GetAxes().e_X]) * m_FirstBe[ coord::GetAxes().e_X].First();

          case e_Natural: BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");
                          BCL_MessageDbg
                          (
                            "argument out of range (X < start : "
                            + util::Format()( x)
                            + " < "
                            + util::Format()( m_Start[ coord::GetAxes().e_X])
                            + "), using linear continuation"
                          );
                          return F( linal::MakeVector( m_Start[ coord::GetAxes().e_X], y)) + ( x - m_Start[ coord::GetAxes().e_X]) * dFdx( linal::MakeVector( m_Start[ coord::GetAxes().e_X], y));
          case e_Periodic: break;
          case e_NotAKnot: break;
        }
      }

      if( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X] < x)
      {
        switch( m_Border[ coord::GetAxes().e_X])
        {
          case e_FirstDer: BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");
                           BCL_MessageDbg
                           (
                             "argument out of range (X > end : "
                             + util::Format()( x)
                             + " > "
                             + util::Format()( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X])
                             + "), using linear continuation"
                           );
                           return F( linal::MakeVector( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X], y)) + ( x - m_Start[ coord::GetAxes().e_X] - ( dimx - 1) * m_Delta[ coord::GetAxes().e_X]) * m_FirstBe[ coord::GetAxes().e_X].Second();

          case e_Natural: BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");
                          BCL_MessageDbg
                          (
                            "argument out of range (X > end : "
                            + util::Format()( x)
                            + " > "
                            + util::Format()( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X])
                            + "), using linear continuation"
                          );
                          return F( linal::MakeVector( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X], y)) + ( x - m_Start[ coord::GetAxes().e_X] - ( dimx - 1) * m_Delta[ coord::GetAxes().e_X]) * dFdx( linal::MakeVector( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X], y));

          case e_Periodic: break;
          case e_NotAKnot: break;
        }
      }

      //check if argument y is in range for non-periodic splines
      if( y < m_Start[ coord::GetAxes().e_Y])
      {
        switch( m_Border[ coord::GetAxes().e_Y])
        {
          case e_FirstDer: BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");
                           BCL_MessageDbg
                           (
                             "argument out of range (y < start : "
                             + util::Format()( y)
                             + " < "
                             + util::Format()( m_Start[ coord::GetAxes().e_Y])
                             + "), using linear continuation"
                           );
                           return F( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y])) + ( y - m_Start[ coord::GetAxes().e_Y]) * m_FirstBe[ coord::GetAxes().e_Y].First();

          case e_Natural: BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");
                          BCL_MessageDbg
                          (
                            "argument out of range (y < start : "
                            + util::Format()( y)
                            + " < "
                            + util::Format()( m_Start[ coord::GetAxes().e_Y])
                            + "), using linear continuation"
                          );
                          return F( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y])) + ( y - m_Start[ coord::GetAxes().e_Y]) * dFdy( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y]));

          case e_Periodic: break;
          case e_NotAKnot: break;
        }
      }

      if( m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y] < y)
      {
        switch( m_Border[ coord::GetAxes().e_Y])
        {
          case e_FirstDer: BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");
                           BCL_MessageDbg
                           (
                             "argument out of range (y > end : "
                             + util::Format()( y)
                             + " > "
                             + util::Format()( m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y])
                             + "), using linear continuation"
                           );
                          return F( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y])) + ( y - m_Start[ coord::GetAxes().e_Y] - ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y]) * m_FirstBe[ coord::GetAxes().e_Y].Second();

          case e_Natural: BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");
                          BCL_MessageDbg
                          (
                            "argument out of range (y > end : "
                            + util::Format()( y)
                            + " > "
                            + util::Format()( m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y])
                            + "), using linear continuation"
                          );
                          return F( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y])) + ( y - m_Start[ coord::GetAxes().e_Y] - ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y]) * dFdy( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y]));

          case e_Periodic: break;
          case e_NotAKnot: break;
        }
      }

      // determine i with m_Start[ coord::GetAxes().e_X]+(i-1)*m_Delta[ coord::GetAxes().e_X] < x < m_Start[ coord::GetAxes().e_X]+i*m_Delta[ coord::GetAxes().e_X] for the correct supporting points
      int i( int( floor( ( x - m_Start[ coord::GetAxes().e_X]) / m_Delta[ coord::GetAxes().e_X]) + 1));

      // determine j with m_Start[ coord::GetAxes().e_Y]+(j-1)*m_Delta[ coord::GetAxes().e_Y] < y < m_Start[ coord::GetAxes().e_Y]+j*m_Delta[ coord::GetAxes().e_Y] for the correct supporting points
      int j( int( floor( ( y - m_Start[ coord::GetAxes().e_Y]) / m_Delta[ coord::GetAxes().e_Y]) + 1));

      // see mkbicubw.f in /home/mj/pspline/pspline
      // slightly corrected formula, since the last two lines of the expression weren't symmetric

//C  on input:  f(1,i,j) = f(x(i),y(j))
//C  on output:  f(1,i,j) unchanged
//C              f(2,i,j) = d2f/dx2(x(i),y(j))
//C              f(3,i,j) = d2f/dy2(x(i),y(j))
//C              f(4,i,j) = d4f/dx2dy2(x(i),y(j))
//C
//C  and the interpolation formula for (x,y) in (x(i),x(i+1))^(y(j),y(j+1))
//C  is:
//C        hx = x(i+1)-x(i)   hy = y(j+1)-y(j)
//C        dxp= (x-x(i))/hx   dxm= 1-dxp     dxp,dxm in (0,1)
//C        dyp= (x-x(i))/hx   dym= 1-dyp     dyp,dym in (0,1)
//C        dx3p = dxp**3-dxp  dx3m = dxm**3-dxm     dxp3,dxm3 in (0,1)
//C
//C   finterp = dxm*(dym*f(1,i,j)+dyp*f(1,i,j+1))
//C            +dxp*(dym*f(1,i+1,j)+dyp*f(1,i+1,j+1))
//C     +1/6*hx**2*
//C            dx3m*(dym*f(2,i,j)+dyp*f(2,i,j+1))
//C           +dx3p*(dym*f(2,i+1,j)+dyp*f(2,i+1,j+1))
//C     +1/6*hy**2*
//C            dxm*(dy3m*f(3,i,j)+dy3p*f(3,i,j+1))
//C           +dxp*(dy3m*f(3,i+1,j)+dy3p*f(3,i+1,j+1))
//C     +1/36*hx**2*hy**2*
//C            dx3m*(dy3m*f(4,i,j)+dy3p*f(4,i,j+1))
//C           +dx3p*(dy3m*f(4,i+1,j)+dy3p*f(4,i+1,j+1))

      const double dxp( ( x - m_Start[ coord::GetAxes().e_X]) / m_Delta[ coord::GetAxes().e_X] - floor( ( x - m_Start[ coord::GetAxes().e_X]) / m_Delta[ coord::GetAxes().e_X]));
      const double dxm( 1 - dxp);
      const double dx3p( ( dxp*dxp*dxp - dxp) * m_DeltaNorm[ coord::GetAxes().e_X]); // =0 at the grid points, adds cubic part of the spline
      const double dx3m( ( dxm*dxm*dxm - dxm) * m_DeltaNorm[ coord::GetAxes().e_X]); // =0 at the grid points, adds cubic part of the spline

      const double dyp( ( y - m_Start[ coord::GetAxes().e_Y]) / m_Delta[ coord::GetAxes().e_Y] - floor( ( y - m_Start[ coord::GetAxes().e_Y]) / m_Delta[ coord::GetAxes().e_Y]));
      const double dym( ( 1 - dyp));
      const double dy3p( ( dyp*dyp*dyp - dyp) * m_DeltaNorm[ coord::GetAxes().e_Y]); // =0 at the grid points, adds cubic part of the spline
      const double dy3m( ( dym*dym*dym - dym) * m_DeltaNorm[ coord::GetAxes().e_Y]); // =0 at the grid points, adds cubic part of the spline

      // generate positive values to prevent some problems with the indices
      while( i < 1) i += dimx;
      while( j < 1) j += dimy;

      return
           dxm * (  dym * m_Values( (  i - 1) % dimx, ( j - 1) % dimy) +  dyp * m_Values( (  i - 1) % dimx, j % dimy))
        +  dxp * (  dym * m_Values(    i      % dimx, ( j - 1) % dimy) +  dyp * m_Values(    i      % dimx, j % dimy))
        + dx3m * (  dym * m_Dsecox( (  i - 1) % dimx, ( j - 1) % dimy) +  dyp * m_Dsecox( (  i - 1) % dimx, j % dimy))
        + dx3p * (  dym * m_Dsecox(    i      % dimx, ( j - 1) % dimy) +  dyp * m_Dsecox(    i      % dimx, j % dimy))
        +  dxm * ( dy3m * m_Dsecoy( (  i - 1) % dimx, ( j - 1) % dimy) + dy3p * m_Dsecoy( (  i - 1) % dimx, j % dimy))
        +  dxp * ( dy3m * m_Dsecoy(    i      % dimx, ( j - 1) % dimy) + dy3p * m_Dsecoy(    i      % dimx, j % dimy))
        + dx3m * ( dy3m * m_Dsecoxy( ( i - 1) % dimx, ( j - 1) % dimy) + dy3p * m_Dsecoxy( ( i - 1) % dimx, j % dimy))
        + dx3p * ( dy3m * m_Dsecoxy(   i      % dimx, ( j - 1) % dimy) + dy3p * m_Dsecoxy(   i      % dimx, j % dimy))
      ;
    }

    //! return partial derivative at certain (x, y) for x
    double BicubicSpline::dFdx( const linal::Vector< double> &ARGUMENTS) const
    {
      BCL_Assert( ARGUMENTS.GetSize() == 2, "number of arguments doesn't match");

      const double x( ARGUMENTS(0));
      const double y( ARGUMENTS(1));

      const int dimx( m_Values.GetNumberRows());
      const int dimy( m_Values.GetNumberCols());

      //check if argument x is in range for non-periodic splines
      if( x < m_Start[ coord::GetAxes().e_X])
        switch( m_Border[ coord::GetAxes().e_X])
        {
          case e_FirstDer: BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return m_FirstBe[ coord::GetAxes().e_X].First();
          case e_Natural : BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return dFdx( linal::MakeVector( m_Start[ coord::GetAxes().e_X], y));
          case e_Periodic: break;
          case e_NotAKnot: break;
        }

      if( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X] < x)
        switch( m_Border[ coord::GetAxes().e_X])
        {
          case e_FirstDer: BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return m_FirstBe[ coord::GetAxes().e_X].Second();
          case e_Natural : BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return dFdx( linal::MakeVector( m_Start[ coord::GetAxes().e_X] + ( dimx-1 ) * m_Delta[ coord::GetAxes().e_X], y));
          case e_Periodic: break;
          case e_NotAKnot: break;
        }

      //check if argument y is in range for non-periodic splines
      if( y < m_Start[ coord::GetAxes().e_Y])
        switch( m_Border[ coord::GetAxes().e_Y])
        {
          case e_FirstDer: BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return dFdx( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y]));
          case e_Natural : BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return dFdx( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y]));
          case e_Periodic: break;
          case e_NotAKnot: break;
        }

      if( m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y] < y)
        switch( m_Border[ coord::GetAxes().e_Y])
        {
          case e_FirstDer: BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return dFdx( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y] + ( dimy-1 ) * m_Delta[ coord::GetAxes().e_Y]));
          case e_Natural : BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return dFdx( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y] + ( dimy-1 ) * m_Delta[ coord::GetAxes().e_Y]));
          case e_Periodic: break;
          case e_NotAKnot: break;
        }

      // determine i with m_Start[ coord::GetAxes().e_X]+(i-1)*m_Delta[ coord::GetAxes().e_X] < x < m_Start[ coord::GetAxes().e_X]+i*m_Delta[ coord::GetAxes().e_X] for the correct supporting points
      int i( int( floor( ( x - m_Start[ coord::GetAxes().e_X]) / m_Delta[ coord::GetAxes().e_X])));
      while( m_Start[ coord::GetAxes().e_X] + i * m_Delta[ coord::GetAxes().e_X] < x) { i++;}

      // determine j with m_Start[ coord::GetAxes().e_Y]+(j-1)*m_Delta[ coord::GetAxes().e_Y] < y < m_Start[ coord::GetAxes().e_Y]+j*m_Delta[ coord::GetAxes().e_Y] for the correct supporting points
      int j( int( floor( ( y - m_Start[ coord::GetAxes().e_Y]) / m_Delta[ coord::GetAxes().e_Y])));
      while( m_Start[ coord::GetAxes().e_Y] + j * m_Delta[ coord::GetAxes().e_Y] < y) { j++;}

      //see F(x, y) for a short explanation of the values
      const double delta_aktx( x-m_Start[ coord::GetAxes().e_X] - ( i - 1) * m_Delta[ coord::GetAxes().e_X]);
      const double delta_akty( y-m_Start[ coord::GetAxes().e_Y] - ( j - 1) * m_Delta[ coord::GetAxes().e_Y]);

      const double dxp( delta_aktx / m_Delta[ coord::GetAxes().e_X]);
      const double dxm( 1 - dxp);

      const double dyp( delta_akty / m_Delta[ coord::GetAxes().e_Y]);
      const double dym( 1 - dyp);
      const double dy3p( ( dyp * dyp * dyp - dyp) * m_DeltaNorm[ coord::GetAxes().e_Y]);
      const double dy3m( ( dym * dym * dym - dym) * m_DeltaNorm[ coord::GetAxes().e_Y]);

      //generate positive values to prevent some problems with the indices
      while( i < 1){ i += dimx;}
      while( j < 1){ j += dimy;}
      BCL_MessageCrt( "Final i orig: " + util::Format()( i));

      return
        -( dym * m_Values( ( i - 1) % dimx, ( j - 1) % dimy) + dyp * m_Values( ( i - 1) % dimx, j % dimy)) / m_Delta[ coord::GetAxes().e_X]
        +( dym * m_Values( i % dimx      , ( j - 1) % dimy) + dyp * m_Values( i % dimx        , j % dimy)) / m_Delta[ coord::GetAxes().e_X]
        - ( 3 * dxm * dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6 * ( dym * m_Dsecox( ( i - 1) % dimx, ( j - 1) % dimy) + dyp * m_Dsecox( ( i - 1) % dimx, j % dimy))
        + ( 3 * dxp * dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6 *( dym*m_Dsecox( i % dimx      , ( j - 1) % dimy) + dyp * m_Dsecox( i % dimx      , j%dimy))
        - ( dy3m * m_Dsecoy( ( i - 1) % dimx, ( j - 1) % dimy) + dy3p * m_Dsecoy( ( i - 1) % dimx, j % dimy)) / m_Delta[ coord::GetAxes().e_X]
        +( dy3m * m_Dsecoy( i % dimx       , ( j - 1) % dimy) + dy3p * m_Dsecoy( i % dimx     , j % dimy)) / m_Delta[ coord::GetAxes().e_X]
        - ( 3 * dxm * dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6 * ( dy3m * m_Dsecoxy( ( i - 1) % dimx, ( j - 1) % dimy) + dy3p * m_Dsecoxy( ( i - 1) % dimx, j % dimy))
        + ( 3 * dxp * dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6 *( dy3m*m_Dsecoxy( i % dimx    , ( j - 1) % dimy) + dy3p * m_Dsecoxy( i % dimx      , j % dimy))
      ;
    }

    //! return partial derivative at certain (x, y) for y
    double BicubicSpline::dFdy( const linal::Vector< double> &ARGUMENTS) const
    {
      BCL_Assert( ARGUMENTS.GetSize() == 2, "number of arguments doesn't match");

      const double x( ARGUMENTS( 0));
      const double y( ARGUMENTS( 1));

      const int dimx( m_Values.GetNumberRows());
      const int dimy( m_Values.GetNumberCols());

      //check if argument x is in range for non-periodic splines
      if( x < m_Start[ coord::GetAxes().e_X])
        switch( m_Border[ coord::GetAxes().e_X])
        {
          case e_FirstDer: BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return dFdy( linal::MakeVector( m_Start[ coord::GetAxes().e_X], y));
          case e_Natural : BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return dFdy( linal::MakeVector( m_Start[ coord::GetAxes().e_X], y));
          case e_Periodic: break;
          case e_NotAKnot: break;
        }

      if( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X] < x)
        switch( m_Border[ coord::GetAxes().e_X])
        {
          case e_FirstDer: BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return dFdy( linal::MakeVector( m_Start[ coord::GetAxes().e_X] + ( dimx-1 ) * m_Delta[ coord::GetAxes().e_X], y));
          case e_Natural : BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return dFdy( linal::MakeVector( m_Start[ coord::GetAxes().e_X] + ( dimx-1 ) * m_Delta[ coord::GetAxes().e_X], y));
          case e_Periodic: break;
          case e_NotAKnot: break;
        }

      //check if argument y is in range for non-periodic splines
      if( y < m_Start[ coord::GetAxes().e_Y])
        switch( m_Border[ coord::GetAxes().e_Y])
        {
          case e_FirstDer: BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return m_FirstBe[ coord::GetAxes().e_Y].First();
          case e_Natural : BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return dFdy( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y]));
          case e_Periodic: break;
          case e_NotAKnot: break;
        }

      if( m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y] < y)
        switch( m_Border[ coord::GetAxes().e_Y])
        {
          case e_FirstDer: BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return m_FirstBe[ coord::GetAxes().e_Y].Second();
          case e_Natural : BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return dFdy( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y] + ( dimy-1 ) * m_Delta[ coord::GetAxes().e_Y]));
          case e_Periodic: break;
          case e_NotAKnot: break;
        }

      //determine i with m_Start[ coord::GetAxes().e_X]+(i-1)*m_Delta[ coord::GetAxes().e_X] < x < m_Start[ coord::GetAxes().e_X]+i*m_Delta[ coord::GetAxes().e_X] for the correct supporting points
      int i( int( floor( ( x - m_Start[ coord::GetAxes().e_X]) / m_Delta[ coord::GetAxes().e_X])));
      while( m_Start[ coord::GetAxes().e_X] + i * m_Delta[ coord::GetAxes().e_X] < x) { i++;}
      if( !i){
        while( m_Start[ coord::GetAxes().e_X] + i * m_Delta[ coord::GetAxes().e_X] > x) { i--;}
        i++;
      }

      //determine j with m_Start[ coord::GetAxes().e_Y]+(j-1)*m_Delta[ coord::GetAxes().e_Y] < y < m_Start[ coord::GetAxes().e_Y]+j*m_Delta[ coord::GetAxes().e_Y] for the correct supporting points
      int j( int( floor( ( y - m_Start[ coord::GetAxes().e_Y]) / m_Delta[ coord::GetAxes().e_Y])));
      while( m_Start[ coord::GetAxes().e_Y] + j * m_Delta[ coord::GetAxes().e_Y] < y) { j++;}
      if( !j){
        while( m_Start[ coord::GetAxes().e_Y] + j * m_Delta[ coord::GetAxes().e_Y] > y) { j--;}
        j++;
      }

      //see F(x, y) for a short explanation of the values
      const double delta_aktx( x - m_Start[ coord::GetAxes().e_X] - ( i - 1) * m_Delta[ coord::GetAxes().e_X]);
      const double delta_akty( y - m_Start[ coord::GetAxes().e_Y] - ( j - 1) * m_Delta[ coord::GetAxes().e_Y]);

      const double dxp( delta_aktx / m_Delta[ coord::GetAxes().e_X]);
      const double dxm( 1 - dxp);
      const double dx3p( ( dxp * dxp * dxp - dxp) * m_DeltaNorm[ coord::GetAxes().e_X]);
      const double dx3m( ( dxm * dxm * dxm - dxm) * m_DeltaNorm[ coord::GetAxes().e_X]);

      const double dyp( delta_akty / m_Delta[ coord::GetAxes().e_Y]);
      const double dym( 1 - dyp);

      //generate positive values to prevent some problems with the indices
      while( i < 1){ i += dimx;}
      while( j < 1){ j += dimy;}

      return
          dxm *( -m_Values( ( i-1)%dimx  , (j-1)%dimy)+m_Values( (i-1)%dimx  , j%dimy))/m_Delta[ coord::GetAxes().e_Y]
        + dxp *( -m_Values( i%dimx      , (j-1)%dimy)+m_Values(i%dimx      , j%dimy))/m_Delta[ coord::GetAxes().e_Y]
        +dx3m *( -m_Dsecox( ( i-1)%dimx  , (j-1)%dimy)+m_Dsecox( (i-1)%dimx  , j%dimy))/m_Delta[ coord::GetAxes().e_Y]
        +dx3p *( -m_Dsecox( i%dimx      , (j-1)%dimy)+m_Dsecox(i%dimx      , j%dimy))/m_Delta[ coord::GetAxes().e_Y]
        + dxm *( -( 3 * dym * dym - 1) * m_Dsecoy( ( i-1)%dimx , ( j - 1)% dimy) +( 3 * dyp * dyp - 1) * m_Dsecoy( ( i-1)%dimx , j % dimy)) * m_Delta[ coord::GetAxes().e_Y]/ 6
        + dxp *( -( 3 * dym * dym - 1) * m_Dsecoy( i%dimx     , ( j - 1)% dimy) +( 3 * dyp * dyp - 1) * m_Dsecoy( i%dimx     , j % dimy)) * m_Delta[ coord::GetAxes().e_Y]/ 6
        +dx3m *( -( 3 * dym * dym - 1) * m_Dsecoxy( ( i-1)%dimx, ( j - 1)% dimy) +( 3 * dyp * dyp - 1) * m_Dsecoxy( ( i-1)%dimx, j % dimy)) * m_Delta[ coord::GetAxes().e_Y]/ 6
        +dx3p *( -( 3 * dym * dym - 1) * m_Dsecoxy( i%dimx    , ( j - 1)% dimy) +( 3 * dyp * dyp - 1) * m_Dsecoxy( i%dimx    , j % dimy)) * m_Delta[ coord::GetAxes().e_Y]/ 6
      ;
    }

    //!return value and derivative at certain (x, y)
    storage::Pair< double, linal::Vector< double> > BicubicSpline::FdF( const linal::Vector< double> &ARGUMENTS) const
    {
      BCL_Assert( ARGUMENTS.GetSize() == 2, "number of arguments doesn't match");

      double x = ARGUMENTS(0);
      double y = ARGUMENTS(1);

      int dimx = m_Values.GetNumberRows();
      int dimy = m_Values.GetNumberCols();

      //auxiliary variables for the function value and the derivatives
      double fvalue( 0), dfdxvalue( 0), dfdyvalue( 0);

      //check if argument is in range for non-periodic splines
      if( ( ( m_Border[ coord::GetAxes().e_X] != e_Periodic) && ( x < m_Start[ coord::GetAxes().e_X] || m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X] < x))
        || ( ( m_Border[ coord::GetAxes().e_Y] != e_Periodic) && ( y < m_Start[ coord::GetAxes().e_Y] || m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y] < y)))
      {
        if( x < m_Start[ coord::GetAxes().e_X])
        {
          BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");
          BCL_MessageDbg
          (
            "argument out of range (X < start : "
            + util::Format()( x)
            + " < "
            + util::Format()( m_Start[ coord::GetAxes().e_X])
            + "), using linear continuation"
          );
          fvalue    = F( linal::MakeVector( m_Start[ coord::GetAxes().e_X], y))+(x-m_Start[ coord::GetAxes().e_X])*dFdx( linal::MakeVector(m_Start[ coord::GetAxes().e_X], y));
          dfdxvalue = dFdx( linal::MakeVector( m_Start[ coord::GetAxes().e_X], y));
          dfdyvalue = dFdy( linal::MakeVector( m_Start[ coord::GetAxes().e_X], y));
        }
        if( x > m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X])
        {
          BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");
          BCL_MessageDbg
          (
            "argument out of range (X > end : "
            + util::Format()( x)
            + " > "
            + util::Format()( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X])
            + "), using linear continuation"
          );
          fvalue    = F(    linal::MakeVector( m_Start[ coord::GetAxes().e_X] + ( dimx-1 ) * m_Delta[ coord::GetAxes().e_X] , y))+(x-m_Start[ coord::GetAxes().e_X] - ( dimx-1 ) * m_Delta[ coord::GetAxes().e_X])*dFdx( linal::MakeVector( m_Start[ coord::GetAxes().e_X] + ( dimx-1 ) * m_Delta[ coord::GetAxes().e_X], y));
          dfdxvalue = dFdx( linal::MakeVector( m_Start[ coord::GetAxes().e_X] + ( dimx-1 ) * m_Delta[ coord::GetAxes().e_X] , y));
          dfdyvalue = dFdy( linal::MakeVector( m_Start[ coord::GetAxes().e_X] + ( dimx-1 ) * m_Delta[ coord::GetAxes().e_X] , y));
        }
        if( y < m_Start[ coord::GetAxes().e_Y])
        {
          BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");
          BCL_MessageDbg
          (
            "argument out of range (Y < start : "
            + util::Format()( y)
            + " < "
            + util::Format()( m_Start[ coord::GetAxes().e_Y])
            + "), using linear continuation"
          );
          fvalue    = F(    linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y]))+(y-m_Start[ coord::GetAxes().e_Y])*dFdy( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y]));
          dfdxvalue = dFdx( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y]));
          dfdyvalue = dFdy( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y]));
        }
        if( y > m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y])
        {
          BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");
          BCL_MessageDbg
          (
            "argument out of range (Y > end : "
            + util::Format()( y)
            + " > "
            + util::Format()( m_Start[ coord::GetAxes().e_Y] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_Y])
            + "), using linear continuation"
          );
          fvalue    = F(    linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y] + ( dimy-1 ) * m_Delta[ coord::GetAxes().e_Y]))+(y-m_Start[ coord::GetAxes().e_Y] - ( dimy-1 ) * m_Delta[ coord::GetAxes().e_Y])*dFdy( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y] + ( dimy-1 ) * m_Delta[ coord::GetAxes().e_Y]));
          dfdxvalue = dFdx( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y] + ( dimy-1 ) * m_Delta[ coord::GetAxes().e_Y]));
          dfdyvalue = dFdy( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y] + ( dimy-1 ) * m_Delta[ coord::GetAxes().e_Y]));
        }
      }else

      {
        //determine i with m_Start[ coord::GetAxes().e_X]+(i-1)*m_Delta[ coord::GetAxes().e_X] < x < m_Start[ coord::GetAxes().e_X]+i*m_Delta[ coord::GetAxes().e_X] for the correct supporting points
        int i( int( floor( ( x - m_Start[ coord::GetAxes().e_X]) / m_Delta[ coord::GetAxes().e_X])));
        while( m_Start[ coord::GetAxes().e_X] + i * m_Delta[ coord::GetAxes().e_X] < x)i++;
        if( !i) {
          while( m_Start[ coord::GetAxes().e_X] + i * m_Delta[ coord::GetAxes().e_X] > x)i--;
          i++;
        }

        //determine j with m_Start[ coord::GetAxes().e_Y]+(j-1)*m_Delta[ coord::GetAxes().e_Y] < y < m_Start[ coord::GetAxes().e_Y]+j*m_Delta[ coord::GetAxes().e_Y] for the correct supporting points
        int j( int( floor( ( y - m_Start[ coord::GetAxes().e_Y]) / m_Delta[ coord::GetAxes().e_Y])));
        while( m_Start[ coord::GetAxes().e_Y] + j * m_Delta[ coord::GetAxes().e_Y] < y)j++;
        if( !j) {
          while( m_Start[ coord::GetAxes().e_Y] + j * m_Delta[ coord::GetAxes().e_Y] > y)j--;
          j++;
        }

        //see method F(x,y) for detailed formula

        double delta_aktx = x-m_Start[ coord::GetAxes().e_X]-(i-1)*m_Delta[ coord::GetAxes().e_X];
        double delta_akty = y-m_Start[ coord::GetAxes().e_Y]-(j-1)*m_Delta[ coord::GetAxes().e_Y];

        double dxp( delta_aktx / m_Delta[ coord::GetAxes().e_X]);
        double dxm( 1 - dxp);
        double dx3p( ( dxp*dxp*dxp - dxp) * m_DeltaNorm[ coord::GetAxes().e_X]);
        double dx3m( ( dxm*dxm*dxm - dxm) * m_DeltaNorm[ coord::GetAxes().e_X]);

        double dyp( delta_akty / m_Delta[ coord::GetAxes().e_Y]);
        double dym( 1 - dyp);
        double dy3p( ( dyp*dyp*dyp - dyp) * m_DeltaNorm[ coord::GetAxes().e_Y]);
        double dy3m( ( dym*dym*dym - dym) * m_DeltaNorm[ coord::GetAxes().e_Y]);

        fvalue =

            dxm*(dym*m_Values( (i-1)%dimx  , (j-1)%dimy)+dyp*m_Values( (i-1)%dimx  , j%dimy))
          + dxp*(dym*m_Values(i%dimx      , (j-1)%dimy)+dyp*m_Values(i%dimx      , j%dimy))
          +dx3m*(dym*m_Dsecox( (i-1)%dimx  , (j-1)%dimy)+dyp*m_Dsecox( (i-1)%dimx  , j%dimy))
          +dx3p*(dym*m_Dsecox(i%dimx      , (j-1)%dimy)+dyp*m_Dsecox(i%dimx      , j%dimy))
          + dxm*(dy3m*m_Dsecoy( (i-1)%dimx , (j-1)%dimy)+dy3p*m_Dsecoy( (i-1)%dimx , j%dimy))
          + dxp*(dy3m*m_Dsecoy(i%dimx     , (j-1)%dimy)+dy3p*m_Dsecoy(i%dimx     , j%dimy))
          +dx3m*(dy3m*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy)+dy3p*m_Dsecoxy( (i-1)%dimx, j%dimy))
          +dx3p*(dy3m*m_Dsecoxy(i%dimx    , (j-1)%dimy)+dy3p*m_Dsecoxy(i%dimx    , j%dimy))
        ;

        dfdxvalue =
          -(dym*m_Values( (i-1)%dimx  , (j-1)%dimy)+dyp*m_Values( (i-1)%dimx  , j%dimy))/m_Delta[ coord::GetAxes().e_X]
          +(dym*m_Values(i%dimx      , (j-1)%dimy)+dyp*m_Values(i%dimx      , j%dimy))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecox( (i-1)%dimx  , (j-1)%dimy)+dyp*m_Dsecox( (i-1)%dimx  , j%dimy))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecox(i%dimx      , (j-1)%dimy)+dyp*m_Dsecox(i%dimx      , j%dimy))
          -(dy3m*m_Dsecoy( (i-1)%dimx , (j-1)%dimy)+dy3p*m_Dsecoy( (i-1)%dimx , j%dimy))/m_Delta[ coord::GetAxes().e_X]
          +(dy3m*m_Dsecoy(i%dimx     , (j-1)%dimy)+dy3p*m_Dsecoy(i%dimx     , j%dimy))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy)+dy3p*m_Dsecoxy( (i-1)%dimx, j%dimy))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxy(i%dimx    , (j-1)%dimy)+dy3p*m_Dsecoxy(i%dimx    , j%dimy))
        ;

        dfdyvalue =
          dxm*(-m_Values( (i-1)%dimx  , (j-1)%dimy)+m_Values( (i-1)%dimx  , j%dimy))/m_Delta[ coord::GetAxes().e_Y]
          + dxp*(-m_Values(i%dimx      , (j-1)%dimy)+m_Values(i%dimx      , j%dimy))/m_Delta[ coord::GetAxes().e_Y]
          +dx3m*(-m_Dsecox( (i-1)%dimx  , (j-1)%dimy)+m_Dsecox( (i-1)%dimx  , j%dimy))/m_Delta[ coord::GetAxes().e_Y]
          +dx3p*(-m_Dsecox(i%dimx      , (j-1)%dimy)+m_Dsecox(i%dimx      , j%dimy))/m_Delta[ coord::GetAxes().e_Y]
          + dxm*(-(3*dym*dym-1)*m_Dsecoy( (i-1)%dimx , (j-1)%dimy)+(3*dyp*dyp-1)*m_Dsecoy( (i-1)%dimx , j%dimy))* m_Delta[ coord::GetAxes().e_Y]/ 6
          + dxp*(-(3*dym*dym-1)*m_Dsecoy(i%dimx     , (j-1)%dimy)+(3*dyp*dyp-1)*m_Dsecoy(i%dimx     , j%dimy))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3m*(-(3*dym*dym-1)*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy)+(3*dyp*dyp-1)*m_Dsecoxy( (i-1)%dimx, j%dimy))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3p*(-(3*dym*dym-1)*m_Dsecoxy(i%dimx    , (j-1)%dimy)+(3*dyp*dyp-1)*m_Dsecoxy(i%dimx    , j%dimy))* m_Delta[ coord::GetAxes().e_Y]/ 6
        ;
      };

      double dfvalues[] = { dfdxvalue, dfdyvalue};
      linal::Vector< double> dfvector( 2, dfvalues);

      return storage::Pair< double, linal::Vector< double> >( fvalue, dfvector);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write BicubicSpline into std::ostream
    std::ostream &BicubicSpline::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write
      OSTREAM << m_Border[ coord::GetAxes().e_X]           << '\n';
      OSTREAM << m_Border[ coord::GetAxes().e_Y]           << '\n';
      OSTREAM << m_Start[ coord::GetAxes().e_X]            << '\n';
      OSTREAM << m_Delta[ coord::GetAxes().e_X]            << '\n';
      OSTREAM << m_Start[ coord::GetAxes().e_Y]            << '\n';
      OSTREAM << m_Delta[ coord::GetAxes().e_Y]            << '\n';
      OSTREAM << m_Values                   << '\n';
      OSTREAM << m_Dsecox                   << '\n';
      OSTREAM << m_Dsecoy                   << '\n';
      OSTREAM << m_Dsecoxy                  << '\n';
      OSTREAM << m_LinCont[ coord::GetAxes().e_X]          << '\n';
      OSTREAM << m_LinCont[ coord::GetAxes().e_Y]          << '\n';
      OSTREAM << m_FirstBe[ coord::GetAxes().e_X].First()  << '\n';
      OSTREAM << m_FirstBe[ coord::GetAxes().e_X].Second() << '\n';
      OSTREAM << m_FirstBe[ coord::GetAxes().e_Y].First()  << '\n';
      OSTREAM << m_FirstBe[ coord::GetAxes().e_Y].Second() << '\n';

      // end
      return OSTREAM;
    }

    //! read BicubicSpline from std::istream
    std::istream &BicubicSpline::Read( std::istream &ISTREAM)
    {
      // read parameters
      int border;
      ISTREAM >> border;
      m_Border[ coord::GetAxes().e_X] = SplineBorderType( border);
      ISTREAM >> border;
      m_Border[ coord::GetAxes().e_Y] = SplineBorderType( border);
      ISTREAM >> m_Start[ coord::GetAxes().e_X];
      ISTREAM >> m_Delta[ coord::GetAxes().e_X];
      ISTREAM >> m_Start[ coord::GetAxes().e_Y];
      ISTREAM >> m_Delta[ coord::GetAxes().e_Y];
      ISTREAM >> m_Values;
      ISTREAM >> m_Dsecox;
      ISTREAM >> m_Dsecoy;
      ISTREAM >> m_Dsecoxy;
      ISTREAM >> m_LinCont[ coord::GetAxes().e_X];
      ISTREAM >> m_LinCont[ coord::GetAxes().e_Y];
      ISTREAM >> m_FirstBe[ coord::GetAxes().e_X].First();
      ISTREAM >> m_FirstBe[ coord::GetAxes().e_X].Second();
      ISTREAM >> m_FirstBe[ coord::GetAxes().e_Y].First();
      ISTREAM >> m_FirstBe[ coord::GetAxes().e_Y].Second();
      m_DeltaNorm[ coord::GetAxes().e_X] = Sqr( m_Delta[ coord::GetAxes().e_X]) / 6.0;
      m_DeltaNorm[ coord::GetAxes().e_Y] = Sqr( m_Delta[ coord::GetAxes().e_Y]) / 6.0;

      // end
      return ISTREAM;
    }

  ////////////////
  // operations //
  ////////////////

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_comparisons.hpp"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Comparisons< double>;
    template class BCL_API Comparisons< float>;
    template class BCL_API Comparisons< int>;
    template class BCL_API Comparisons< unsigned int>;
    template class BCL_API Comparisons< unsigned long>;
    template class BCL_API Comparisons< unsigned long long>;
    template class BCL_API Comparisons< bool>;
    template class BCL_API Comparisons< char>;

  } // namespace math

  namespace util
  {

    template class BCL_API Enumerate< ShPtr< BinaryFunctionInterface<             double,             double, bool> >, math::Comparisons<             double> >;
    template class BCL_API Enumerate< ShPtr< BinaryFunctionInterface<              float,              float, bool> >, math::Comparisons<              float> >;
    template class BCL_API Enumerate< ShPtr< BinaryFunctionInterface<                int,                int, bool> >, math::Comparisons<                int> >;
    template class BCL_API Enumerate< ShPtr< BinaryFunctionInterface<       unsigned int,       unsigned int, bool> >, math::Comparisons<       unsigned int> >;
    template class BCL_API Enumerate< ShPtr< BinaryFunctionInterface<      unsigned long,      unsigned long, bool> >, math::Comparisons<      unsigned long> >;
    template class BCL_API Enumerate< ShPtr< BinaryFunctionInterface< unsigned long long, unsigned long long, bool> >, math::Comparisons< unsigned long long> >;
    template class BCL_API Enumerate< ShPtr< BinaryFunctionInterface<               bool,               bool, bool> >, math::Comparisons<               bool> >;
    template class BCL_API Enumerate< ShPtr< BinaryFunctionInterface<               char,               char, bool> >, math::Comparisons<               char> >;

  } // namespace util
} // namespace bcl
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
#include "math/bcl_math_contingency_matrix_measures.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

  //////////
  // data //
  //////////

    // separate, anonymous namespace to prevent exporting symbols
    namespace
    {
      // add each of the possible instances to the enumerated instances
      util::ObjectInterface *AddInstances()
      {
        // keep a pointer to the last created instance
        util::ObjectInterface *last_instance( NULL);
        for( size_t measure( 0); measure < ContingencyMatrixMeasures::s_NumberMeasures; ++measure)
        {
          last_instance =
            util::Enumerated< util::FunctionInterfaceSerializable< ContingencyMatrix, double> >::AddInstance
            (
              new ContingencyMatrixMeasures( static_cast< ContingencyMatrixMeasures::Measure>( measure))
            );
        }
        return last_instance;
      }
    }

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ContingencyMatrixMeasures::s_Instances( AddInstances());

    //! @brief Get the name, description, and function of the given measure
    //! @param MEASURE the measure of interest
    //! @return the short name or abbreviation for the class
    const storage::Triplet< std::string, std::string, function::MemberConst< ContingencyMatrix, double> > &
      ContingencyMatrixMeasures::GetMeasureInfo( const Measure &MEASURE)
    {
      typedef storage::Triplet< std::string, std::string, function::MemberConst< ContingencyMatrix, double> > t_Info;
      static const t_Info s_info[ s_NumberMeasures + 1] =
      {
        t_Info( "TPR", "True positive rate = TP / P", &math::ContingencyMatrix::GetTruePositiveRate),
        t_Info( "FPR", "False positive rate = FP / P", &math::ContingencyMatrix::GetFalsePositiveRate),
        t_Info( "FNR", "False negative rate = FN / N", &math::ContingencyMatrix::GetFalseNegativeRate),
        t_Info( "TNR", "True negative rate = TN / N", &math::ContingencyMatrix::GetTrueNegativeRate),
        t_Info( "Accuracy", "(TP + TN) / T", &math::ContingencyMatrix::GetAccuracy),
        t_Info( "Specificity", "= TNR = TN / N", &math::ContingencyMatrix::GetTrueNegativeRate),
        t_Info( "Precision", "= TP / PP (predicted positives)", &math::ContingencyMatrix::GetPrecision),
        t_Info( "Recall", "= TPR = TP / P", &math::ContingencyMatrix::GetRecall),
        t_Info( "PPV", "Positive predictive value = Precision = TP / PP (predicted positives)", &math::ContingencyMatrix::GetPositivePredictiveValue),
        t_Info( "Ideal-PPV", "Max possible PPV | precision = Min(1,P/PP).  Not appropriate for plotting against measures that consider FP or FN!", &math::ContingencyMatrix::GetIdealPositivePredictiveValue),
        t_Info( "Ideal-PPV_FPRelative", "P/(P+FP).  Appropriate for plotting against measures that consider FP or FN!", &math::ContingencyMatrix::GetIdealPositivePredictiveValueConsideringFalsePositives),
        t_Info( "NPV", "Negative predictive value = TN / PN (Predicted negatives)", &math::ContingencyMatrix::GetNegativePredictiveValue),
        t_Info( "FDR", "False discovery rate = 1 - precision = FP / PP", &math::ContingencyMatrix::GetFalseDiscoveryRate),
        t_Info( "HitRate", "fraction predicted positive = PP / T", &math::ContingencyMatrix::GetFractionPredictedPositives),
        t_Info( "MCC", "Matthews Correlation Coefficient = (TP * TN - FP * FN) / Sqrt( P * N * PP * PN)", &math::ContingencyMatrix::GetMatthewsCorrelationCoefficient),
        t_Info( "Enrichment", "Ratio of precision to positive rate = Precision / ( P / T)", &math::ContingencyMatrix::GetEnrichment),
        t_Info( "InformationGainRatio", "Information gain ratio, accounts for Positive/Negative split (range 0-1)", &math::ContingencyMatrix::GetInformationGainRatio),
        t_Info( "Cutoff", "The cutoff value at which the contingency matrix was created", NULL),
        t_Info( "LocalPPV", "Local-PPV. Best estimate of actual predictive value for predictions in a given range", NULL),
        t_Info( "", "", &math::ContingencyMatrix::GetTruePositiveRate)
      };
      return s_info[ MEASURE];
    };

    //! @brief Measure as string
    //! @param MEASURE the measure
    //! @return the string for the measure
    const std::string &ContingencyMatrixMeasures::GetMeasureName( const Measure &MEASURE)
    {
      return GetMeasureInfo( MEASURE).First();
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor sets everything to 0
    ContingencyMatrixMeasures::ContingencyMatrixMeasures( const Measure &MEASURE) :
      m_Measure( MEASURE),
      m_OptimizationDirection( true)
    {
      function::MemberConst< ContingencyMatrix, double> function( GetMeasureInfo( m_Measure).Third());
      if( function.IsDefined())
      {
        ContingencyMatrix good( 99, 1, 1, 99), bad( 1, 99, 99, 1);
        m_OptimizationDirection = function( good) > function( bad);
      }
    }

    //! @brief virtual copy constructor
    ContingencyMatrixMeasures *ContingencyMatrixMeasures::Clone() const
    {
      return new ContingencyMatrixMeasures( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ContingencyMatrixMeasures::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the measure that this contingency matrix calculator calculates
    const ContingencyMatrixMeasures::Measure &ContingencyMatrixMeasures::GetMeasure() const
    {
      return m_Measure;
    }

    //! @brief determine the optimization direction for this particular measure
    //! @return true if larger values for this metric are better
    const bool &ContingencyMatrixMeasures::GetOptimizationParity() const
    {
      return m_OptimizationDirection;
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &ContingencyMatrixMeasures::GetAlias() const
    {
      return GetMeasureInfo( m_Measure).First();
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator() taking a double based on the internal measure
    //! @param MATRIX contingency matrix of interest
    //! @return returns a double based on the selected measure
    double ContingencyMatrixMeasures::operator()( const ContingencyMatrix &MATRIX) const
    {
      return GetMeasureInfo( m_Measure).Third()( MATRIX);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ContingencyMatrixMeasures::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( GetMeasureInfo( m_Measure).Second());
      return serializer;
    }

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    const std::string &GetNamespaceIdentifier()
    {
      static const std::string *s_namespace_name( new std::string( util::ExtractNamespaceIdentifier( __PRETTY_FUNCTION__)));
      return *s_namespace_name;
    }

    //! returns X to the power of EXPONENT analogous to pow but for size_t
    // fastest possible x ^ n Russian peasant algorithm, Integer needs to be an integral type
    template<>
    size_t Pow< size_t>( const size_t &VALUE, const size_t &EXPONENT)
    {
      // initialize p to one
      size_t p( 1);

      // Create local variables to operate on the data
      size_t value( VALUE);
      size_t exponent( EXPONENT);

      if( p != value && exponent)
      {
        // if exponent is an odd number set p equal to value
        if( exponent & 1)
        {
          p = value;
        }
        // This loop computes the value
        while( exponent)
        {
          value *= value;
          exponent >>= 1;
          if( exponent & 1)
          {
            p *= value;
          }
        } // end while
      } // end if
      return p;
    }

    //! this function returns a weight according to an angle between between 0 .. pi :
    //! 0 .. pi/3 = 1; pi / 3 .. 2 * pi / 3: cosfunction from 1 to 0; 2 * pi / 3 .. 1 = 0
    double WeightBetweenZeroAndPi_ThreeSections( const double &ANGLE_RAD)
    {
      // lower boundary - return
      static const double lower_angle(     g_Pi / 3);
      static const double upper_angle( 2 * g_Pi / 3);

      //if smaller than lower angle return 1
      if( ANGLE_RAD <= lower_angle)
      {
        return double( 1);
      }

      //if larger than upper angle return 1
      if( ANGLE_RAD >= upper_angle)
      {
        return double( 0);
      }

      //otherwise apply sine transition
      return ( cos( ( ANGLE_RAD - lower_angle) * 3) + 1) / 2;
    }

    //! this function returns a weight according to an angle between between 0 .. pi : cosine function from 1 to 0;
    double WeightBetweenZeroAndPi( const double &ANGLE_RAD)
    {
      //otherwise apply sine transition
      return ( cos( ANGLE_RAD) + 1) / 2;
    }

    //! @brief compute the binomial coefficient n | k
    //! @param N size of set from which to choose
    //! @param K the number to choose
    //! @return the number of ways to choose K unique elements from a set of size N
    size_t BinomialCoefficient( const size_t &N, const size_t &K)
    {
      // handle the case where K > N-K
      if( K > N - K)
      {
        return BinomialCoefficient( N, N - K);
      }
      size_t combinations( 1);
      for( size_t numerator( N), denominator( 1); denominator <= K; ++denominator, --numerator)
      {
        combinations = combinations * numerator / denominator;
      }
      return combinations;
    }

    //! @brief compute the factorial
    //! @param N size of set from which to choose
    //! @return the number of permutations (orderings) of a set of size N = product( 1...N)
    uint64_t Factorial( const size_t &N)
    {
      // 21! is larger than 1 << 64, so only store 0-20!
      static const size_t max_factorial( 20);
      static const uint64_t s_factorials[ max_factorial + 1] =
      {
        UINT64_C( 1),
        UINT64_C( 1),
        UINT64_C( 2),
        UINT64_C( 6),
        UINT64_C( 24),
        UINT64_C( 120),
        UINT64_C( 720),
        UINT64_C( 5040),
        UINT64_C( 40320),
        UINT64_C( 362880),
        UINT64_C( 3628800),
        UINT64_C( 39916800),
        UINT64_C( 479001600),
        UINT64_C( 6227020800),
        UINT64_C( 87178291200),
        UINT64_C( 1307674368000),
        UINT64_C( 20922789888000),
        UINT64_C( 355687428096000),
        UINT64_C( 6402373705728000),
        UINT64_C( 121645100408832000),
        UINT64_C( 2432902008176640000)
      };
      if( N > max_factorial)
      {
        // return the max 64 bit integer
        return std::numeric_limits< uint64_t>::max();
      }
      return s_factorials[ N];
    }

    //! @brief error function
    //! http://en.wikipedia.org/wiki/Error_function
    //! http://homepages.physik.uni-muenchen.de/~Winitzki/erf-approx.pdf
    //! only accurate to better than 4 * 10**(-4) in relative error
    //! @param x argument to function
    //! @return the integral
    double Erf( const double x)
    {
      const double a( 0.140012288686666606); // approximation for the factor used in the calculation, further information see link above
      const double xsqr( Sqr( x)); // used frequently in the equation
      const double integral( Sqrt( 1.0 - exp( -xsqr * ( 4.0 / g_Pi + a * xsqr) / ( 1.0 + a * xsqr))));
      if( x >= 0.0)
      {
        return integral;
      }
      else//  if ( x < 0)
      {
        return -integral;
      }
    }

    //! @brief complimentary error function 1-erf(x)
    double Erfc( const double x)
    {
      return 1.0 - Erf( x);
    }

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_cubic_spline.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_inversion_gauss_jordan.h"
#include "math/bcl_math_smooth_data.h"
#include "storage/bcl_storage_pair.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> CubicSpline::s_Instance
    (
      GetObjectInstances().AddInstance( new CubicSpline())
    );

  ////////////////
  // operations //
  ////////////////

    //! @brief train CubicSpline
    //! @param BORDER determines the behavior of the spline at the borders (natural, first derivative, periodic)
    //! @param START the start of the interval the spline is defined on
    //! @param DELTA the distance between two support points of the spline
    //! @param RESULTS the function values at the support points of the spline
    //! @param FIRSTBE values for the first order derivative at begin and end of spline (only FIRSTDER)
    //! @return trained spline
    CubicSpline &CubicSpline::Train
    (
      const SplineBorderType BORDER,
      const double START,
      const double DELTA,
      const linal::VectorConstInterface< double> &RESULTS,
      const storage::Pair< double, double> &FIRSTBE
    )
    {
      // check, if the points are given in positive direction
      BCL_Assert( DELTA > 0, "deltax <= 0 not supported");

      // determine value for dimension of x
      const int dim( RESULTS.GetSize());

      // assigning values
      m_Border = BORDER;
      m_Start  = START;
      m_Delta  = DELTA;

      BCL_Assert( BORDER != e_NotAKnot, "The specified boundary condition is not yet supported");

      // auxiliary variables
      const double delta1( DELTA / 6), delta2( 2 * DELTA / 3);

      m_Values = RESULTS;

      linal::Matrix< double> coeffs( dim, dim);
      linal::Vector< double> derivs( dim, 0.0);

      // train once for the values of fxx
      // those values are equivalent for every type of spline considered here
      for( int i( 1); i < dim - 1; ++i)
      {
        coeffs( i, i - 1) = delta1;
        coeffs( i, i    ) = delta2;
        coeffs( i, i + 1) = delta1;

        derivs( i) = ( m_Values( i + 1) - 2 * m_Values( i) + m_Values( i - 1)) / DELTA;
      }

      // setting the second order derivative on start and end to 0, "natural" cubic spline
      if( m_Border == e_Natural)
      {
        coeffs(     0,     0) = 1;
        coeffs( dim-1, dim-1) = 1;
        // natural borders yields a simple tridiagonal matrix that can be readily solved using this specialty
        // function
        m_Dsecox = linal::MatrixInversionInterface< double>::SolveTridiagonalMatrix( coeffs, derivs);
        return *this;
      }

      // periodic, the function starts over after reaching the ends, continuously differentiable everywhere
      if( m_Border == e_Periodic)
      {
        coeffs(     0, dim - 1) = delta1;
        coeffs(     0,       0) = delta2;
        coeffs(     0,       1) = delta1;
        derivs(     0)          = ( m_Values( 1) - 2 * m_Values( 0) + m_Values( dim-1)) / DELTA;

        coeffs( dim - 1, dim - 2) = delta1;
        coeffs( dim - 1, dim - 1) = delta2;
        coeffs( dim - 1,       0) = delta1;
        derivs( dim - 1)          = ( m_Values( 0) - 2 * m_Values( dim-1) + m_Values( dim-2)) / DELTA;
      }

      // set the first order derivative at x_0 to first_start and at x_dim-1 to first_end
      if( m_Border == e_FirstDer)
      {
        coeffs(       0,       0) = -delta2/2;
        coeffs(       0,       1) = -delta1;
        derivs(       0)          = FIRSTBE.First() - ( m_Values( 1) - m_Values( 0)) / DELTA;

        coeffs( dim - 1, dim - 1) = delta2 / 2;
        coeffs( dim - 1, dim - 2) = delta1;
        derivs( dim - 1)          = FIRSTBE.Second() - ( m_Values( dim - 1) - m_Values( dim - 2)) / DELTA;
      }

      // computation of the second order derivatives in every given point of the spline
      // The spline matrix is always very well-conditioned, so pivoting is a waste of time, so skip pivoting
      // in the GJ solver by passing it false
      m_Dsecox = linal::MatrixInversionGaussJordan< double>( coeffs, false).Solve( derivs);

      return *this;
    }

    //! @brief train CubicSpline with pre processed data
    //! @param BORDER determines the behavior of the spline at the borders (natural, first derivative, periodic)
    //! @param START the start of the interval the spline is defined on
    //! @param DELTA the distance between two support points of the spline
    //! @param RESULTS the function values at the support points of the spline
    //! @param FIRSTBE values for the first order derivative at begin and end of spline (only FIRSTDER)
    //! @param PREPROC determines the degree of smoothing, should be between 0 and 1 to give the part/influence
    //! of the actual value to the pre processed value
    //! @return trained spline
    CubicSpline &CubicSpline::TrainWithPreprocessing
    (
      const SplineBorderType BORDER,
      const double START,
      const double DELTA,
      const linal::VectorConstInterface< double> &RESULTS,
      const storage::Pair< double, double> &FIRSTBE, const double PREPROC
    )
    {
      // passing a smoothed vector to the Train function
      return Train( BORDER, START, DELTA, SmoothData::SmoothVector( RESULTS, PREPROC, true), FIRSTBE);
    }

    //! @brief return derivative at certain ARGUMENT
    //! @param ARGUMENT x value
    //! @return derivative at ARGUMENT
    double CubicSpline::dF( const double &ARGUMENT) const
    {

      // number of grid points
      const int dim( m_Values.GetSize());

      int i( int( floor( ( ARGUMENT - m_Start) / m_Delta)) + 1);

      // not within supporting points - left
      if( i < 1)
      {
        // if the spline is periodic, adjust i to be positive ( > 0) and within range
        if( m_Border == e_Periodic)
        {
          // bring close to actual range
          i %= dim;

          // if between end and start
          if( i == 0)
          {
            // see Numerical recipes in C++, pages 116-118
            double dxp( fmod( ARGUMENT - m_Start, m_Delta) / m_Delta); // relative offset from the beginning of the actual interval
            if( dxp < 0.0)
            {
              dxp += 1.0;
            }

            // return
            return Derivative( dim - 1, 0, dxp);
          }
          // generate positive index value ( > 0)
          while( i < 1)
          {
            i += dim;
          }
        }
        else
        {
          return Derivative( 0, 1, double( 0));
        }
      }
      // not within supporting points - right
      else if( i >= dim)
      {
        const int end( dim - 1);

        // if the spline is periodic, adjust i to be positive ( > 0)
        if( m_Border == e_Periodic)
        {
          // generate index value within range
          i %= dim;

          // special case, where interpolation happens between end and beginning
          if( i == 0)
          {
            // see Numerical recipes in C++, pages 116-118
            const double dxp( fmod( ARGUMENT - m_Start, m_Delta) / m_Delta); // relative offset from the beginning of the actual interval

            // return
            return Derivative( end, 0, dxp);
          }
        }
        else
        {
          // derivative at last supporting points
          return Derivative( end - 1, end, 1.0);
        }
      }

      double fn_mod( fmod( ARGUMENT - m_Start, m_Delta));
      double dxp( fn_mod / m_Delta);

      // see Numerical recipes in C++, pages 116-118
      //double dxp( math::Fmod( ARGUMENT - m_Start, m_Delta) / m_Delta); // delta_akt is an offset from the beginning of the actual interval\n";

      // for numerical reason the fmod function sometimes give a value less than m_Delta
      // divison will have worked correctly, so need to reset dxp to 0.0.
      if( Absolute( fn_mod - m_Delta) < std::numeric_limits< double>::epsilon())
      {
        dxp = 0.0;
      }
      else if( dxp < 0.0)
      {
        dxp += 1.0;
      }

      return Derivative( i - 1, i, dxp);
    }

    //! @brief return derivative and value at certain ARGUMENT
    //! @param ARGUMENT x value
    //! @return value and derivative at ARGUMENT
    storage::Pair< double, double> CubicSpline::FdF( const double &ARGUMENT) const
    {
      return storage::Pair< double, double>( operator()( ARGUMENT), dF( ARGUMENT));
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a t_ResultType object
    //! @param ARGUMENT Argument to be used to evaluate the function
    //! @return function value of the given argument
    double CubicSpline::operator()( const double &ARGUMENT) const
    {
      // number of grid points
      const int dim( m_Values.GetSize());

      // determine i with m_Start+(i-1)*m_Delta < ARGUMENT < m_Start+i*m_Delta for the correct supporting points
      int i( int( floor( ( ARGUMENT - m_Start) / m_Delta)) + 1);

      // not within supporting points - left
      if( i < 1)
      {
        // if the spline is periodic, adjust i to be positive ( > 0) and within range
        if( m_Border == e_Periodic)
        {
          // bring close to actual range
          i %= dim;

          // if between end and start
          if( i == 0)
          {
            // see Numerical recipes in C++, pages 116-118
            double dxp( fmod( ARGUMENT - m_Start, m_Delta) / m_Delta); // relative offset from the beginning of the actual interval
            if( dxp < 0.0)
            {
              dxp += 1.0;
            }

            // return
            return Function( dim - 1, 0, dxp);
          }
          // generate positive index value ( > 0)
          while( i < 1)
          {
            i += dim;
          }
        }
        else
        {
          return Function( 0, 1, double( 0)) + ( ARGUMENT - m_Start) * Derivative( 0, 1, double( 0));
        }
      }
      // not within supporting points - right
      else if( i >= dim)
      {
        const int end( dim - 1);

        // if the spline is periodic, adjust i to be positive ( > 0)
        if( m_Border == e_Periodic)
        {
          // generate index value within range
          i %= dim;

          // special case, where interpolation happens between end and beginning
          if( i == 0)
          {
            // see Numerical recipes in C++, pages 116-118
            const double dxp( fmod( ARGUMENT - m_Start, m_Delta) / m_Delta); // relative offset from the beginning of the actual interval

            // return
            return Function( end, 0, dxp);
          }
        }
        else
        {
          // derivative at last supporting points
          return Function( end - 1, end, 1.0) + ( ARGUMENT - m_Start - ( dim - 1) * m_Delta) * Derivative( end - 1, end, 1.0);
        }
      }

      // see Numerical recipes in C++, pages 116-118
      double dxp( fmod( ARGUMENT - m_Start, m_Delta) / m_Delta); // delta_akt is an offset from the beginning of the actual interval\n";
      if( dxp < 0.0)
      {
        dxp += 1.0;
      }

      return Function( i - 1, i, dxp);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read CubicSpline from std::istream
    //! @param ISTREAM the stream to read the spline from
    //! @return istream after reading
    std::istream &CubicSpline::Read( std::istream &ISTREAM)
    {
      // read parameters
      int border;
      ISTREAM >> border;
      m_Border = SplineBorderType( border);
      io::Serialize::Read( m_Start , ISTREAM);
      io::Serialize::Read( m_Delta , ISTREAM);
      io::Serialize::Read( m_Values, ISTREAM);
      io::Serialize::Read( m_Dsecox, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write CubicSpline into std::ostream
    //! @param OSTREAM the stream to write the spline to
    //! @param INDENT indentation of the spline
    //! @return ostream after writing
    std::ostream &CubicSpline::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write parameters
      io::Serialize::Write( m_Border, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Start , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Delta , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Values, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Dsecox, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  ///////////////////////
  // helper  functions //
  ///////////////////////

    //! @brief calculate derivative between two cells
    //! @param INDEX_LEFT index of left grid point
    //! @param INDEX_RIGHT index of right grid point
    //! @param DXP relative distance from left grid point, must be element [0, 1]
    //! @return derivative depending on relative distance DXP
    double CubicSpline::Derivative( const int INDEX_LEFT, const int INDEX_RIGHT, const double DXP) const
    {
      // see Numerical recipes in C++, pages 116-118
      // relative distance from right grid point
      const double dxm( 1 - DXP);

      return
          ( m_Values( INDEX_RIGHT) - m_Values( INDEX_LEFT)) / m_Delta
        - ( 3 * dxm * dxm - 1) / 6 * m_Delta * m_Dsecox( INDEX_LEFT)
        + ( 3 * DXP * DXP - 1) / 6 * m_Delta * m_Dsecox( INDEX_RIGHT);
    }

    //! @brief calculate function between two cells
    //! @param INDEX_LEFT index of left grid point
    //! @param INDEX_RIGHT index of right grid point
    //! @param DXP relative distance from left grid point, must be element [0, 1]
    //! @return function depending on relative distance DXP
    double CubicSpline::Function( const int INDEX_LEFT, const int INDEX_RIGHT, const double DXP) const
    {
      // see Numerical recipes in C++, pages 116-118
      // relative distance from right grid point
      const double dxm( 1 - DXP);
      const double dx3p( ( DXP * DXP * DXP - DXP) * Sqr( m_Delta) / 6); // =0 at the gridpoints, adds cubic part of the spline
      const double dx3m( ( dxm * dxm * dxm - dxm) * Sqr( m_Delta) / 6); // =0 at the gridpoints, adds cubic part of the spline

      return
          dxm * m_Values( INDEX_LEFT) + DXP * m_Values( INDEX_RIGHT)
        + dx3m * m_Dsecox( INDEX_LEFT) + dx3p * m_Dsecox( INDEX_RIGHT);
    }

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_cubic_spline_damped.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_inversion_gauss_jordan.h"
#include "math/bcl_math_running_average.h"
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_logger_interface.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> CubicSplineDamped::s_Instance
    (
      GetObjectInstances().AddInstance( new CubicSplineDamped())
    );

  ////////////////
  // operations //
  ////////////////

    //! @brief train CubicSplineDamped
    //! @param X abscissa Values of dataset
    //! @param Y ordinate values of dataset
    //! @return trained spline
    CubicSplineDamped &CubicSplineDamped::Train
    (
      const linal::VectorConstInterface< double> &X,
      const linal::VectorConstInterface< double> &Y,
      const double &DY_START,
      const double &DY_END
    )
    {
      m_ConstantDelta = false;
      m_Delta = 0.0;
      BCL_Assert( X.GetSize() == Y.GetSize(), "Train called with non-equal sized vectors!");
      // test whether the x-range is unique.
      size_t n_duplicates( 0);
      if( X.GetSize() > size_t( 1))
      {
        double prev( X( 0));
        for( size_t i( 1), n_values( X.GetSize()); i < n_values; ++i)
        {
          if( X( i) == prev)
          {
            ++n_duplicates;
          }
          BCL_Assert( X( i) >= prev, "cubic spline requires sorted X vector! Given vector was: " + util::Format()( X));
          prev = X( i);
        }
      }
      if( !n_duplicates)
      {
        m_X = X;
        m_Y = Y;
      }
      else
      {
        // duplicates were found. Copy unique m_X, average Y across duplicate m_X
        const size_t n_uniq( X.GetSize() - n_duplicates);
        m_X = linal::Vector< double>( n_uniq);
        m_Y = linal::Vector< double>( n_uniq);
        double prev( std::numeric_limits< double>::quiet_NaN());
        for( size_t i( 0), n_values( X.GetSize()), uniq_i( 0); i < n_values;)
        {
          if( X( i) == prev)
          {
            // average Y across all shared X values
            RunningAverage< double> ave;
            ave += Y( uniq_i - 1);
            ave += Y( i);
            while( ++i < n_values && X( i) == prev)
            {
              ave += Y( i);
            }
            m_Y( uniq_i - 1) = ave.GetAverage();
          }
          else
          {
            prev = X( i);
            m_X( uniq_i) = X( i);
            m_Y( uniq_i) = Y( i);
            ++i;
            ++uniq_i;
          }
        }
      }
      m_Delta = ( m_X.Last() - m_X.First()) / double( std::max( m_X.GetSize(), size_t( 2)) - 1);

      // determine derivatives
      const size_t size( m_X.GetSize()), sizem1( size - 1);
      m_ConstantDelta = EqualWithinMachineTolerance( m_X, linal::FillVector( size, m_X.First(), m_Delta));
      m_dY = linal::Vector< double>( size, double( 0.0));

      if( size <= size_t( 1))
      {
        m_A = m_B = linal::Vector< double>();
        return *this;
      }
      m_A = m_B = linal::Vector< double>( sizem1, double( 0.0));

      linal::Vector< double> secants( m_dY), intervals( m_dY);
      for( size_t i( 0); i < sizem1; ++i)
      {
        intervals( i) = m_X( i + 1) - m_X( i);
        secants( i) = ( m_Y( i + 1) - m_Y( i)) / intervals( i);
      }

      // boundary condition; specified derivative at RHS
      if( size > 2)
      {
        intervals( size - 1) = intervals( size - 2);
        secants( size - 1) = util::IsDefined( DY_END) ? DY_END : secants( size - 2);
      }
      else
      {
        secants( size - 1) = util::IsDefined( DY_END) ? DY_END : 0.0;
      }
      for( size_t i( 1); i < size; ++i)
      {
        const bool going_down( secants( i - 1) <= 0.0);
        if( going_down != ( secants( i) <= 0.0))
        {
          // local maxima or minima in the input data / derivative must be 0
          continue;
        }
        // compute the secant through i-1 - i+1
        const double secant_ddy
        (
          ( secants( i - 1) * intervals( i) + secants( i) * intervals( i - 1))
          /
          ( intervals( i) + intervals( i - 1))
        );
        // derivative magnitude is the smaller absolute of the twice the adjacent secants, or the average secant
        m_dY( i) =
          (
            going_down
            ? std::max( secant_ddy, 2.0 * std::max( secants( i - 1), secants( i)))
            : std::min( secant_ddy, 2.0 * std::min( secants( i - 1), secants( i)))
          );
      }
      if( !util::IsDefined( DY_START))
      {
        // compute the secant through 0 - 1
        // equation 24
        const double secant_ddy
        (
          ( secants( 0) * intervals( 0) * ( 2.0 + intervals( 1) / intervals( 0)) - secants( 1))
          /
          ( intervals( 0) + intervals( 1))
        );
        const bool going_down( secants( 0) <= 0.0);
        if( going_down != ( secant_ddy < 0.0))
        {
          m_dY( 0) = 0.0;
        }
        else
        {
          m_dY( 0) =
          (
            going_down
            ? std::max( secant_ddy, 2.0 * secants( 0))
            : std::min( secant_ddy, 2.0 * secants( 0))
          );
        }
      }
      else
      {
        m_dY( 0) = DY_START;
      }
      if( !util::IsDefined( DY_END) && size > 2)
      {
        // equation 25
        const double interval_left( intervals( size - 3)), interval_right( intervals( size - 2));
        const double secants_final( secants( size - 2));
        const double secant_ddy
        (
          ( secants_final * interval_right * ( 2.0 + interval_left / interval_right) - secants( size - 3))
          /
          ( interval_right + interval_left)
        );
        const bool going_down( secants_final <= 0.0);
        if( going_down != ( secant_ddy < 0.0))
        {
          m_dY( size - 1) = 0.0;
        }
        else
        {
          m_dY( size - 1) =
          (
            going_down
            ? std::max( secant_ddy, 2.0 * secants_final)
            : std::min( secant_ddy, 2.0 * secants_final)
          );
        }
      }
      for( size_t i( 0); i < sizem1; ++i)
      {
        m_B( i) = ( 3.0 * secants( i) - 2.0 * m_dY( i) - m_dY( i + 1)) / intervals( i);
        m_A( i) = ( m_dY( i) + m_dY( i + 1) - 2.0 * secants( i)) / ( intervals( i) * intervals( i));
      }
      return *this;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator taking an value in the range of the trained cubic spline and returning what the function value
    //! @brief would be if there were a function
    //! @param ARGUMENT Argument to be used to evaluate the function
    //! @return function value of the given argument
    double CubicSplineDamped::operator()( const double &ARGUMENT) const
    {
      const size_t length( m_X.GetSize());
      if( length == size_t( 1))
      {
        // 1 point; no spline, assume constant value
        return m_Y( 0);
      }

      if( ARGUMENT <= m_X( 0))
      {
        return m_Y( 0) + m_dY( 0) * ( m_X( 0) - ARGUMENT);
      }
      else if( ARGUMENT >= m_X( length - 1))
      {
        return m_Y( length - 1) + m_dY( length - 1) * ( ARGUMENT - m_X( length - 1));
      }
      else if( !util::IsDefined( ARGUMENT))
      {
        return util::GetUndefined< double>();
      }

      size_t index( 0);

      if( !m_ConstantDelta)
      {
        // find the interval which contains ARGUMENT ( the sampled point)
        while( index < length && ARGUMENT > m_X( index))
        {
          ++index;
        }

        // handle case point as at a gridpoint
        if( index < length && ARGUMENT == m_X( index))
        {
          return m_Y( index);
        }

        --index;
        if( index >= m_X.GetSize())
        {
          BCL_MessageStd( "BAD INDEX A: X: " + util::Format()( ARGUMENT) + " " + util::Format()( m_X( 0)));
        }
      }
      else
      {
        // determine i with m_Start+(i-1)*m_Delta < ARGUMENT < m_Start+i*m_Delta for the correct supporting points
        index = int( floor( ( ARGUMENT - m_X( 0)) / m_Delta));
        if( index >= m_X.GetSize())
        {
          BCL_MessageStd( "BAD INDEX B: X: " + util::Format()( ARGUMENT) + " " + util::Format()( m_X( 0)));
        }
      }

      // Horner's rule to calculate cubic polynomial
      const double x( ARGUMENT - m_X( index));
      return ( ( m_A( index) * x + m_B( index)) * x + m_dY( index)) * x + m_Y( index);
    }

    //! @brief train CubicSpline, requiring a maximally smooth, piecewise monotonic solution
    //! @param START x-value for the start
    //! @param DELTA difference between consecutive X-points
    //! @param Y ordinate values of dataset
    //! @param DY_START first order derivative at X(0). if nan/undefined, will be the same as first order derivative at X(1)
    //! @param DY_END   first order derivative at X_last. if nan/undefined, will be the same as first order derivative at X(last-1)
    //! @return *this
    //! Reference: M. Steffen. "A simple method for monotonic interpolation in one dimension" Astronomy and
    //! Astrophysics. 1990
    //! Note that it is not required that the Y be monotonic. All this function does is ensure that each cubic piece
    //! of the function is monotonic over the given interval and keeps the y' continuous. y'' will not, in general, be
    //! continuous if the spline is trained with this function
    CubicSplineDamped &CubicSplineDamped::Train
    (
      const double &START,
      const double &DELTA,
      const linal::VectorConstInterface< double> &Y,
      const double &DY_START,
      const double &DY_END
    )
    {
      return this->Train( linal::FillVector( Y.GetSize(), START, DELTA), Y, DY_START, DY_END);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read CubicSplineDamped from std::istream
    //! @param ISTREAM the stream to read the spline from
    //! @return istream after reading
    std::istream &CubicSplineDamped::Read( std::istream &ISTREAM)
    {
      // read parameters
      io::Serialize::Read( m_X , ISTREAM);
      io::Serialize::Read( m_Y , ISTREAM);
      io::Serialize::Read( m_dY, ISTREAM);
      io::Serialize::Read( m_A, ISTREAM);
      io::Serialize::Read( m_B, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write CubicSplineDamped into std::ostream
    //! @param OSTREAM the stream to write the spline to
    //! @param INDENT indentation of the spline
    //! @return ostream after writing
    std::ostream &CubicSplineDamped::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write parameters
      io::Serialize::Write( m_X, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Y, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_dY, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_A, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_B, OSTREAM, INDENT) << '\n';

      // end
      return OSTREAM;
    }

  ///////////////////////
  // helper  functions //
  ///////////////////////

    //! @brief return derivative at certain ARGUMENT
    //! @param ARGUMENT x value
    //! @return derivative at ARGUMENT
    const double CubicSplineDamped::dF( const double &ARGUMENT) const
    {
      size_t length( m_X.GetSize());

      // if outside left boundary
      if( ARGUMENT < m_X( 0))
      {
        return m_dY( 0);
      }
      // if outside right boundary
      else if( ARGUMENT >= m_X( length - 1))
      {
        return m_dY( length - 1);
      }

      // locate left index
      size_t left_index( 0);
      if( !m_ConstantDelta)
      {
        for( size_t i( 0); m_X( i) <= ARGUMENT; ++i)
        {
          left_index = i;
        }
      }
      else
      {
        left_index = int( floor( ( ARGUMENT - m_X( 0)) / m_Delta));
      }

      // Set the offset from the beginning of the actual interval
      return Derivative( left_index, ( ARGUMENT - m_X( left_index)));
    }

    storage::Pair< double, double> CubicSplineDamped::FdF( const double &ARGUMENT) const
    {
      return storage::Pair< double, double>( operator()( ARGUMENT), dF( ARGUMENT));
    }

    //! @brief calculate derivative between two cells
    //! @param INDEX_LEFT index of left grid point
    //! @param DXP relative distance from left grid point, must be element [0, 1]
    //! @return derivative depending on relative distance DXP
    double CubicSplineDamped::Derivative( const int INDEX_LEFT, const double DXP) const
    {
      return m_dY( INDEX_LEFT) + DXP * ( 2.0 * m_B( INDEX_LEFT) + 3.0 * m_A( INDEX_LEFT) * DXP);
    }

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_cubic_spline_variable_delta.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_inversion_gauss_jordan.h"
#include "math/bcl_math_running_average.h"
#include "storage/bcl_storage_pair.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> CubicSplineVariableDelta::s_Instance
    (
      GetObjectInstances().AddInstance( new CubicSplineVariableDelta())
    );

  ////////////////
  // operations //
  ////////////////

    //! @brief train CubicSplineVariableDelta
    //! @param X abscissa Values of dataset
    //! @param Y ordinate values of dataset
    //! @return trained spline
    CubicSplineVariableDelta &CubicSplineVariableDelta::Train
    (
      const SplineBorderType BORDER,
      const linal::VectorConstInterface< double> &X,
      const linal::VectorConstInterface< double> &Y,
      const storage::Pair< double, double> &FIRSTBE
    )
    {
      m_Border = BORDER;

      // test whether the x-range is unique.
      size_t n_duplicates( 0);
      if( X.GetSize())
      {
        double prev( X( 0));
        for( size_t i( 1), n_values( X.GetSize()); i < n_values; ++i)
        {
          if( X( i) == prev)
          {
            ++n_duplicates;
          }
          BCL_Assert( X( i) >= prev, "cubic spline requires sorted X vector!");
          prev = X( i);
        }
      }
      if( !n_duplicates)
      {
        m_X = X;
        m_Y = Y;
      }
      else
      {
        // duplicates were found. Copy unique m_X, average Y across duplicate m_X
        const size_t n_uniq( X.GetSize() - n_duplicates);
        m_X = linal::Vector< double>( n_uniq);
        m_Y = linal::Vector< double>( n_uniq);
        double prev( std::numeric_limits< double>::quiet_NaN());
        for( size_t i( 0), n_values( X.GetSize()), uniq_i( 0); i < n_values;)
        {
          if( X( i) == prev)
          {
            // average Y across all shared X values
            RunningAverage< double> ave;
            ave += Y( uniq_i - 1);
            ave += Y( i);
            while( ++i < n_values && X( i) == prev)
            {
              ave += Y( i);
            }
            m_Y( uniq_i - 1) = ave.GetAverage();
          }
          else
          {
            prev = X( i);
            m_X( uniq_i) = X( i);
            m_Y( uniq_i) = Y( i);
            ++i;
            ++uniq_i;
          }
        }
      }

      BCL_Assert( BORDER == e_NotAKnot || BORDER == e_Natural, "The specified boundary condition is not yet supported");

      // Compute 2nd derivatives at each point
      GenerateSecondDerivative();

      return *this;
    }

    //! @brief train CubicSpline for Constant Delta.  This overloads the train function with the same input from the
    //         original implementation of cublic spline
    //! @param BORDER determines the behavior of the spline at the borders (natural, first derivative, periodic)
    //! @param START the start of the interval the spline is defined on
    //! @param DELTA the distance between two support points of the spline
    //! @param RESULTS the function values at the support points of the spline
    //! @param FIRSTBE values for the first order derivative at begin and end of spline (only FIRSTDER)
    //! @return trained spline
    CubicSplineVariableDelta &CubicSplineVariableDelta::Train
    (
      const SplineBorderType BORDER,
      const double START,
      const double DELTA,
      const linal::Vector< double> &RESULTS,
      const storage::Pair< double, double> &FIRSTBE
    )
    {
      m_X = linal::FillVector( RESULTS.GetSize(), START, DELTA);
      m_Y = RESULTS;
      m_Border = BORDER;

      // Compute 2nd derivatives at each point
      GenerateSecondDerivative();

      return *this;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator taking an value in the range of the trained cubic spline and returning what the function value
    //! @brief would be if there were a function
    //! @param ARGUMENT Argument to be used to evaluate the function
    //! @return function value of the given argument
    double CubicSplineVariableDelta::operator()( const double &ARGUMENT) const
    {
      const size_t length( m_X.GetSize());
      if( length == size_t( 1))
      {
        // 1 point; no spline, assume constant value
        return m_Y( 0);
      }

      // !Todo should we extend the spline to extrapolation in addition to interpolation
      BCL_Assert( ARGUMENT >= m_X( 0) && ARGUMENT <= m_X( length - 1), "The input is outside the trained spline");

      double A0, A1, A2, A3;
      size_t index( 0);

      // find the interval which contains ARGUMENT ( the sampled point)
      while( index < length && ARGUMENT > m_X( index))
      {
        ++index;
      }

      // handle case point as at a gridpoint. This is necessary only due to the bounds assertion above
      if( index < length && ARGUMENT == m_X( index))
      {
        return m_Y( index);
      }

      --index;

      const double dx( m_X( index + 1) - m_X( index));
      const double dy( m_Y( index + 1) - m_Y( index));

      A0 = m_Y( index);
      A1 = dy / dx - ( dx / 6.0) * ( m_SecondDerivative( index + 1) + 2.0 * m_SecondDerivative( index));
      A2 = m_SecondDerivative( index) / 2.0;
      A3 = ( m_SecondDerivative( index + 1) - m_SecondDerivative( index)) / ( 6.0 * dx);

      // Horner's rule to calculate cubic polynomial
      const double x( ARGUMENT - m_X( index));
      return ( ( A3 * x + A2) * x + A1) * x + A0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read CubicSplineVariableDelta from std::istream
    //! @param ISTREAM the stream to read the spline from
    //! @return istream after reading
    std::istream &CubicSplineVariableDelta::Read( std::istream &ISTREAM)
    {
      // read parameters
      io::Serialize::Read( m_X , ISTREAM);
      io::Serialize::Read( m_Y , ISTREAM);
      io::Serialize::Read( m_SecondDerivative, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write CubicSplineVariableDelta into std::ostream
    //! @param OSTREAM the stream to write the spline to
    //! @param INDENT indentation of the spline
    //! @return ostream after writing
    std::ostream &CubicSplineVariableDelta::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write parameters
      io::Serialize::Write( m_X, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Y, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SecondDerivative, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  ///////////////////////
  // helper  functions //
  ///////////////////////

    //! @brief use m_x and m_y and initinal yd vectors to set m_SecondDerivative
    void CubicSplineVariableDelta::GenerateSecondDerivative()
    {
      size_t vector_size( m_X.GetSize());
      if( vector_size < size_t( 3))
      {
        // handle 0-2 points
        m_SecondDerivative = linal::Vector< double>( vector_size, double( 0.0));
        return;
      }

      linal::Vector< double> yd( vector_size);
      double h0, h1, r0, r1;

      // initialize first row data

      // delta x0
      h0 = m_X( 1) - m_X( 0);

      // delta x1
      h1 = m_X( 2) - m_X( 1);

      // delta y0
      r0 = ( m_Y( 1) - m_Y( 0)) / h0;

      // delta y1
      r1 = ( m_Y( 2) - m_Y( 1)) / h1;

      // This matrix will be filled with the coefficients of the derivative equations
      linal::Matrix< double> coeffs( vector_size, vector_size);

      if( m_Border == e_NotAKnot)
      {
        coeffs( 0, 0) = -h1;
        coeffs( 0, 1) = h0 + h1;
        coeffs( 0, 2) = -h0;
        yd( 0) = ( 0.0);
      }
      else // ( m_Border == e_Natural)
      {
        // Set B Value for initial point
         coeffs( 0, 0) = 1.0;
         yd( 0) = ( 0.0);
      }

      // initialize interior data
      for( size_t i( 1); i < vector_size - 1; i++)
      {
        h0 = m_X( i) - m_X( i - 1);
        h1 = m_X( i + 1) - m_X( i);
        r0 = ( m_Y( i) - m_Y( i - 1)) / h0;
        r1 = ( m_Y( i + 1) - m_Y( i)) / h1;

        // Set A Value
        coeffs( i, i - 1) = h0;

        // Set B Value
        coeffs( i, i) = 2 * ( h0 + h1);

        // Set C Value
        coeffs( i, i + 1) = h1;

        yd( i) = 6 * ( r1 - r0);
      }

      // initialize last row data

      if( m_Border == e_NotAKnot)
      {
        // Set A Value
        coeffs( vector_size - 1, vector_size - 3) = -h1;

        // Set B Value
        coeffs( vector_size - 1, vector_size - 2) = h0 + h1;

        // Set C Value
        coeffs( vector_size - 1, vector_size - 1) = -h0;

        yd( vector_size - 1) = 0.0;
        m_SecondDerivative = linal::MatrixInversionGaussJordan< double>( coeffs, false).Solve( yd);
      }
      else // if ( m_Border == e_Natural)
      {
        // Set C Value
        coeffs( vector_size - 1, vector_size - 1) = 1.0;
        yd( vector_size - 1) = 0.0;

        // computation of the second order derivatives in every given point of the spline
        // natural borders yields a simple tridiagonal matrix that can be readily solved using this specialty
        // function
        m_SecondDerivative = linal::MatrixInversionInterface< double>::SolveTridiagonalMatrix( coeffs, yd);
      }
    }

    //! @brief return derivative at certain ARGUMENT
    //! @param ARGUMENT x value
    //! @return derivative at ARGUMENT
    const double CubicSplineVariableDelta::dF( const double &ARGUMENT) const
    {

      size_t length( m_X.GetSize());

      // if outside left boundary
      if( ARGUMENT < m_X( 0))
      {
        return Derivative( 0, 1, double( 0.0));
      }
      // if outside right boundary
      else if( ARGUMENT >= m_X( length - 1))
      {
        const int end( length - 1);
        return Derivative( end - 1, end, 1.0);
      }

      // locate left index
      size_t i( 0);
      int left_index( 0);

      while( m_X( i) <= ARGUMENT)
      {
        left_index = i;
        ++i;
      }

      // Assign right index
      int right_index( left_index + 1);

      double left_border( m_X( left_index));
      double right_border( m_X( right_index));
      double delta( right_border - left_border);

      // Set the offset from the beginning of the actual interval
      double dxp( ( ARGUMENT - m_X( left_index)) / delta);

      return Derivative( left_index, right_index, dxp);
    }

    storage::Pair< double, double> CubicSplineVariableDelta::FdF( const double &ARGUMENT) const
    {
      return storage::Pair< double, double>( operator()( ARGUMENT), dF( ARGUMENT));
    }

    //! @brief calculate derivative between two cells
    //! @param INDEX_LEFT index of left grid point
    //! @param INDEX_RIGHT index of right grid point
    //! @param DXP relative distance from left grid point, must be element [0, 1]
    //! @return derivative depending on relative distance DXP
    double CubicSplineVariableDelta::Derivative( const int INDEX_LEFT, const int INDEX_RIGHT, const double DXP) const
    {
      const double dxm( 1 - DXP);
      const double delta( m_X( INDEX_RIGHT) - m_X( INDEX_LEFT));

      return
          ( m_Y( INDEX_RIGHT) - m_Y( INDEX_LEFT)) / delta
        - ( 3.0 * dxm * dxm - 1) / 6.0 * delta * m_SecondDerivative( INDEX_LEFT)
        + ( 3.0 * DXP * DXP - 1) / 6.0 * delta * m_SecondDerivative( INDEX_RIGHT);
    }

  } // namespace math
} // namespace bcl
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

// include for this class
#include <math/bcl_math_discrete_set_selector.h>
// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_polynomial.h"
#include "random/bcl_random_uniform_distribution.h"

namespace bcl
{
  namespace math
  {

    //! @brief default constructor; yields uniform probability distribution
    DiscreteSetSelector::DiscreteSetSelector() :
      m_Weights(),
      m_ProbabilityFunction()
    {
      Polynomial constant_polynomial;
      linal::Vector< double> const_vector( 1);
      const_vector( 0) = 1.0;
      constant_polynomial.SetCoefficients( const_vector);
      m_ProbabilityFunction = util::CloneToShPtr( constant_polynomial);
    }

    //! @brief constructor with a probability function and intial values
    //! @param PROBABILITY_FUNCTION the probability function to use to select from values given to the object
    //! @param VALUES the values to add initially
    DiscreteSetSelector::DiscreteSetSelector
    (
      const FunctionInterfaceSerializable< double, double> &PROBABILITY_FUNCTION,
      const linal::VectorConstInterface< double> &VALUES
    ) :
      m_Weights(),
      m_ProbabilityFunction( util::CloneToShPtr( PROBABILITY_FUNCTION))
    {
      Prepare( VALUES);
    }

    //! @brief copy constructor
    DiscreteSetSelector::DiscreteSetSelector( const DiscreteSetSelector &OTHER) :
      m_Weights( OTHER.m_Weights),
      m_ProbabilityFunction( util::CloneToShPtr( *OTHER.m_ProbabilityFunction))
    {
    }

    //! @brief assignment operator
    DiscreteSetSelector &DiscreteSetSelector::operator =( const DiscreteSetSelector &OTHER)
    {
      m_Weights = OTHER.m_Weights;
      m_ProbabilityFunction = util::CloneToShPtr( *OTHER.m_ProbabilityFunction);
      return *this;
    }

    //! @brief clone function
    //! @return a copy of this class
    DiscreteSetSelector *DiscreteSetSelector::Clone() const
    {
      return new DiscreteSetSelector( *this);
    }

    //! @brief the name of this class
    //! @return the name of the class
    const std::string &DiscreteSetSelector::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief set up the selector using the given values (resets everything beforehand)
    //! @param @VALUES a vector of states that will be converted into weights
    void DiscreteSetSelector::Prepare( const linal::VectorConstInterface< double> &VALUES)
    {
      BCL_Assert
      (
        m_ProbabilityFunction.IsDefined(),
        "Prepare was called before DiscreteSetSelector probability function was defined"
      );

      size_t n_values( VALUES.GetSize());
      m_Weights = linal::Vector< double>( n_values);
      // VALUES is empty, do nothing
      if( !n_values)
      {
        return;
      }

      // Initially set the weights vector to the raw probabilities, then scale it accoring to the sum of the values that
      // were given
      for( size_t i( 0); i < n_values; ++i)
      {
        m_Weights( i) = ( *m_ProbabilityFunction)( VALUES( i));
      }

      // if all values are zero then this should be a normalized uniform distribution, otherwise just normalize the weights
      double weights_sum( m_Weights.Sum());
      if( weights_sum == 0.0)
      {
        m_Weights = 1.0 / n_values;
      }
      else
      {
        m_Weights /= weights_sum;
      }
    }

    //! @brief Reset the class
    void DiscreteSetSelector::Reset()
    {
      m_Weights = linal::Vector< double>();
    }

    //! @brief operator() selects a state from the distrubution
    //! @param ARGUMENT a number between 0.0 and 1.0 that is used for member selection
    //! @return the index of the state that was chosen, or undefined if there were no weights selected
    size_t DiscreteSetSelector::operator()( const double &ARGUMENT) const
    {
      size_t result( util::GetUndefined< size_t>());

      double arg( ARGUMENT);
      if( !util::IsDefined( arg))
      {
        arg = random::GetGlobalRandom().Double();
      }

      if( m_Weights.GetSize())
      {
        // The objective sum should be between 0 and 1
        double sum( std::max( std::min( arg, 1.0), 0.0));

        // Add weight values until they exceed the objective sum
        for( result = 0; result < m_Weights.GetSize(); ++result)
        {
          sum -= m_Weights( result);
          if( sum <= 0.0)
          {
            break;
          }
        }
        result = std::min( result, m_Weights.GetSize() - 1);
      }
      return result;
    }

    //! @brief get the weights vector
    //! @return linal::vector containing weights
    const linal::Vector< double> &DiscreteSetSelector::GetWeights() const
    {
      return m_Weights;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DiscreteSetSelector::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return output stream which was written to
    std::ostream &DiscreteSetSelector::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_divide_equals.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API DivideEquals< double>;
    template class BCL_API DivideEquals< float>;

  } // namespace math

} // namespace bcl
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
#include "math/bcl_math_gaussian_function.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_range.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> GaussianFunction::s_Instance
    (
      util::Enumerated< FunctionInterfaceSerializable< double, double> >::AddInstance( new GaussianFunction())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    GaussianFunction::GaussianFunction() :
      m_Mean( 1),
      m_Sigma( 1),
      m_OneOverTwoSquareSigma( 1 / ( 2 * Sqr( 1)))
    {
    }

    //! @brief construct from mean and sigma
    GaussianFunction::GaussianFunction
    (
      const double MEAN,
      const double SIGMA
    ) :
      m_Mean( MEAN),
      m_Sigma( SIGMA),
      m_OneOverTwoSquareSigma( 1.0 / ( 2.0 * Sqr( SIGMA)))
    {
    }

    //! @brief Clone function
    //! @return pointer to new GaussianFunction
    GaussianFunction *GaussianFunction::Clone() const
    {
      return new GaussianFunction( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &GaussianFunction::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief access to mean
    //! @return the mean of the gaussian
    double GaussianFunction::GetMean() const
    {
      return m_Mean;
    }

    //! @brief access to sigma
    //! @return the sigma of the guassian
    double GaussianFunction::GetSigma() const
    {
      return m_Sigma;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &GaussianFunction::GetAlias() const
    {
      static const std::string s_Name( "GaussianFunction");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get symmetric range relative to the mean for a given multiple of sigma
    //! @param MULTIPLE_SIGMA n multiple of sigma
    //! @return Range[ mean - n * sigma, mean + n * sigma]
    Range< double> GaussianFunction::GetRange( const size_t MULTIPLE_SIGMA) const
    {
      const double half_width( MULTIPLE_SIGMA * m_Sigma);
      return Range< double>( m_Mean - half_width, m_Mean + half_width);
    }

    //! @brief get normalization factor
    //! @return a factor so that integral f(x) -inf to +inf = 1
    double GaussianFunction::GetNormalizationParameter() const
    {
      return double( 1) / ( m_Sigma * Sqrt( 2 * g_Pi));
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator f( x) = y = e^( -(x-mean)^2 / ( 2*sigma^2))
    //! @param ARGUMENT the x
    //! @return f( x)
    double GaussianFunction::operator()( const double &ARGUMENT) const
    {
      return exp( -( Sqr( ARGUMENT - m_Mean) * m_OneOverTwoSquareSigma));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &GaussianFunction::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Mean, ISTREAM);
      io::Serialize::Read( m_Sigma, ISTREAM);

      // calculate members
      m_OneOverTwoSquareSigma = 1 / ( 2 * Sqr( m_Sigma));

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &GaussianFunction::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Mean, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_Sigma, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer GaussianFunction::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Represents a Gaussian function. A function of the form f(x) = exp( -( x - mean)^2 / ( 2*sigma^2))"
      );

      parameters.AddInitializer
      (
        "mean",
        "mean of the gaussian distribution",
        io::Serialization::GetAgent( &m_Mean),
        "1"
      );

      parameters.AddInitializer
      (
        "sigma",
        "standard deviation or sigma of the distribution",
        io::Serialization::GetAgent( &m_Sigma),
        "1"
      );

      return parameters;
    }

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    bool GaussianFunction::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      m_OneOverTwoSquareSigma = double( 1.0) / double( 2.0 * m_Sigma * m_Sigma);
      return true;
    }

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_gnuplot.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "coord/bcl_coord_axes.h"
#include "linal/bcl_linal_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    //! @brief flag for setting number of pixels in x direction
    //! @return flag for setting number of pixels in x direction
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagXPixels()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "gnuplot_pixels_x",
          "The number of pixels in x direction.",
          command::Parameter( "number_pixels", "size_t which is the number of pixels", "900")
        )
      );

      return s_flag;
    }

    //! @brief flag for setting number of pixels in y direction
    //! @return flag for setting number of pixels in y direction
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagYPixels()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "gnuplot_pixels_y",
          "The number of pixels in y direction.",
          command::Parameter( "number_pixels", "size_t which is the number of pixels", "600")
        )
      );

      return s_flag;
    }

    //! @brief flag for setting the title of the plot
    //! @return flag for setting the title of the plot
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagTitle()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "gnuplot_title",
          "The title of the gnuplot.",
          command::Parameter( "title", "string which is the title of the plot", "YourPlotTitleHere")
        )
      );

      return s_flag;
    }

    //! @brief flag for determining whether the key is shown or not in the plot
    //! @return flag for determining whether the key is shown or not in the plot
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagSetKey()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "gnuplot_set_key",
          "Determines whether the key is shown or not in the plot."
        )
      );

      return s_flag;
    }

    //! @brief flag for setting the label on the x axis
    //! @return flag for setting the label on the x axis
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagXLabel()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "gnuplot_x_label",
          "The label of the x axis.",
          command::Parameter( "label", "string which is the label of the x axis", "x-axis")
        )
      );

      return s_flag;
    }

    //! @brief flag for setting the label on the y axis
    //! @return flag for setting the label on the y axis
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagYLabel()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "gnuplot_y_label",
          "The label of the y axis.",
          command::Parameter( "label", "string which is the label of the y axis", "y-axis")
        )
      );

      return s_flag;
    }

    //! @brief flag for setting the maximum value on the x axis
    //! @return flag for setting the maximum value on the x axis
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagXMax()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "gnuplot_max_x",
          "The maximum value shown on the x axis.",
          command::Parameter( "max_value", "double which is the maximum value shown on the x axis", "25.0")
        )
      );

      return s_flag;
    }

    //! @brief flag for setting the maximum value on the y axis
    //! @return flag for setting the maximum value on the y axis
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagYMax()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "gnuplot_max_y",
          "The maximum value shown on the y axis.",
          command::Parameter( "max_value", "double which is the maximum value shown on the y axis", "25.0")
        )
      );

      return s_flag;
    }

    //! @brief flag for setting the minimum value on the x axis
    //! @return flag for setting the minimum value on the x axis
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagXMin()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "gnuplot_min_x",
          "The minimum value shown on the x axis.",
          command::Parameter( "min_value", "double which is the minimum value shown on the x axis", "0.0")
        )
      );

      return s_flag;
    }

    //! @brief flag for setting the minimum value on the y axis
    //! @return flag for setting the minimum value on the y axis
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagYMin()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "gnuplot_min_y",
          "The minimum value shown on the y axis.",
          command::Parameter( "min_value", "double which is the minimum value shown on the y axis", "0.0")
        )
      );

      return s_flag;
    }

    //! @brief flag for setting series names to be plotted
    //! @return flag for setting series names to be plotted
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagSeriesNames()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagDynamic
        (
          "gnuplot_series_names",
          "The series names to be plotted. Zero or more can be passed",
          command::Parameter( "series_name", "string which is a series names to be plotted", "series_a")
        )
      );

      return s_flag;
    }

    //! @brief flag for specifying filenames of tables to be included in distribution plot
    //! @return flag for specifying filenames of tables to be included in distribution plot
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagTableInputFilenames()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagDynamic
        (
          "gnuplot_input_table_filenames",
          "The filenames of tables to be included in distribution plot. Zero or more can be passed",
          command::Parameter( "table_filenames", "string which filename of a table to be plotted", "table_a")
        )
      );

      return s_flag;
    }

    //! @brief flag for specifying columns of the tables to be used for the plot by column name for x axis
    //! @return flag for specifying columns of the tables to be used for the plot by column name for x axis
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagTableColumnsX()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagDynamic
        (
          "gnuplot_table_columns_x",
          "The columns of the tables to be used for the plot by column name for x axis. Zero or more can be passed."
          " If more than one is passed they will be summed to create one value per row (datapoint).",
          command::Parameter( "column_name", "string which is the column names for plotting", "RMSD100")
        )
      );

      return s_flag;
    }

    //! @brief flag for specifying columns of the tables to be used for the plot by column name for y axis
    //! @return flag for specifying columns of the tables to be used for the plot by column name for y axis
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagTableColumnsY()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagDynamic
        (
          "gnuplot_table_columns_y",
          "The columns of the tables to be used for the plot by column name for y axis. Zero or more can be passed."
          " If more than one is passed they will be summed to create one value per row (datapoint).",
          command::Parameter( "column_name", "string which is the column names for plotting", "RMSD100")
        )
      );

      return s_flag;
    }

    //! @brief flag for specifying the name of the outputted gnuplot file
    //! @return flag for specifying the name of the outputted gnuplot file
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagGnuplotOutputFilename()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "gnuplot_output_filename",
          "The name of the outputted gnuplot file.",
          command::Parameter
          (
            "output_filename", "string which is the name of the outputted gnuplot file", "output.gnuplot"
          )
        )
      );

      return s_flag;
    }

    //! @brief flag for specifying the interval of tics on the x axis
    //! @return flag for specifying the interval of tics on the x axis
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagTicIntervalX()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "gnuplot_tics_x",
          "The interval of tics on the x axis.",
          command::Parameter( "tick_interval", "size_t which is the interval of tics on the x axis", "1")
        )
      );

      return s_flag;
    }

    //! @brief flag for specifying the interval of tics on the y axis
    //! @return flag for specifying the interval of tics on the y axis
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagTicIntervalY()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "gnuplot_tics_y",
          "The interval of tics on the y axis.",
          command::Parameter( "tick_interval", "size_t which is the interval of tics on the y axis", "1")
        )
      );

      return s_flag;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Gnuplot::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief string for axis
    //! @param AXIS the axis
    //! @return the string for the axis
    const std::string &Gnuplot::GetAxisString( const coord::Axis &AXIS)
    {
      static const std::string s_axis_string[] = { "x", "y", "cb"};
      return s_axis_string[ AXIS.GetIndex()];
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief tics from a vector of binning
    //! @param BINNING vector of bin centers
    //! @param BINS_PER_TIC number of bins that should be associated with one tic
    //! @param FORMAT format converting double to string
    //! @return vector of tic labels
    storage::Vector< std::string> Gnuplot::TicsFromBinning
    (
      const linal::Vector< double> &BINNING,
      const size_t BINS_PER_TIC,
      const util::Format &FORMAT
    )
    {
      storage::Vector< std::string> tics;
      for( const double *ptr( BINNING.Begin()), *ptr_end( BINNING.End()); ptr < ptr_end; ptr += BINS_PER_TIC)
      {
        tics.PushBack( FORMAT( *ptr));
      }

      // end
      return tics;
    }

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_gnuplot_heatmap.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"
#include "math/bcl_math_bicubic_spline.h"
#include "math/bcl_math_cubic_spline.h"
#include "math/bcl_math_cubic_spline_damped.h"
#include "math/bcl_math_histogram.h"
#include "math/bcl_math_histogram_2d.h"
#include "math/bcl_math_limits.h"
#include "math/bcl_math_running_average.h"
#include "util/bcl_util_logger_interface.h"
#include "util/bcl_util_si_ptr_vector.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    //! @brief GetEnumDescriptor provides the name of ENUM
    //! @param ENUM - the enum for which a name is desired
    const std::string &GnuplotHeatmap::GetPaletteTypeString( const PaletteTypes &ENUM)
    {
      static const std::string s_descriptors[] =
      {
        "Default",
        "GreeRedViolet",
        "GreyScale",
        "BlackBlueVioletYellowWhite",
        "BlueGreenYellowRed",
        "WhiteBlueGreenRed",
        "BlueGreenRedWhite",
        "None",
        GetStaticClassName< PaletteTypes>()
      };

      return s_descriptors[ ENUM];
    }

    //! @brief Gnuplot Palette String provides the string to set the palette
    //! @param ENUM - the enum for which the palette is desired
    const std::string &GnuplotHeatmap::GetPaletteString( const PaletteTypes &ENUM)
    {
      static const std::string s_palette_strings[] =
      {
        "rgbformulae 22, 13, -31",
        "rgbformulae 3, 11, 6 # green-red-violet",
        "grey negative",
        "rgbformulae 30,31,32 # color printable on gray (black-blue-violet-yellow-white)",
        "rgbformulae 33,13,10 # rainbow (blue-green-yellow-red)",
        "defined (0 1 1 1, 0.00001 0 0 1, 1 0 1 0, 2 1 0 0)",
        "defined (0 0 0 1, 1 0 1 0, 2 1 0 0, 2.00001 1 1 1)",
        "#None",
        "#no palette set"
      };

      return s_palette_strings[ ENUM];
    }

  //////////
  // data //
  //////////

    //! symbol that represents Angstrom
    const std::string GnuplotHeatmap::s_AngstromSymbolGnuplot = "\\305";

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    GnuplotHeatmap::GnuplotHeatmap() :
      m_Data( 0, 0, 0.0),
      m_MinZ( util::GetUndefined< double>()),
      m_MaxZ( util::GetUndefined< double>()),
      m_Title(  ""),
      m_LabelX( ""),
      m_LabelY( ""),
      m_LabelZ( ""),
      m_TicsX(    ),
      m_TicsY(    ),
      m_TicsXCenter( true),
      m_TicsYCenter( true),
      m_BinsPerTicX( 1),
      m_BinsPerTicY( 1),
      m_PixelX( 1080),
      m_PixelY(  800),
      m_Ratio( util::GetUndefined< double>()),
      m_RotateXTics( 0.0),
      m_Filename( "heatmap"),
      m_Font( "Arial"),
      m_FontSize( 12),
      m_BoxCoords(),
      m_BoxCoordsSpecifications(),
      m_ShowColorBox( true),
      m_WritePreHeader( true),
      m_WriteHeader( true),
      m_WriteBoxes( true),
      m_WriteData( true),
      m_Palette( e_Default),
      m_ShowXTics( true),
      m_NoMirrorTics( false),
      m_TopMargin   ( util::GetUndefinedDouble()),
      m_BottomMargin( util::GetUndefinedDouble()),
      m_RightMargin ( util::GetUndefinedDouble()),
      m_LeftMargin  ( util::GetUndefinedDouble())
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief set from a matrix
    //! @param MATRIC the matrix of values to use
    void GnuplotHeatmap::SetFromMatrix( const linal::Matrix< double> &MATRIX)
    {
      m_Data = MATRIX;
    }

    //! @brief set heatmap from histogram 1D
    //! @param HISTOGRAM1D 1d histogram
    //! @param VERTICAL true - vertical, false horizontal
    //! @param CENTER_TICS center the tics (false, between two bins)
    void GnuplotHeatmap::SetFromHistogram( const Histogram &HISTOGRAM, const bool VERTICAL, const bool CENTER_TICS)
    {
      SetFromHistograms( util::SiPtrVector< const Histogram>::Create( HISTOGRAM, HISTOGRAM), VERTICAL, CENTER_TICS);
    }

    //! @brief set heatmap from histograms
    //! @param HISTOGRAMS siptr vector of histograms
    //! @param VERTICAL true - vertical, false horizontal
    //! @param CENTER_TICS center the tics (false, between two bins)
    void GnuplotHeatmap::SetFromHistograms( const util::SiPtrVector< const Histogram> &HISTOGRAMS, const bool VERTICAL, const bool CENTER_TICS)
    {
      if( HISTOGRAMS.IsEmpty())
      {
        return;
      }

      const Histogram &template_hist( *HISTOGRAMS.FirstElement());
      const linal::Vector< double> binning
      (
        linal::FillVector< double>
        (
          template_hist.GetNumberOfBins() + size_t( !CENTER_TICS),
          template_hist.GetBoundaries().First() + ( CENTER_TICS ? 0.5 * template_hist.GetBinSize() : 0.0),
          template_hist.GetBinSize()
        )
      );
      const storage::Vector< std::string> tics( Gnuplot::TicsFromBinning( binning, 1, util::Format().W( 4)));

      if( VERTICAL)
      {
        m_Data = linal::Matrix< double>( template_hist.GetNumberOfBins(), HISTOGRAMS.GetSize(), 0.0);
        SetTicsY( tics, CENTER_TICS, 1);
      }
      else
      {
        m_Data = linal::Matrix< double>( HISTOGRAMS.GetSize(), template_hist.GetNumberOfBins(), 0.0);
        SetTicsX( tics, CENTER_TICS, 1);
      }

      size_t count( 0);
      // iterate over all histograms
      for
      (
        util::SiPtrVector< const Histogram>::const_iterator itr( HISTOGRAMS.Begin()), itr_end( HISTOGRAMS.End());
        itr != itr_end;
        ++itr, ++count
      )
      {
        VERTICAL ? m_Data.ReplaceCol( count, ( *itr)->GetHistogram()) : m_Data.ReplaceRow( count, ( *itr)->GetHistogram());
      }
    }

    //! @brief set heatmap from unary function
    //! @param FUNCTION the function to use
    //! @param NUMBER_FUNCTION_VALUES number of function values to be calculated
    //! @param START first argument
    //! @param DELTA difference between two consecutive arguments
    //! @param VERTICAL true - vertical, false - horizontal
    //! @param CENTER_TICS center the tics (false, between two bins)
    //! @param NORMALIZE if true, sum of the function values will be normalized to one
    //! @param TICS_PER_DELTA number of delta-units per x-axis tic
    void GnuplotHeatmap::SetFromFunction
    (
      const FunctionInterfaceSerializable< double, double> &FUNCTION,
      const size_t NUMBER_FUNCTION_VALUES,
      const double START,
      const double DELTA,
      const bool VERTICAL,
      const bool CENTER_TICS,
      const bool NORMALIZE,
      const size_t TICS_PER_DELTA
    )
    {
      SetFromFunctions
      (
        util::SiPtrVector< const FunctionInterfaceSerializable< double, double> >::Create( FUNCTION, FUNCTION),
        NUMBER_FUNCTION_VALUES, START, DELTA, VERTICAL, CENTER_TICS, NORMALIZE, TICS_PER_DELTA
      );
    }

    //! @brief set heatmap from cubic spline
    //! @param CUBIC_SPLINE the spline to use
    //! @param VERTICAL true - vertical, false horizontal
    //! @param CENTER_TICS center the tics (false, between two bins)
    void GnuplotHeatmap::SetFromCubicSpline
    (
      const CubicSpline &CUBIC_SPLINE,
      const bool VERTICAL,
      const bool CENTER_TICS
    )
    {
      SetFromFunction
      (
        CUBIC_SPLINE,
        CUBIC_SPLINE.GetValues().GetSize(),
        CUBIC_SPLINE.GetStart(),
        CUBIC_SPLINE.GetDelta(),
        VERTICAL,
        CENTER_TICS,
        false
      );
    }

    //! @brief set heatmap from cubic spline
    //! @param CUBIC_SPLINE the spline to use
    //! @param VERTICAL true - vertical, false horizontal
    //! @param CENTER_TICS center the tics (false, between two bins)
    void GnuplotHeatmap::SetFromCubicSpline
    (
      const CubicSplineDamped &CUBIC_SPLINE,
      const bool VERTICAL,
      const bool CENTER_TICS
    )
    {
      SetFromFunction
      (
        CUBIC_SPLINE,
        CUBIC_SPLINE.GetXValues().GetSize() * 4,
        CUBIC_SPLINE.GetXValues().First(),
        CUBIC_SPLINE.GetDelta() / 4.0,
        VERTICAL,
        CENTER_TICS,
        false,
        4
      );
    }

    //! @brief set heatmap from unary functions
    //! @param FUNCTIONS the functions to use
    //! @param NUMBER_FUNCTION_ARGUMENTS number of function values to be calculated
    //! @param START first argument
    //! @param DELTA difference between two consecutive arguments
    //! @param VERTICAL true - vertical, false horizontal
    //! @param CENTER_TICS center the tics (false, between two bins)
    //! @param NORMALIZE_EACH_FUNCTION if true, sum of the values for each function will be normalized to one
    //! @param TICS_PER_DELTA number of delta-units per x-axis tic
    void GnuplotHeatmap::SetFromFunctions
    (
      const util::SiPtrVector< const FunctionInterfaceSerializable< double, double> > &FUNCTIONS,
      const size_t NUMBER_FUNCTION_VALUES,
      const double START,
      const double DELTA,
      const bool VERTICAL,
      const bool CENTER_TICS,
      const bool NORMALIZE_EACH_FUNCTION,
      const size_t TICS_PER_DELTA
    )
    {
      if( FUNCTIONS.IsEmpty())
      {
        return;
      }

      const linal::Vector< double> binning
      (
        linal::FillVector< double>
        (
          NUMBER_FUNCTION_VALUES + size_t( !CENTER_TICS), START - ( CENTER_TICS ? 0.0 : 0.5 * DELTA), DELTA
        )
      );
      const storage::Vector< std::string> tics( TicsFromBinning( binning, TICS_PER_DELTA, util::Format()));

      if( VERTICAL)
      {
        m_Data = linal::Matrix< double>( NUMBER_FUNCTION_VALUES, FUNCTIONS.GetSize(), 0.0);
        SetTicsY( tics, CENTER_TICS, TICS_PER_DELTA);
      }
      else
      {
        m_Data = linal::Matrix< double>( FUNCTIONS.GetSize(), NUMBER_FUNCTION_VALUES, 0.0);
        SetTicsX( tics, CENTER_TICS, TICS_PER_DELTA);
      }

      size_t count( 0);
      // iterate over all histograms
      for
      (
        util::SiPtrVector< const FunctionInterfaceSerializable< double, double> >::const_iterator
          itr( FUNCTIONS.Begin()), itr_end( FUNCTIONS.End());
        itr != itr_end;
        ++itr, ++count
      )
      {
        const FunctionInterfaceSerializable< double, double> &function( **itr);
        // create the row from function values in the spline
        linal::Vector< double> val( NUMBER_FUNCTION_VALUES);
        double *val_ptr( val.Begin());
        for( size_t bin( 0); bin < NUMBER_FUNCTION_VALUES; ++bin, ++val_ptr)
        {
          *val_ptr = function( START + bin * DELTA);
        }
        if( NORMALIZE_EACH_FUNCTION)
        {
          val.SetToSum( 1.0);
        }
        VERTICAL ? m_Data.ReplaceCol( count, val) : m_Data.ReplaceRow( count, val);
      }
    }

    //! @brief set heatmap from cubic splines
    //! @param CUBIC_SPLINES the splines to use
    //! @param VERTICAL true - vertical, false horizontal
    //! @param CENTER_TICS center the tics (false, between two bins)
    void GnuplotHeatmap::SetFromCubicSplines
    (
      const util::SiPtrVector< const CubicSpline> &CUBIC_SPLINES,
      const bool VERTICAL,
      const bool CENTER_TICS
    )
    {
      // check if splines are defined
      if( CUBIC_SPLINES.IsEmpty() && !CUBIC_SPLINES.IsDefined())
      {
        return;
      }

      const CubicSpline &template_spline( *CUBIC_SPLINES.FirstElement());
      SetFromFunctions
      (
        util::SiPtrVector< const FunctionInterfaceSerializable< double, double> >( CUBIC_SPLINES),
        template_spline.GetValues().GetSize() * 4,
        template_spline.GetStart(),
        template_spline.GetDelta() / 4.0,
        VERTICAL,
        CENTER_TICS,
        false,
        4
      );
    }

    //! @brief set heatmap from cubic splines
    //! @param CUBIC_SPLINES the splines to use
    //! @param VERTICAL true - vertical, false horizontal
    //! @param CENTER_TICS center the tics (false, between two bins)
    void GnuplotHeatmap::SetFromCubicSplines
    (
      const util::SiPtrVector< const CubicSplineDamped> &CUBIC_SPLINES,
      const bool VERTICAL,
      const bool CENTER_TICS
    )
    {
      // check if splines are defined
      if( CUBIC_SPLINES.IsEmpty() && !CUBIC_SPLINES.IsDefined())
      {
        return;
      }

      RunningAverage< double> min_val;
      RunningAverage< double> max_size;
      for( auto itr( CUBIC_SPLINES.Begin()), itr_end( CUBIC_SPLINES.End()); itr != itr_end; ++itr)
      {
        max_size +=( *itr)->GetXValues().GetSize();
        min_val += ( *itr)->GetXValues().First();
      }
      const CubicSplineDamped &template_spline( *CUBIC_SPLINES.FirstElement());
      SetFromFunctions
      (
        util::SiPtrVector< const FunctionInterfaceSerializable< double, double> >( CUBIC_SPLINES),
        size_t( max_size.GetAverage()) * 4,
        min_val.GetAverage(),
        template_spline.GetDelta() / 4,
        VERTICAL,
        CENTER_TICS,
        false,
        4
      );
    }

    //! @brief set heatmap from histogram 2D
    //! @param HISTOGRAM2D 2d histogram
    //! @param CENTER_TICS_X center the tics x axis (false, between two bins)
    //! @param CENTER_TICS_Y center the tics y axis (false, between two bins)
    void GnuplotHeatmap::SetFromHistogram
    (
      const Histogram2D &HISTOGRAM2D,
      const bool CENTER_TICS_X,
      const bool CENTER_TICS_Y
    )
    {
      m_Data = HISTOGRAM2D.GetHistogram();

      // binning tics
      const linal::Vector< double> binning_x
      (
        linal::FillVector< double>
        (
          HISTOGRAM2D.GetNumberOfBinsX() + size_t( !CENTER_TICS_X),
          HISTOGRAM2D.GetBoundariesX().First() + ( CENTER_TICS_X ? 0.5 * HISTOGRAM2D.GetBinSizeXY().First() : 0.0),
          HISTOGRAM2D.GetBinSizeXY().First()
        )
      );
      const linal::Vector< double> binning_y
      (
        linal::FillVector< double>
        (
          HISTOGRAM2D.GetNumberOfBinsY() + size_t( !CENTER_TICS_Y),
          HISTOGRAM2D.GetBoundariesY().First() + ( CENTER_TICS_Y ? 0.5 * HISTOGRAM2D.GetBinSizeXY().Second() : 0.0),
          HISTOGRAM2D.GetBinSizeXY().Second()
        )
      );
      const storage::Vector< std::string> xtics( TicsFromBinning( binning_x, 1, util::Format().W( 4)));
      const storage::Vector< std::string> ytics( TicsFromBinning( binning_y, 1, util::Format().W( 4)));
      SetTicsX( xtics, CENTER_TICS_X, 1);
      SetTicsY( ytics, CENTER_TICS_Y, 1);
    }

    //! @brief set heatmap from binary function
    //! @param BINARY_FUNCTION the function to use
    //! @param NUMBER_FUNCTION_VALUES_X number of function values to be calculated in x
    //! @param START_X first argument x
    //! @param DELTA_X difference between two consecutive arguments x
    //! @param NUMBER_FUNCTION_VALUES_Y number of function values to be calculated in y
    //! @param START_Y first argument y
    //! @param DELTA_Y difference between two consecutive arguments y
    //! @param CENTER_TICS_X center the tics x axis (false, between two bins)
    //! @param CENTER_TICS_Y center the tics y axis (false, between two bins)
    void GnuplotHeatmap::SetFromBinaryFunction
    (
      const BinaryFunctionInterface< double, double, double> &BINARY_FUNCTION,
      const size_t NUMBER_FUNCTION_VALUES_X,
      const double START_X,
      const double DELTA_X,
      const size_t NUMBER_FUNCTION_VALUES_Y,
      const double START_Y,
      const double DELTA_Y,
      const bool CENTER_TICS_X,
      const bool CENTER_TICS_Y
    )
    {
      // data
      m_Data = linal::Matrix< double>( NUMBER_FUNCTION_VALUES_Y, NUMBER_FUNCTION_VALUES_X, double( 0.0));

      // binning
      const linal::Vector< double> binning_x
      (
        linal::FillVector< double>
        (
          NUMBER_FUNCTION_VALUES_X, START_X - ( CENTER_TICS_X ? 0.0 : 0.5 * DELTA_X), DELTA_X
        )
      );
      const linal::Vector< double> binning_y
      (
        linal::FillVector< double>
        (
          NUMBER_FUNCTION_VALUES_Y, START_Y - ( CENTER_TICS_Y ? 0.0 : 0.5 * DELTA_Y), DELTA_Y
        )
      );

      for( size_t count_x( 0); count_x < NUMBER_FUNCTION_VALUES_Y; ++count_x)
      {
        for( size_t count_y( 0); count_y < NUMBER_FUNCTION_VALUES_X; ++count_y)
        {
          const double val_Y( START_Y + count_x * DELTA_Y);
          const double val_X( START_X + count_y * DELTA_X);

          m_Data( count_x, count_y) = BINARY_FUNCTION( val_X, val_Y);
        }
      }

      SetTicsX( TicsFromBinning( binning_x, 1, util::Format().FFP( 1)), CENTER_TICS_X, 1);
      SetTicsY( TicsFromBinning( binning_y, 1, util::Format().FFP( 2)), CENTER_TICS_Y, 1);
    }

    //! @brief set heatmap from bicubic spline
    //! @param BICUBIC_SPLINE the bicubic spline to use
    //! @param CENTER_TICS_X center the tics x axis (false, between two bins)
    //! @param CENTER_TICS_Y center the tics y axis (false, between two bins)
    void GnuplotHeatmap::SetFromBicubicSpline
    (
      const BicubicSpline &BICUBIC_SPLINE,
      const bool CENTER_TICS_X,
      const bool CENTER_TICS_Y
    )
    {
      SetFromBinaryFunction
      (
        BICUBIC_SPLINE,
        BICUBIC_SPLINE.GetValues().GetNumberRows() * 4.0, BICUBIC_SPLINE.GetStartX(), BICUBIC_SPLINE.GetDeltaX() / 4.0,
        BICUBIC_SPLINE.GetValues().GetNumberCols() * 4.0, BICUBIC_SPLINE.GetStartY(), BICUBIC_SPLINE.GetDeltaY() / 4.0,
        CENTER_TICS_X,
        CENTER_TICS_Y
      );
    }

    //! @brief set title and labels
    //! @param TITLE title for the map
    //! @param LABEL_X label for the x axis
    //! @param LABEL_Y label for the y axis
    //! @param LABEL_Z label for the z axis
    void GnuplotHeatmap::SetTitleAndLabel
    (
      const std::string &TITLE,
      const std::string &LABEL_X,
      const std::string &LABEL_Y,
      const std::string &LABEL_Z
    )
    {
      m_Title = TITLE;
      m_LabelX = LABEL_X;
      m_LabelY = LABEL_Y;
      m_LabelZ = LABEL_Z;
    }

    //! @brief set the min and max of z
    //! @param MIN_Z minimum on z - if undefined, pick the lowest automatically
    //! @param MAX_Z maximum on z - if undefined, pick the highest automatically
    void GnuplotHeatmap::SetMinMaxZ( const double MIN_Z, const double MAX_Z)
    {
      m_MinZ = MIN_Z;
      m_MaxZ = MAX_Z;
    }

    //! @brief set tics for x axis - the two outside tics are placed at either the outer boundary or on the two outer bins
    //! @param TICS_X vector of strings
    //! @param CENTER_TICS center the tics (false, between two bins)
    //! @param NTH_BIN every nth bin a tick is added - same number of bins to left and right are required
    //! @return true if the number of tics required matches the number of bins, so that tics are symmetric
    bool GnuplotHeatmap::SetTicsX( const storage::Vector< std::string> &TICS_X, const bool CENTER_TICS, const size_t NTH_BIN)
    {
      // empty data or tics
      if( m_Data.IsEmpty() || TICS_X.IsEmpty() || NTH_BIN == 0)
      {
        BCL_MessageCrt( "m_Data.IsEmpty() || TICS_X.IsEmpty() || NTH_BIN == 0");
        return false;
      }

      if( CENTER_TICS)
      {
        if( m_Data.GetNumberCols() - 1 != ( TICS_X.GetSize() - 1) * NTH_BIN)
        {
          BCL_MessageCrt
          (
            "m_Data.GetNumberCols() - 1 != ( TICS_X.GetSize() - 1) * NTH_BIN but center tics"
          );
          BCL_MessageCrt
          (
            "m_Data.GetNumberCols() - 1 " + util::Format()( m_Data.GetNumberCols())
          );
          BCL_MessageCrt
          (
            "( TICS_X.GetSize() - 1) " + util::Format()( ( TICS_X.GetSize() - 1))
          );
          BCL_MessageCrt
          (
            "NTH_BIN " + util::Format()( NTH_BIN)
          );
          return false;
        }
      }
      else
      {
        if( m_Data.GetNumberCols() != ( TICS_X.GetSize() - 1) * NTH_BIN)
        {
          BCL_MessageCrt
          (
            "m_Data.GetNumberCols() != ( TICS_X.GetSize() - 1) * NTH_BIN but not center tics"
          );
          BCL_MessageCrt
          (
            "m_Data.GetNumberCols() " + util::Format()( m_Data.GetNumberCols())
          );
          BCL_MessageCrt
          (
            "( TICS_X.GetSize() - 1) " + util::Format()( ( TICS_X.GetSize() - 1))
          );
          BCL_MessageCrt
          (
            "NTH_BIN " + util::Format()( NTH_BIN)
          );
          return false;
        }
      }

      m_TicsX       = TICS_X;
      m_TicsXCenter = CENTER_TICS;
      m_BinsPerTicX = NTH_BIN;

      return true;
    }

    //! @brief set tics for y axis - the two outside tics are placed at either the outer boundary or on the two outer bins
    //! @param TICS_Y vector of strings
    //! @param CENTER_TICS center the tics (false, to the top and bottom of the bin)
    //! @param NTH_BIN every nth bin a tick is added - same number of bins to left and right are required
    //! @return true if the number of tics required matches the number of bins, so that tics are symmetric
    bool GnuplotHeatmap::SetTicsY( const storage::Vector< std::string> &TICS_Y, const bool CENTER_TICS, const size_t NTH_BIN)
    {
      // empty data or tics
      if( m_Data.IsEmpty() || TICS_Y.IsEmpty() || NTH_BIN == 0)
      {
        return false;
      }

      if( CENTER_TICS)
      {
        if( m_Data.GetNumberRows() - 1 != ( TICS_Y.GetSize() - 1) * NTH_BIN)
        {
          return false;
        }
      }
      else
      {
        if( m_Data.GetNumberRows() != ( TICS_Y.GetSize() - 1) * NTH_BIN)
        {
          return false;
        }
      }

      m_TicsY       = TICS_Y;
      m_TicsYCenter = CENTER_TICS;
      m_BinsPerTicY = NTH_BIN;

      return true;
    }

    //! @brief set the pixel size and ratio
    //! @param PIXEL_X number of pixel of x
    //! @param PIXEL_Y number of pixel of y
    //! @param RATIO ratio for plot y / x - if undefined, no roatio will be set
    void GnuplotHeatmap::SetPixelAndRatio( const size_t PIXEL_X, const size_t PIXEL_Y, const double RATIO)
    {
      m_PixelX = PIXEL_X;
      m_PixelY = PIXEL_Y;
      m_Ratio  = RATIO;
    }

    //! @brief set the font
    //! @param FONT string of font name
    //! @param FONT_SIZE size of font in pixel
    void GnuplotHeatmap::SetFont( const std::string &FONT, const size_t FONT_SIZE)
    {
      m_Font = FONT;
      m_FontSize = FONT_SIZE;
    }

    //! @brief sets the coordinates for boxes
    //! @param BOX_COORDS the lower left, upper right coordinates for each box [(x,y),(x,y)]
    //! @param X_RANGE the range in the x direction of the plot - the total size of the x direction
    //! @param Y_RANGE the range in the y direction of the plot -  the total size of the y direction
    //! @param X_MIN the minimum value in the x direction of the heat map
    //! @param Y_MIN the minimum value in the y direction of the heat map
    void GnuplotHeatmap::SetBoxes
    (
      const storage::Vector< storage::VectorND< 2, storage::VectorND< 2, double> > > &BOX_COORDS,
      const double X_RANGE, const double Y_RANGE, const double X_MIN, const double Y_MIN
    )
    {
      BCL_Assert( !m_Data.IsEmpty(), "set the gnuplot before setting boxes");

      storage::Vector< storage::VectorND< 2, storage::VectorND< 2, double> > > scaled_coords;

      const double x_range( m_Data.GetNumberCols());
      const double y_range( m_Data.GetNumberRows());
      static const double s_x_start( -0.5);
      static const double s_y_start( -0.5);

      // iterate over the boxes to scale the coordinates into the graph
      for
      (
        storage::Vector< storage::VectorND< 2, storage::VectorND< 2, double> > >::const_iterator
          box_itr( BOX_COORDS.Begin()), box_itr_end( BOX_COORDS.End());
        box_itr != box_itr_end;
        ++box_itr
      )
      {
        const double lower_left_x( box_itr->First().First());
        const double lower_left_y( box_itr->First().Second());
        const double scaled_lower_left_x( ( lower_left_x - X_MIN) / X_RANGE * x_range + s_x_start);
        const double scaled_lower_left_y( ( lower_left_y - Y_MIN) / Y_RANGE * y_range + s_y_start);
        const storage::VectorND< 2, double> scaled_lower_left( scaled_lower_left_x, scaled_lower_left_y);

        const double upper_right_x( box_itr->Second().First());
        const double upper_right_y( box_itr->Second().Second());
        const double scaled_upper_right_x( ( upper_right_x - X_MIN) / X_RANGE * x_range + s_x_start);
        const double scaled_upper_right_y( ( upper_right_y - Y_MIN) / Y_RANGE * y_range + s_y_start);
        BCL_MessageStd( "upper_right_y " + util::Format()( upper_right_y));
        const storage::VectorND< 2, double> scaled_upper_right( scaled_upper_right_x, scaled_upper_right_y);

        scaled_coords.PushBack( storage::VectorND< 2, storage::VectorND< 2, double> >( scaled_lower_left, scaled_upper_right));
      }

      m_BoxCoords = scaled_coords;
    }

    //! @brief sets the specifications for boxes
    //! @param SPECIFICATIONS strings which specify the format of the boxes
    void GnuplotHeatmap::SetBoxesSpecifications( const storage::Vector< std::string> &SPECIFICATIONS)
    {
      m_BoxCoordsSpecifications = SPECIFICATIONS;
    }

    //! @brief if set, then when WriteScript is called, the WritePreHeader string will be written. True by default
    void GnuplotHeatmap::SetWritePreHeader( const bool SET)
    {
      m_WritePreHeader = SET;
    }

    //! @brief if set, then when WriteScript is called, the WriteHeader    string will be written. True by default
    void GnuplotHeatmap::SetWriteHeader( const bool SET)
    {
      m_WriteHeader = SET;
    }

    //! @brief if set, then when WriteScript is called, the WriteBoxes     string will be written. True by default
    void GnuplotHeatmap::SetWriteBoxes( const bool SET)
    {
      m_WriteBoxes = SET;
    }

    //! @brief if set, then when WriteScript is called, the WriteData      string will be written. True by default
    void GnuplotHeatmap::SetWriteData( const bool SET)
    {
      m_WriteData = SET;
    }

    //! @brief set, the palette to use
    //! @param PALETTE the color palette used by the splot
    void GnuplotHeatmap::SetPalette( const PaletteTypes PALETTE)
    {
      m_Palette = PALETTE;
    }

    //! @brief if set, then the x tics will be shown
    void GnuplotHeatmap::SetShowXTics( const bool SET)
    {
      m_ShowXTics = SET;
    }

    //! @brief if set, then the the tics won't be mirrored on the opposite side of the axis
    void GnuplotHeatmap::SetNoMirrorTics( const bool SET)
    {
      m_NoMirrorTics = SET;
    }

    //! @brief set the ratio to which the actual plot area extends to the edge of the png
    //! @param TOP    ratio to which this plot extends to the top    edge of the png
    //! @param BOTTOM ratio to which this plot extends to the bottom edge of the png
    //! @param RIGHT  ratio to which this plot extends to the right  edge of the png
    //! @param LEFT   ratio to which this plot extends to the left   edge of the png
    void GnuplotHeatmap::SetMargins( const double TOP, const double BOTTOM, const double RIGHT, const double LEFT)
    {
      m_TopMargin    = TOP;
      m_BottomMargin = BOTTOM;
      m_RightMargin  = RIGHT;
      m_LeftMargin   = LEFT;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief multiply the heatmap
    //! @details can be used to get similar scaling between different heat maps
    //! @param FACTOR
    void GnuplotHeatmap::Muliply( const double FACTOR)
    {
      m_Data *= FACTOR;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write unchangeable information necessary for the plot
    //! @param OSTREAM stream to be written to
    //! @param TITLE title for the heat map
    std::ostream &GnuplotHeatmap::WritePreHeader( std::ostream &OSTREAM) const
    {
      // write comments
      OSTREAM << "# BCL generated heatmap\n";
      OSTREAM << "set terminal png enhanced transparent ";
      if( !m_Font.empty())
      {
        OSTREAM << "font \"" << m_Font << ',' << m_FontSize << "\" ";
      }
      OSTREAM << "size " << m_PixelX << ',' << m_PixelY << '\n';
      OSTREAM << "set output \"" << m_Filename << ".png\"\n";
      OSTREAM << "set encoding iso_8859_1\n";

      return OSTREAM;
    }

    //! @brief write a choice of palettes
    //! @param OSTREAM stream to be written to
    std::ostream &GnuplotHeatmap::WritePalettes( std::ostream &OSTREAM) const
    {
      OSTREAM << "set palette " << GetPaletteString( m_Palette) << '\n';

      // end
      return OSTREAM;
    }

    //! @brief write a header
    //! @param OSTREAM stream to be written to
    //! @param TITLE title for the heat map
    std::ostream &GnuplotHeatmap::WriteHeader( std::ostream &OSTREAM) const
    {
      OSTREAM << "set view map\n";
      if( !m_Title.empty())
      {
        OSTREAM << "set title \"" << m_Title << "\"\n";
      }
      else
      {
        OSTREAM << "unset title\n";
      }
      OSTREAM << "unset key\n\n";

      if( util::IsDefined( m_Ratio))
      {
        OSTREAM << "set size ratio " << m_Ratio << '\n';
      }

      // axes
      WriteAxisInfo( OSTREAM, coord::GetAxes().e_X, m_LabelX,   -0.5, m_Data.GetNumberCols() - 0.5, m_TicsX, m_TicsXCenter, m_BinsPerTicX);
      WriteAxisInfo( OSTREAM, coord::GetAxes().e_Y, m_LabelY,   -0.5, m_Data.GetNumberRows() - 0.5, m_TicsY, m_TicsYCenter, m_BinsPerTicY);
      WriteAxisInfo( OSTREAM, coord::GetAxes().e_Z, m_LabelZ, m_MinZ,                       m_MaxZ, storage::Vector< std::string>(), false, 0);

      // true if the color box should not be shown
      if( !m_ShowColorBox)
      {
        OSTREAM << "unset colorbox\n";
      }

      if( !m_ShowXTics)
      {
        OSTREAM << "set noxtics\n";
      }
      if( !m_NoMirrorTics)
      {
        OSTREAM << "\n";
      }

      // true if the user set the margins
      if
      (
        util::IsDefined( m_TopMargin)   && util::IsDefined( m_BottomMargin) &&
        util::IsDefined( m_RightMargin) && util::IsDefined( m_LeftMargin)
      )
      {
        // set the margins
        OSTREAM << "set lmargin screen " << m_LeftMargin   << "\n";
        OSTREAM << "set rmargin screen " << m_RightMargin  << "\n";
        OSTREAM << "set bmargin screen " << m_BottomMargin << "\n";
        OSTREAM << "set tmargin screen " << m_TopMargin    << "\n";
      }

//      OSTREAM << "set colorbox horizontal user origin 0.1,0.9 size 0.8,0.01\n";

      // end
      return OSTREAM;
    }

    //! @brief write the whole information for an axis
    //! @param OSTREAM the stream to write to
    //! @param AXIS the axis that is written
    //! @param LABEL the axis label - if empty, it will be disabled
    //! @param MIN the minimal value for this axis
    //! @param MAX the max value on this axis
    //! @param TICS tics for the axis
    //! @param CENTER_TICS center the tics (false, to the top and bottom of the bin)
    //! @param NTH_BIN every nth bin a tick is added - same number of bins to left and right are required
    std::ostream &GnuplotHeatmap::WriteAxisInfo
    (
      std::ostream &OSTREAM,
      const coord::Axis &AXIS,
      const std::string &LABEL,
      const double MIN,
      const double MAX,
      const storage::Vector< std::string> &TICS,
      const bool CENTER_TICS,
      const size_t NTH_BIN
    ) const
    {
      // string for that axis
      const std::string &axis_string( Gnuplot::GetAxisString( AXIS));

      // label
      if( !LABEL.empty())
      {
        OSTREAM << "set " << axis_string << "label \"" << LABEL << "\"";
        // Move axis label further away if tick labels are rotated
        OSTREAM << ( ( AXIS == coord::GetAxes().e_X && m_RotateXTics != 0) ? " offset 0,-1\n" : "\n");
      }

      // range
      OSTREAM << "set " << axis_string << "range [ " << ( util::IsDefined( MIN) ? util::Format()( MIN) : "*") << " : " << ( util::IsDefined( MAX) ? util::Format()( MAX) : "*") << " ]\n";

      // tics
      if( TICS.IsEmpty())
      {
        if( AXIS != coord::GetAxes().e_Z)
        {
          OSTREAM << "set no" << axis_string << "tics\n";
        }
        else
        {
          OSTREAM << "#set cbtics 1\n";
          OSTREAM << "#set format cb \"%3.1f\"\n";
        }
      }
      else
      {
        OSTREAM << "set " << axis_string << "tics ";

        if( m_NoMirrorTics)
        {
          OSTREAM << "nomirror ";
        }

        if( AXIS == coord::GetAxes().e_X)
        {
          OSTREAM << "rotate by " << m_RotateXTics << ' ';
        }

        OSTREAM << "(";
        double current_bin( 0.0);
        if( !CENTER_TICS)
        {
          current_bin -= 0.5;
        }

        // prevent endless loop
        if( NTH_BIN == 0)
        {
          return OSTREAM;
        }

        for
        (
          storage::Vector< std::string>::const_iterator itr( TICS.Begin()), itr_end( TICS.End());
          itr != itr_end;
          current_bin += NTH_BIN
        )
        {
//          OSTREAM << "\"" << *itr << "     \" " << current_bin;
          OSTREAM << "\" " << *itr << " \" " << current_bin;
          ++itr;
          if( itr == itr_end)
          {
            break;
          }
          OSTREAM << ", ";
        }
        OSTREAM << ")\n";
      }

      // end
      return OSTREAM;
    }

    //! @brief write any boxes to the script for display
    //! @param OSTREAM the stream the boxes will be written to
    std::ostream &GnuplotHeatmap::WriteBoxes( std::ostream &OSTREAM) const
    {
      static const std::string s_box_color( "black");
      static const std::string s_fill_style( "empty");
      static const size_t s_border_color( 0);
      static const double s_border_linewidth( 0.9);
      static const std::string s_border_command
      (
        "border " + util::Format()( s_border_color) + " linewidth " + util::Format()( s_border_linewidth)
      );
      static const std::string s_default_specifications
      (
        " front fillcolor rgb \"" + s_box_color + "\" fs " + s_fill_style + " " + s_border_command
      );

      size_t counter( 0);
      // iterate through the boxes
      for
      (
        storage::Vector< storage::VectorND< 2, storage::VectorND< 2, double> > >::const_iterator
          box_itr( m_BoxCoords.Begin()), box_itr_end( m_BoxCoords.End());
        box_itr != box_itr_end;
        ++box_itr, ++counter
      )
      {
        const double lower_left_x( box_itr->First().First());
        const double lower_left_y( box_itr->First().Second());
        const double upper_right_x( box_itr->Second().First());
        const double upper_right_y( box_itr->Second().Second());

        std::string specifications( s_default_specifications);
        if( m_BoxCoordsSpecifications.GetSize() > counter)
        {
          specifications = m_BoxCoordsSpecifications( counter);
        }

        OSTREAM << "set object rect from " << lower_left_x << "," << lower_left_y << " to "
          << upper_right_x << "," << upper_right_y
          << " " << specifications << "\n";
      }

      return OSTREAM;
    }

    //! @brief write the data
    //! @param OSTREAM the stream to write to
    //! @return the stream written to
    std::ostream &GnuplotHeatmap::WriteData( std::ostream &OSTREAM) const
    {
      OSTREAM << "splot '-' using 1:2:3 with image\n";
      OSTREAM << "# number x values " << m_Data.GetNumberCols() << '\n';
      OSTREAM << "# number y values " << m_Data.GetNumberRows() << '\n';

      // ptr to first element
      const double *ptr( m_Data.Begin());

      // iterate over rows
      for( size_t row( 0); row != m_Data.GetNumberRows(); ++row)
      {
        for( size_t col( 0); col != m_Data.GetNumberCols(); ++col, ++ptr)
        {
          OSTREAM << col << '\t' << row << " \t" << *ptr << '\n';
        }
        OSTREAM << '\n';
      }

      // indicate end of data for heatmap
      OSTREAM << "e\n";

      // end
      return OSTREAM;
    }

    //! @brief write the gnuplot script
    //! @param OSTREAM the stream to write to
    void GnuplotHeatmap::WriteScript( std::ostream &OSTREAM) const
    {
      if( m_WritePreHeader)
      {
        WritePreHeader( OSTREAM);
      }
      if( m_WriteHeader)
      {
        WriteHeader(    OSTREAM);
      }
      if( m_Palette != e_None)
      {
        WritePalettes( OSTREAM);
      }
      if( m_WriteBoxes)
      {
        WriteBoxes(     OSTREAM);
      }
      if( m_WriteData)
      {
        WriteData(      OSTREAM);
      }
    }

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_gnuplot_multiplot.h"

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
    const util::SiPtr< const util::ObjectInterface> GnuplotMultiplot::s_Instance
    (
      GetObjectInstances().AddInstance( new GnuplotMultiplot())
    );
  
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    GnuplotMultiplot::GnuplotMultiplot() :
      m_Plots(),
      m_Rows( 1),
      m_Cols( 1),
      m_Filename( "heatmap"),
      m_Font( "Arial"),
      m_FontSize( 12),
      m_PixelX( 600),
      m_PixelY( 400),
      m_Ratio( util::GetUndefinedDouble()),
      m_WritePreHeader( true),
      m_WritePalettes( true)
    {
    }

    //! @brief Clone function
    //! @return pointer to new GnuplotMultiplot
    GnuplotMultiplot *GnuplotMultiplot::Clone() const
    {
      return new GnuplotMultiplot( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GnuplotMultiplot::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief set tics for x axis - the two outside tics are placed at either the outer boundary or on the two outer bins
    //! @param TICS_X vector of strings
    //! @param CENTER_TICS center the tics (false, between two bins)
    //! @param NTH_BIN every nth bin a tick is added - same number of bins to left and right are required
    //! @return true if the number of tics required matches the number of bins, so that tics are symmetric
    bool GnuplotMultiplot::SetTicsX
    (
      const storage::Vector< std::string> &TICS_X, const bool CENTER_TICS, const size_t NTH_BIN
    )
    {
      return true;
    }

    //! @brief set tics for y axis - the two outside tics are placed at either the outer boundary or on the two outer bins
    //! @param TICS_Y vector of strings
    //! @param CENTER_TICS center the tics (false, to the top and bottom of the bin)
    //! @param NTH_BIN every nth bin a tick is added - same number of bins to left and right are required
    //! @return true if the number of tics required matches the number of bins, so that tics are symmetric
    bool GnuplotMultiplot::SetTicsY
    (
      const storage::Vector< std::string> &TICS_Y, const bool CENTER_TICS, const size_t NTH_BIN
    )
    {
      return true;
    }

    //! @brief set the pixel size and ratio
    //! @param PIXEL_X number of pixel of x
    //! @param PIXEL_Y number of pixel of y
    //! @param RATIO ratio for plot y / x
    void GnuplotMultiplot::SetPixelAndRatio( const size_t PIXEL_X, const size_t PIXEL_Y, const double RATIO)
    {
      m_PixelX = PIXEL_X;
      m_PixelY = PIXEL_Y;
      m_Ratio  = RATIO;
    }

    //! @brief set the rotation of the xtics
    //! @param DEGREE_ROTATION degree of rotation - 0 non
    void GnuplotMultiplot::SetRotationXTics( const double DEGREE_ROTATION)
    {
      return;
    }

    //! @brief set the filename of the generated plot
    //! @param FILENAME filename
    void GnuplotMultiplot::SetFilename( const std::string &FILENAME)
    {
      m_Filename = FILENAME;
    }

    //! @brief set the font
    //! @param FONT string of font name
    //! @param FONT_SIZE size of font in pixel
    void GnuplotMultiplot::SetFont( const std::string &FONT, const size_t FONT_SIZE)
    {
      m_Font = FONT;
      m_FontSize = FONT_SIZE;
    }

    //! @brief if set, then when WriteScript is called, the WritePreHeader string will be written
    void GnuplotMultiplot::SetWritePreHeader( const bool SET)
    {
      m_WritePreHeader = SET;
    }

    //! @brief if set, then when WriteScript is called, the WriteHeader string will be written
    void GnuplotMultiplot::SetWriteHeader( const bool SET)
    {
      return;
    }

    //! @brief if set, then when WriteScript is called, the WriteBoxes string will be written
    void GnuplotMultiplot::SetWriteBoxes( const bool SET)
    {
      return;
    }

    //! @brief if set, then when WriteScript is called, the WriteData string will be written
    void GnuplotMultiplot::SetWriteData( const bool SET)
    {
      return;
    }

    //! @brief set the row and column layout of the multiplot
    void GnuplotMultiplot::SetRowsCols( size_t ROWS, size_t COLS)
    {
      m_Rows = ROWS;
      m_Cols = COLS;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief add a plot to the multiplot
    //! @param PLOT the plot which will be added to the multiplot
    void GnuplotMultiplot::Insert( const util::ShPtr< Gnuplot> &PLOT)
    {
      m_Plots.PushBack( PLOT);
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write unchangeable information necessary for the plot
    //! @param OSTREAM stream to be written to
    //! @param TITLE title for the heat map
    std::ostream &GnuplotMultiplot::WritePreHeader( std::ostream &OSTREAM) const
    {
      // write comments
      OSTREAM << "# BCL generated heatmap\n";
      OSTREAM << "set terminal png transparent ";
      if( !m_Font.empty())
      {
        OSTREAM << "font \"" << m_Font << "\" " << m_FontSize;
      }
      OSTREAM << " size " << m_PixelX << ',' << m_PixelY << '\n';
      if( util::IsDefined( m_Ratio))
      {
        OSTREAM << "set size ratio " << m_Ratio << '\n';
      }
      OSTREAM << "set output \"" << m_Filename << ".png\"\n";
      OSTREAM << "set encoding iso_8859_1\n";
      OSTREAM << "set multiplot layout " << m_Rows << "," << m_Cols << "\n";

      return OSTREAM;
    }

    //! @brief write a header
    //! @param OSTREAM stream to be written to
    std::ostream &GnuplotMultiplot::WriteHeader( std::ostream &OSTREAM) const
    {
      return OSTREAM;
    }

    //! @brief write the whole information for an axis
    //! @param OSTREAM the stream to write to
    //! @param AXIS the axis that is written
    //! @param LABEL the axis label - if empty, it will be disabled
    //! @param MIN the minimal value for this axis
    //! @param MAX the max value on this axis
    //! @param TICS tics for the axis
    //! @param CENTER_TICS center the tics (false, to the top and bottom of the bin)
    //! @param NTH_BIN every nth bin a tick is added - same number of bins to left and right are required
    std::ostream &GnuplotMultiplot::WriteAxisInfo
    (
      std::ostream &OSTREAM,
      const coord::Axis &AXIS,
      const std::string &LABEL,
      const double MIN,
      const double MAX,
      const storage::Vector< std::string> &TICS,
      const bool CENTER_TICS,
      const size_t NTH_BIN
    ) const
    {
      return OSTREAM;
    }

    //! @brief write the data
    //! @param OSTREAM the stream to write to
    //! @return the stream written to
    std::ostream &GnuplotMultiplot::WriteData( std::ostream &OSTREAM) const
    {
      return OSTREAM;
    }

    //! @brief write the gnuplot script
    //! @param OSTREAM the stream to write to
    void GnuplotMultiplot::WriteScript( std::ostream &OSTREAM) const
    {
      if( m_WritePreHeader)
      {
        WritePreHeader( OSTREAM);
      }

      // iterate through the plots to print them all out
      for
      (
        util::ShPtrList< Gnuplot>::const_iterator plot_itr( m_Plots.Begin()), plot_itr_end( m_Plots.End());
        plot_itr != plot_itr_end; ++plot_itr
      )
      {
        ( *plot_itr)->WriteScript( OSTREAM);
      }
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &GnuplotMultiplot::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &GnuplotMultiplot::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_histogram_2d.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_operations.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

  //////////
  // data //
  //////////

    //! @brief string to indicate the center
    //! @return string to indicate the center
    const std::string &Histogram2D::GetCenterString()
    {
      static std::string s_center_string( "center");
      return s_center_string;
    }

    //! @brief string to indicate counts
    //! @return string to indicate counts
    const std::string &Histogram2D::GetCountsString()
    {
      static std::string s_counts( "counts");
      return s_counts;
    }

    //! @brief string to indicate the left boundary in output
    //! @return string to indicate the left boundary in output
    const std::string &Histogram2D::GetLeftBoundaryString()
    {
      static std::string s_left_boundary_string( "...<");
      return s_left_boundary_string;
    }

    //! @brief string to indicate a bin
    //! @return string to indicate a bin
    const std::string &Histogram2D::GetBinString()
    {
      static std::string s_bin_string( "<..>");
      return s_bin_string;
    }

    //! @brief string to indicate the right boundary in output
    //! @return string to indicate the right boundary in output
    const std::string &Histogram2D::GetRightBoundaryString()
    {
      static std::string s_right_boundary_string( ">...");
      return s_right_boundary_string;
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Histogram2D::s_Instance
    (
      GetObjectInstances().AddInstance( new Histogram2D())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor - initializes everything to 0
    Histogram2D::Histogram2D() :
      m_LowerUpperBoundariesX( 0.0, 0.0),
      m_LowerUpperBoundariesY( 0.0, 0.0),
      m_BinSizeXY( 0.0, 0.0),
      m_Histogram( 0, 0)
    {
    }

    //! construct Histogram from X and Y starting values MIN_X_Y, from Pair BINSIZE_X_Y and from NUMBER_OF_BINS_X_Y
    Histogram2D::Histogram2D
    (
      const storage::VectorND< 2, double> &MIN_X_Y,
      const storage::VectorND< 2, double> &BINSIZE_X_Y,
      const storage::VectorND< 2, size_t> &NUMBER_OF_BINS_X_Y
    ) :
      m_LowerUpperBoundariesX( MIN_X_Y.First(), MIN_X_Y.First() + BINSIZE_X_Y.First() * NUMBER_OF_BINS_X_Y.First()),
      m_LowerUpperBoundariesY( MIN_X_Y.Second(), MIN_X_Y.Second() + BINSIZE_X_Y.Second() * NUMBER_OF_BINS_X_Y.Second()),
      m_BinSizeXY( BINSIZE_X_Y),
      m_Histogram( NUMBER_OF_BINS_X_Y.Second() + 2, NUMBER_OF_BINS_X_Y.First() + 2)
    {
    }

    //! copy constructor
    Histogram2D *Histogram2D::Clone() const
    {
      return new Histogram2D( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Histogram2D::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! get number of bins X
    size_t Histogram2D::GetNumberOfBinsX() const
    {
      return m_Histogram.GetNumberCols() - 2;
    }

    //! get number of bins Y
    size_t Histogram2D::GetNumberOfBinsY() const
    {
      return m_Histogram.GetNumberRows() - 2;
    }

    //! get binsize XY
    const storage::VectorND< 2, double> &Histogram2D::GetBinSizeXY() const
    {
      return m_BinSizeXY;
    }

    //! GetBoundaries X
    const storage::VectorND< 2, double> &Histogram2D::GetBoundariesX() const
    {
      return m_LowerUpperBoundariesX;
    }

    //! GetBoundaries Y
    const storage::VectorND< 2, double> &Histogram2D::GetBoundariesY() const
    {
      return m_LowerUpperBoundariesY;
    }

    //! GetBoundaries counts X
    storage::VectorND< 2, linal::Vector< double> > Histogram2D::GetBoundariesCountsX() const
    {
      return storage::VectorND< 2, linal::Vector< double> >( m_Histogram.GetCol( 0), m_Histogram.GetCol( m_Histogram.GetNumberCols() - 1));
    }

    //! GetBoundaries counts Y
    storage::VectorND< 2, linal::Vector< double> > Histogram2D::GetBoundariesCountsY() const
    {
      return storage::VectorND< 2, linal::Vector< double> >( m_Histogram.GetRow( 0), m_Histogram.GetRow( m_Histogram.GetNumberRows() - 1));
    }

    //! GetHistogram
    linal::Matrix< double> Histogram2D::GetHistogram() const
    {
      return m_Histogram.CreateSubMatrix( m_Histogram.GetNumberRows() - 2, m_Histogram.GetNumberCols() - 2, 1, 1);
    }

    //! get sum of all counts
    double Histogram2D::GetSumOfAllCounts() const
    {
      return m_Histogram.Sum();
    }

    //! GetBinning
    storage::VectorND< 2, linal::Vector< double> > Histogram2D::GetBinningXY() const
    {
      storage::VectorND< 2, linal::Vector< double> > binning( linal::Vector< double>( m_Histogram.GetNumberCols() - 2, double( 0.0)), linal::Vector< double>( m_Histogram.GetNumberRows() - 2, double( 0.0)));

      //binning in X direction
      double i( 0.5);
      for( double *ptr( binning.First().Begin()), *ptr_end( binning.First().End()); ptr != ptr_end; ++ptr, ++i)
      {
        *ptr = m_LowerUpperBoundariesX.First() + i * m_BinSizeXY.First();
      }

      //binning in Y direction
      i = 0.5;
      for( double *ptr( binning.Second().Begin()), *ptr_end( binning.Second().End()); ptr != ptr_end; ++ptr, ++i)
      {
        *ptr = m_LowerUpperBoundariesY.First() + i * m_BinSizeXY.Second();
      }

      //return
      return binning;
    }

  ////////////////
  // operations //
  ////////////////

    //! calculate the Histogram2D from a VALUES_VECTOR
    Histogram2D &Histogram2D::CalculateHistogram
    (
      const storage::Vector< storage::VectorND< 2, double> > &VALUES_VECTOR
    )
    {
      //reset all counts to zero
      Reset();

      //count
      for
      (
        std::vector< storage::VectorND< 2, double> >::const_iterator
          itr( VALUES_VECTOR.Begin()), itr_end( VALUES_VECTOR.End());
        itr != itr_end; ++itr
      )
      {
        PushBack( *itr);
      }

      return *this;
    }

    //! reset all counts
    void Histogram2D::Reset()
    {
      m_Histogram.SetZero();
    }

    //! pushback a pair of values to the right position in the histogram
    void Histogram2D::PushBack
    (
      const storage::VectorND< 2, double> &PAIR_OF_VALUES,
      const double &WEIGHT
    )
    {
      if( !util::IsDefined( PAIR_OF_VALUES.First()) || !util::IsDefined( PAIR_OF_VALUES.Second())) return;

      size_t bin_x( 0), bin_y( 0);

      //determine field in x-direction
      if( PAIR_OF_VALUES.First() < m_LowerUpperBoundariesX.First())
      {
        bin_x = 0;
      }
      else if( PAIR_OF_VALUES.First() >= m_LowerUpperBoundariesX.Second())
      {
        bin_x = m_Histogram.GetNumberCols() - 1;
      }
      else
      {
        bin_x = size_t( ( PAIR_OF_VALUES.First() - m_LowerUpperBoundariesX.First()) / m_BinSizeXY.First() + 1);
      }

      //determine field in y-direction
      if( PAIR_OF_VALUES.Second() < m_LowerUpperBoundariesY.First())
      {
        bin_y = 0;
      }
      else if( PAIR_OF_VALUES.Second() >= m_LowerUpperBoundariesY.Second())
      {
        bin_y = m_Histogram.GetNumberRows() - 1;
      }
      else
      {
        bin_y = size_t( ( PAIR_OF_VALUES.Second() - m_LowerUpperBoundariesY.First()) / m_BinSizeXY.Second() + 1);
      }

      m_Histogram( bin_y, bin_x) += WEIGHT;
    }

    //! @brief combine this with a given histogram2d by adding up all counts
    //! all parameters have to be identical for this operation to work
    //! @param HISTOGRAM_2D histogram to be added to this one
    //! @return true if it was successful
    bool Histogram2D::Combine( const Histogram2D &HISTOGRAM_2D)
    {
      // allow combining, if this is empty
      if( IsEmpty())
      {
        *this = HISTOGRAM_2D;
        return true;
      }

      // check that the parameters agree
      if
      (
           m_LowerUpperBoundariesX.First()  != HISTOGRAM_2D.m_LowerUpperBoundariesX.First()
        || m_LowerUpperBoundariesX.Second() != HISTOGRAM_2D.m_LowerUpperBoundariesX.Second()
        || m_LowerUpperBoundariesY.First()  != HISTOGRAM_2D.m_LowerUpperBoundariesY.First()
        || m_LowerUpperBoundariesY.Second() != HISTOGRAM_2D.m_LowerUpperBoundariesY.Second()
        || m_BinSizeXY.First() != HISTOGRAM_2D.m_BinSizeXY.First()
        || m_BinSizeXY.Second() != HISTOGRAM_2D.m_BinSizeXY.Second()
      )
      {
        return false;
      }

      // add the histogram
      m_Histogram += HISTOGRAM_2D.m_Histogram;

      // end
      return true;
    }

    //! checks if there is any count
    bool Histogram2D::IsEmpty() const
    {
      return GetSumOfAllCounts() == 0;
    }

    //! @brief normalize the counts in the histogram
    //! each bin is divided by the total number of counts (even the boundary counts)
    void Histogram2D::Normalize()
    {
      // create conts double "total_counts" and initialize with result of "GetSumOfAllCounts"
      const double total_counts( GetSumOfAllCounts());

      // if total_counts is 0 then return since no normalization needed
      if( !total_counts)
      {
        return;
      }

      // normalize
      m_Histogram /= total_counts;
    }

    //! @brief normalize the row counts in the histogram
    //! @detail each cell is divided by the maximum number of counts in the row
    void Histogram2D::NormalizeRows()
    {
      for( size_t row_index( 0), row_count( m_Histogram.GetNumberRows()); row_index < row_count; ++row_index)
      {
        // copy row from histogram and calculate the count
        linal::Vector< double> row( m_Histogram.GetRow( row_index));
        const double max_count( row.Max());

        // iterate over each element dividing by the total count of this column
        for( linal::Vector< double>::iterator itr( row.Begin()), itr_end( row.End()); itr != itr_end; ++itr)
        {
          *itr /= max_count;
        }

        // replace the row in the histogram with the column-normalized one
        m_Histogram.ReplaceRow( row_index, row);
      }
    }

    //! @brief normalize the column counts in the histogram
    //! each cell is divided by the total number of counts in the column
    void Histogram2D::NormalizeY()
    {
      for( size_t col( 0), col_count( m_Histogram.GetNumberCols()); col < col_count; ++col)
      {
        // copy column from histogram and calculate the count
        linal::Vector< double> col_vec( m_Histogram.GetCol( col));
        const double total_count_in_current_col( col_vec.Sum());

        // if the count if larger than 1, normalize this column
        if( total_count_in_current_col > 1)
        {
          // iterate over each element dividing by the total count of this column
          for( linal::Vector< double>::iterator itr( col_vec.Begin()), itr_end( col_vec.End()); itr != itr_end; ++itr)
          {
            *itr /= total_count_in_current_col;
          }

          // replace the row in the histogram with the column-normalized one
          m_Histogram.ReplaceCol( col, col_vec);
        }
      }
    }

    //! @brief Normalize by another histogram 2d, which represents background counts
    //! @param BACKGROUND background counts histogram
    void Histogram2D::NormalizeByBackground( const Histogram2D &HIST)
    {
      BCL_Assert
      (
        m_LowerUpperBoundariesX == HIST.m_LowerUpperBoundariesX
        && m_LowerUpperBoundariesY == HIST.m_LowerUpperBoundariesY
        && m_BinSizeXY == HIST.m_BinSizeXY,
        "Normalization could not occur because histograms were of different size or binning"
      );
      for( size_t i( 0), nr( HIST.m_Histogram.GetNumberRows()); i < nr; ++i)
      {
        auto this_row( m_Histogram.GetRow( i));
        auto bg_row( HIST.m_Histogram.GetRow( i));
        auto row_threshold( std::max( this_row.Sum() / 100.0 / double( bg_row.GetSize()), 10.5));
        auto itr_bg( bg_row.Begin());
        for
        (
          auto itr_this( this_row.Begin()), itr_this_end( this_row.End());
          itr_this != itr_this_end;
          ++itr_this, ++itr_bg
        )
        {
          // normalize each non-zero bin
          if( *itr_this > row_threshold || *itr_this < -row_threshold)
          {
            *itr_this /= *itr_bg;
          }
          else
          {
            *itr_this = 0.0;
          }
        }
      }
    }

    //! @brief calculate the logarithm of the column counts + 1 in the histogram
    void Histogram2D::Log()
    {
      for( size_t col( 0), col_count( m_Histogram.GetNumberCols()); col < col_count; ++col)
      {
        linal::Vector< double> col_vec( m_Histogram.GetCol( col));
        for( linal::Vector< double>::iterator itr( col_vec.Begin()), itr_end( col_vec.End()); itr != itr_end; ++itr)
        {
          *itr = std::log10( *itr + 1.0);
        }
        m_Histogram.ReplaceCol( col, col_vec);
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write into std::ostream x and y values horizontally
    std::ostream &Histogram2D::Write
    (
      std::ostream &OSTREAM,
      const util::Format &FORMAT_BINNING,
      const util::Format &FORMAT_VALUES
    ) const
    {
      //write header line
      OSTREAM << GetCenterString() << '\t' << GetLeftBoundaryString() << '\t';
      for( size_t i( 0); i < m_Histogram.GetNumberCols() - 2; ++i)
      {
        OSTREAM << GetBinString() << '\t';
      }
      OSTREAM << GetRightBoundaryString() << '\n' << '\t' << GetCenterString() << '\t' << GetCountsString() << '\t';

      //write x values
      //write lower boundary X
      OSTREAM << FORMAT_BINNING( m_LowerUpperBoundariesX.First()) << '\t';

      //write middle point of each bin in X direction
      double i( 0.5);
      for( i = 0.5; i < m_Histogram.GetNumberCols() - 2; ++i)
      {
        OSTREAM << FORMAT_BINNING( m_LowerUpperBoundariesX.First() + i * m_BinSizeXY.First()) << '\t';
      }

      //write upper boundary X
      OSTREAM << FORMAT_BINNING( m_LowerUpperBoundariesX.Second()) << '\n';

      //write lower boundary Y
      OSTREAM << '\t' << GetLeftBoundaryString() << '\t' << FORMAT_BINNING( m_LowerUpperBoundariesY.First()) << '\t';

      //write lower boundary counts
      for( const double *ptr( m_Histogram[ 0]); ptr != m_Histogram[ 0] + m_Histogram.GetNumberCols(); ++ptr)
      {
        OSTREAM << FORMAT_VALUES( *ptr) << '\t';
      }

      OSTREAM << '\n';

      //write all counts
      i = 0.5;
      for( size_t row( 1); row < m_Histogram.GetNumberRows() - 1; ++row, ++i)
      {
        OSTREAM << '\t' << GetBinString() << '\t' << FORMAT_BINNING( m_LowerUpperBoundariesY.First() + i * m_BinSizeXY.Second()) << '\t';
        for( const double *ptr( m_Histogram[ row]); ptr != m_Histogram[ row] + m_Histogram.GetNumberCols(); ++ptr)
        {
          OSTREAM << FORMAT_VALUES( *ptr) << '\t';
        }
        OSTREAM << '\n';
      }

      //write upper boundary Y
      OSTREAM << '\t' << GetRightBoundaryString() << '\t' << FORMAT_BINNING( m_LowerUpperBoundariesY.Second()) << '\t';

      //write upper boundary counts
      for
      (
        const double *ptr( m_Histogram[ m_Histogram.GetNumberRows() - 1]);
        ptr != m_Histogram[ m_Histogram.GetNumberRows() - 1] + m_Histogram.GetNumberCols(); ++ptr
      )
      {
        OSTREAM << FORMAT_VALUES( *ptr) << '\t';
      }

      OSTREAM << '\n';

      return OSTREAM;
    }

    //! read Histogram2D from std::istream
    std::istream &Histogram2D::Read( std::istream &ISTREAM)
    {
      // check whether file is valid Histogram2D file
      std::string identify;
      ISTREAM >> identify;
      BCL_Assert( identify == GetCenterString(), "unexpected string. Should be \"" + GetCenterString() + "\" instead of " + identify);
      ISTREAM >> identify;
      BCL_Assert( identify == GetLeftBoundaryString(), "unexpected string. Should be \"" + GetLeftBoundaryString() + "\" instead of " + identify);

      size_t number_cols( 1), number_rows( 1);
      //count number of cols
      while( !ISTREAM.eof() && identify != GetRightBoundaryString())
      {
        ISTREAM >> identify;
        number_cols++;
        if( identify != GetRightBoundaryString())
        {
          BCL_Assert( identify == GetBinString(), "unexpected string. Should be \"" + GetBinString() + "\"");
        }
      }

      ISTREAM >> identify;
      BCL_Assert( identify == GetCenterString(), "unexpected string. Should be \"" + GetCenterString() + "\"");
      ISTREAM >> identify;
      BCL_Assert( identify == GetCountsString(), "unexpected string. Should be \"" + GetCountsString() + "\"");

      //read lower boundary X
      ISTREAM >> m_LowerUpperBoundariesX.First();

      //read upper boundary X
      while( !ISTREAM.eof())
      {
        ISTREAM >> identify;
        if( identify != GetLeftBoundaryString())
        {
          m_LowerUpperBoundariesX.Second() = util::ConvertStringToNumericalValue< double>( identify);
        }
        else
        {
           break;
        }
      }

      //read lower boundary Y
      ISTREAM >> m_LowerUpperBoundariesY.First();

      //read values
      storage::Vector< double> values;

      std::string tmp;
      while( !ISTREAM.eof())
      {
        do
        {
          ISTREAM >> tmp;
          if( util::IsNumerical( tmp))
          {
            values.PushBack( util::ConvertStringToNumericalValue< double>( tmp));
          }

        } while( !ISTREAM.eof() && tmp != GetBinString() && tmp != GetRightBoundaryString());
        number_rows++;
        if( tmp != GetRightBoundaryString())
        {
          ISTREAM >> tmp;
          continue;
        }
        ISTREAM >> m_LowerUpperBoundariesY.Second();
        break;
      }

      for( size_t i = 0; i < number_cols && !ISTREAM.eof(); ++i)
      {
        ISTREAM >> tmp;
        values.PushBack( util::ConvertStringToNumericalValue< double>( tmp));
      }

      m_Histogram = linal::Matrix< double>( number_rows, number_cols, values);

      m_BinSizeXY.First() = ( m_LowerUpperBoundariesX.Second() - m_LowerUpperBoundariesX.First()) / double( number_cols - 2);
      m_BinSizeXY.Second() = ( m_LowerUpperBoundariesY.Second() - m_LowerUpperBoundariesY.First()) / double( number_rows - 2);

      //return
      return ISTREAM;
    }

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_histogram_3d.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_operations.h"
#include "linal/bcl_linal_vector_operations.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

  //////////
  // data //
  //////////

    //! @brief string to indicate the center
    //! @return string to indicate the center
    const std::string &Histogram3D::GetCenterString()
    {
      static std::string s_center_string( "center");
      return s_center_string;
    }

    //! @brief string to indicate counts
    //! @return string to indicate counts
    const std::string &Histogram3D::GetCountsString()
    {
      static std::string s_counts( "counts");
      return s_counts;
    }

    //! @brief string to indicate the left boundary in output
    //! @return string to indicate the left boundary in output
    const std::string &Histogram3D::GetLeftBoundaryString()
    {
      static std::string s_left_boundary_string( "...<");
      return s_left_boundary_string;
    }

    //! @brief string to indicate a bin
    //! @return string to indicate a bin
    const std::string &Histogram3D::GetBinString()
    {
      static std::string s_bin_string( "<..>");
      return s_bin_string;
    }

    //! @brief string to indicate the right boundary in output
    //! @return string to indicate the right boundary in output
    const std::string &Histogram3D::GetRightBoundaryString()
    {
      static std::string s_right_boundary_string( ">...");
      return s_right_boundary_string;
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Histogram3D::s_Instance
    (
      GetObjectInstances().AddInstance( new Histogram3D())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor - initializes everything to 0
    Histogram3D::Histogram3D() :
      m_LowerUpperBoundariesX( 0.0, 0.0),
      m_LowerUpperBoundariesY( 0.0, 0.0),
      m_LowerUpperBoundariesZ( 0.0, 0.0),
      m_BinSizeXYZ( 0.0, 0.0, 0.0),
      m_Histogram( 0, 0, 0, 0.0)
    {
    }

    //! construct Histogram from X and Y starting values MIN_X_Y, from Pair BINSIZE_X_Y and from NUMBER_OF_BINS_X_Y
    Histogram3D::Histogram3D
    (
      const storage::VectorND< 3, double> &MIN_X_Y_Z,
      const storage::VectorND< 3, double> &BINSIZE_X_Y_Z,
      const storage::VectorND< 3, size_t> &NUMBER_OF_BINS_X_Y_Z,
      const double &INITIAL_VAL
    ) :
      m_LowerUpperBoundariesX( MIN_X_Y_Z.First(), MIN_X_Y_Z.First() + BINSIZE_X_Y_Z.First() * NUMBER_OF_BINS_X_Y_Z.First()),
      m_LowerUpperBoundariesY( MIN_X_Y_Z.Second(), MIN_X_Y_Z.Second() + BINSIZE_X_Y_Z.Second() * NUMBER_OF_BINS_X_Y_Z.Second()),
      m_LowerUpperBoundariesZ( MIN_X_Y_Z.Third(), MIN_X_Y_Z.Third() + BINSIZE_X_Y_Z.Third() * NUMBER_OF_BINS_X_Y_Z.Third()),
      m_BinSizeXYZ( BINSIZE_X_Y_Z),
      m_Histogram( NUMBER_OF_BINS_X_Y_Z.First(), NUMBER_OF_BINS_X_Y_Z.Second(), NUMBER_OF_BINS_X_Y_Z.Third(), INITIAL_VAL)
    {
      SetupBinning();
    }

    //! copy constructor
    Histogram3D *Histogram3D::Clone() const
    {
      return new Histogram3D( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Histogram3D::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! get number of bins X
    size_t Histogram3D::GetNumberOfBinsX() const
    {
      return m_Histogram.NumberLayers();
    }

    //! get number of bins Y
    size_t Histogram3D::GetNumberOfBinsY() const
    {
      return m_Histogram.GetNumberRows();
    }

    //! get number of bins Y
    size_t Histogram3D::GetNumberOfBinsZ() const
    {
      return m_Histogram.GetNumberCols();
    }

    //! get binsize XY
    const storage::VectorND< 3, double> &Histogram3D::GetBinSizeXYZ() const
    {
      return m_BinSizeXYZ;
    }

    //! GetBoundaries X
    const storage::VectorND< 2, double> &Histogram3D::GetBoundariesX() const
    {
      return m_LowerUpperBoundariesX;
    }

    //! GetBoundaries Y
    const storage::VectorND< 2, double> &Histogram3D::GetBoundariesY() const
    {
      return m_LowerUpperBoundariesY;
    }

    //! GetBoundaries Y
    const storage::VectorND< 2, double> &Histogram3D::GetBoundariesZ() const
    {
      return m_LowerUpperBoundariesZ;
    }

    //! GetHistogram
    const Tensor< double> &Histogram3D::GetHistogram() const
    {
      return m_Histogram;
    }

    //! GetHistogram
    Tensor< double> &Histogram3D::GetChangeableHistogram()
    {
      return m_Histogram;
    }

    //! get sum of all counts
    double Histogram3D::GetSumOfAllCounts() const
    {
      return m_Histogram.GetValues().Sum();
    }

    //! GetBinning
    storage::VectorND< 3, linal::Vector< double> > Histogram3D::GetBinningXYZ() const
    {
      return storage::VectorND< 3, linal::Vector< double> >( m_BinningX, m_BinningY, m_BinningZ);
    }

  ////////////////
  // operations //
  ////////////////

    //! reset all counts
    void Histogram3D::Reset()
    {
      m_Histogram.GetValues() = 0.0;
    }

    //! pushback a pair of values to the right position in the histogram
    void Histogram3D::PushBack
    (
      const double &X,
      const double &Y,
      const double &Z,
      const double &WEIGHT
    )
    {
      if( !WEIGHT || !util::IsDefined( X) || !util::IsDefined( Y) || !util::IsDefined( Z))
      {
        return;
      }

      size_t bin_x( 0), bin_y( 0), bin_z( 0);

      //determine field in x-direction
      if( X < m_LowerUpperBoundariesX.First())
      {
        bin_x = 0;
      }
      else if( X >= m_LowerUpperBoundariesX.Second())
      {
        bin_x = m_Histogram.NumberLayers() - 1;
      }
      else
      {
        bin_x = size_t( ( X - m_LowerUpperBoundariesX.First()) / m_BinSizeXYZ.First());
      }

      //determine field in y-direction
      if( Y < m_LowerUpperBoundariesY.First())
      {
        bin_y = 0;
      }
      else if( Y >= m_LowerUpperBoundariesY.Second())
      {
        bin_y = m_Histogram.GetNumberRows() - 1;
      }
      else
      {
        bin_y = size_t( ( Y - m_LowerUpperBoundariesY.First()) / m_BinSizeXYZ.Second());
      }

      //determine field in z-direction
      if( Z < m_LowerUpperBoundariesZ.First())
      {
        bin_z = 0;
      }
      else if( Z >= m_LowerUpperBoundariesZ.Second())
      {
        bin_z = m_Histogram.GetNumberCols() - 1;
      }
      else
      {
        bin_z = size_t( ( Z - m_LowerUpperBoundariesZ.First()) / m_BinSizeXYZ.Third());
      }

      m_Histogram( bin_x, bin_y, bin_z) += WEIGHT;
    }

    //! Use tri-linear interpolation to obtain a value at an arbitrary point
    //! @param X, Y, Z the coordinates of interest
    double Histogram3D::Interpolate( const double &X, const double &Y, const double &Z) const
    {
      if( !util::IsDefined( X) || !util::IsDefined( Y) || !util::IsDefined( Z))
      {
        return util::GetUndefined< double>();
      }

      size_t bin_x_lo( 0), bin_x_hi( 0), bin_y_lo( 0), bin_y_hi( 0), bin_z_lo( 0), bin_z_hi( 0);

      //determine field in x-direction
      if( X < m_BinningX( 0))
      {
        bin_x_lo = 0;
        bin_x_hi = 1;
      }
      else if( X >= m_BinningX.Last())
      {
        bin_x_lo = m_Histogram.NumberLayers() - 2;
        bin_x_hi = m_Histogram.NumberLayers() - 1;
      }
      else
      {
        bin_x_lo = size_t( ( X - m_BinningX( 0)) / m_BinSizeXYZ.First());
        bin_x_hi = std::min( bin_x_lo + 1, m_Histogram.NumberLayers() - 1);
      }

      //determine field in y-direction
      if( Y < m_BinningY( 0))
      {
        bin_y_lo = 0;
        bin_y_hi = 1;
      }
      else if( Y >= m_BinningY.Last())
      {
        bin_y_lo = m_Histogram.GetNumberRows() - 2;
        bin_y_hi = m_Histogram.GetNumberRows() - 1;
      }
      else
      {
        bin_y_lo = size_t( ( Y - m_BinningY( 0)) / m_BinSizeXYZ.Second());
        bin_y_hi = std::min( bin_y_lo + 1, m_Histogram.GetNumberRows() - 1);
      }

      //determine field in z-direction
      if( Z < m_BinningZ( 0))
      {
        bin_z_lo = 0;
        bin_z_hi = 1;
      }
      else if( Z >= m_BinningZ.Last())
      {
        bin_z_lo = m_Histogram.GetNumberCols() - 2;
        bin_z_hi = m_Histogram.GetNumberCols() - 1;
      }
      else
      {
        bin_z_lo = size_t( ( Z - m_BinningZ( 0)) / m_BinSizeXYZ.Third());
        bin_z_hi = std::min( bin_z_lo + 1, m_Histogram.GetNumberCols() - 1);
      }

      const double x_to_lo( ( X - m_BinningX( bin_x_lo)));
      const double x_to_hi( ( m_BinningX( bin_x_hi) - X));
      const double y_to_lo( ( Y - m_BinningY( bin_y_lo)));
      const double y_to_hi( ( m_BinningY( bin_y_hi) - Y));
      const double z_to_lo( ( Z - m_BinningZ( bin_z_lo)));
      const double z_to_hi( ( m_BinningZ( bin_z_hi) - Z));

      const double c00
      (
        m_Histogram( bin_x_lo, bin_y_lo, bin_z_lo) * x_to_hi
        + m_Histogram( bin_x_hi, bin_y_lo, bin_z_lo) * x_to_lo
      );
      const double c01
      (
        m_Histogram( bin_x_lo, bin_y_lo, bin_z_hi) * x_to_hi
        + m_Histogram( bin_x_hi, bin_y_lo, bin_z_hi) * x_to_lo
      );
      const double c10
      (
        m_Histogram( bin_x_lo, bin_y_hi, bin_z_lo) * x_to_hi
        + m_Histogram( bin_x_hi, bin_y_hi, bin_z_lo) * x_to_lo
      );
      const double c11
      (
        m_Histogram( bin_x_lo, bin_y_hi, bin_z_hi) * x_to_hi
        + m_Histogram( bin_x_hi, bin_y_hi, bin_z_hi) * x_to_lo
      );

      const double c0( c00 * y_to_hi + c10 * y_to_lo);
      const double c1( c01 * y_to_hi + c11 * y_to_lo);

      return ( c0 * z_to_hi + c1 * z_to_lo) / ( m_BinSizeXYZ.First() * m_BinSizeXYZ.Second() * m_BinSizeXYZ.Third());
    }

    //! Get the value for the nearest bin
    //! @param X, Y, Z the coordinates of interest
    double Histogram3D::Value( const double &X, const double &Y, const double &Z) const
    {
      if( !util::IsDefined( X) || !util::IsDefined( Y) || !util::IsDefined( Z))
      {
        return util::GetUndefined< double>();
      }

      size_t bin_x( 0), bin_y( 0), bin_z( 0);

      //determine field in x-direction
      if( X < m_BinningX( 0) || !m_BinSizeXYZ.First())
      {
        bin_x = 0;
      }
      else if( X >= m_LowerUpperBoundariesX.Second())
      {
        bin_x = m_Histogram.NumberLayers() - 1;
      }
      else
      {
        bin_x = size_t( ( X - m_LowerUpperBoundariesX.First()) / m_BinSizeXYZ.First());
      }

      //determine field in y-direction
      if( Y < m_BinningY( 0) || !m_BinSizeXYZ.Second())
      {
        bin_y = 0;
      }
      else if( Y >= m_LowerUpperBoundariesY.Second())
      {
        bin_y = m_Histogram.GetNumberRows() - 1;
      }
      else
      {
        bin_y = size_t( ( Y - m_LowerUpperBoundariesY.First()) / m_BinSizeXYZ.Second());
      }

      //determine field in z-direction
      if( Z < m_BinningZ( 0) || !m_BinSizeXYZ.Third())
      {
        bin_z = 0;
      }
      else if( Z >= m_LowerUpperBoundariesZ.Second())
      {
        bin_z = m_Histogram.GetNumberCols() - 1;
      }
      else
      {
        bin_z = size_t( ( Z - m_LowerUpperBoundariesZ.First()) / m_BinSizeXYZ.Third());
      }
      return m_Histogram( bin_x, bin_y, bin_z);
    }

    //! @brief combine this with a given histogram2d by adding up all counts
    //! all parameters have to be identical for this operation to work
    //! @param HISTOGRAM_3D histogram to be added to this one
    //! @return true if it was successful
    bool Histogram3D::Combine( const Histogram3D &HISTOGRAM_3D)
    {
      // allow combining, if this is empty
      if( IsEmpty())
      {
        *this = HISTOGRAM_3D;
        return true;
      }

      // check that the parameters agree
      if
      (
           m_LowerUpperBoundariesX  != HISTOGRAM_3D.m_LowerUpperBoundariesX
        || m_LowerUpperBoundariesY  != HISTOGRAM_3D.m_LowerUpperBoundariesY
        || m_LowerUpperBoundariesZ  != HISTOGRAM_3D.m_LowerUpperBoundariesZ
        || m_BinSizeXYZ             != HISTOGRAM_3D.m_BinSizeXYZ
      )
      {
        return false;
      }

      // add the histogram
      linal::VectorReference< double> ref( m_Histogram.GetValues());
      ref += HISTOGRAM_3D.m_Histogram.GetValues();

      // end
      return true;
    }

    //! checks if there is any count
    bool Histogram3D::IsEmpty() const
    {
      return GetSumOfAllCounts() == 0;
    }

    //! @brief normalize the counts in the histogram
    //! each bin is divided by the total number of counts (even the boundary counts)
    void Histogram3D::Normalize()
    {
      m_Histogram.GetValues().SetToSum( 1.0);
    }

    //! @brief Normalize by another histogram 3d, which represents background counts
    //! @param BACKGROUND background counts histogram
    void Histogram3D::NormalizeByBackground( const Histogram3D &HIST)
    {
      BCL_Assert
      (
        m_LowerUpperBoundariesX == HIST.m_LowerUpperBoundariesX
        && m_LowerUpperBoundariesY == HIST.m_LowerUpperBoundariesY
        && m_LowerUpperBoundariesZ == HIST.m_LowerUpperBoundariesZ
        && m_BinSizeXYZ == HIST.m_BinSizeXYZ,
        "Normalization could not occur because histograms were of different size or binning"
      );
      auto itr_bg( HIST.m_Histogram.GetValues().Begin());
      for
      (
        auto itr_this( m_Histogram.Begin()), itr_this_end( m_Histogram.End());
        itr_this != itr_this_end;
        ++itr_this, ++itr_bg
      )
      {
        // normalize each non-zero bin
        if( *itr_this && *itr_bg)
        {
          *itr_this /= *itr_bg;
        }
        else
        {
          *itr_this = 0.0;
        }
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write into std::ostream x and y values horizontally
    std::ostream &Histogram3D::Write
    (
      std::ostream &OSTREAM,
      const util::Format &FORMAT_BINNING,
      const util::Format &FORMAT_VALUES
    ) const
    {
      // Write as a simple X Y Z Value file (readable in VisIt)
      OSTREAM << "X Y Z value\n";
      for( size_t x_i( 0), x_mx( m_Histogram.NumberLayers()); x_i < x_mx; ++x_i)
      {
        const double x( m_LowerUpperBoundariesX.First() + m_BinSizeXYZ.First() * ( 0.5 + double( x_i)));
        for( size_t y_i( 0), y_mx( m_Histogram.GetNumberRows()); y_i < y_mx; ++y_i)
        {
          const double y( m_LowerUpperBoundariesY.First() + m_BinSizeXYZ.Second() * ( 0.5 + double( y_i)));
          for( size_t z_i( 0), z_mx( m_Histogram.GetNumberCols()); z_i < z_mx; ++z_i)
          {
            const double z( m_LowerUpperBoundariesZ.First() + m_BinSizeXYZ.Third() * ( 0.5 + double( z_i)));
            OSTREAM << FORMAT_BINNING( x) << ' ' << FORMAT_BINNING( y) << ' ' << FORMAT_BINNING( z) << ' '
                    << FORMAT_VALUES( m_Histogram( x_i, y_i, z_i)) << '\n';
          }
        }
      }

      return OSTREAM;
    }

    //! read Histogram3D from std::istream
    std::istream &Histogram3D::Read( std::istream &ISTREAM)
    {
      storage::Vector< double> values;
      storage::Vector< double> x_values;
      storage::Vector< double> y_values;
      storage::Vector< double> z_values;
      std::string tmp;
      std::getline( ISTREAM, tmp);
      if( tmp.empty())
      {
        std::getline( ISTREAM, tmp);
      }
      BCL_Assert( util::StartsWith( tmp, "X Y Z "), "Histogram 3D corrupted. Should have started with X Y Z");
      storage::Vector< double> tmp_values;
      while( std::getline( ISTREAM, tmp))
      {
        tmp_values = util::SplitStringToNumerical< double>( tmp, " ");
        if( tmp_values.IsEmpty())
        {
          break;
        }
        BCL_Assert
        (
          tmp_values.GetSize() == 4,
          "Histogram3D Corrupted. Should have had 4 values on each line, instead had: " + tmp
        );
        values.PushBack( tmp_values( 3));
        if( x_values.IsEmpty())
        {
          x_values.PushBack( tmp_values( 0));
          y_values.PushBack( tmp_values( 1));
          z_values.PushBack( tmp_values( 2));
        }
        else if( tmp_values( 0) > x_values.LastElement())
        {
          x_values.PushBack( tmp_values( 0));
        }
        else if( tmp_values( 1) > y_values.LastElement())
        {
          y_values.PushBack( tmp_values( 1));
        }
        else if( tmp_values( 2) > z_values.LastElement())
        {
          z_values.PushBack( tmp_values( 2));
        }
      }
      m_BinSizeXYZ( 0) = x_values.GetSize() > size_t( 1)
                         ? ( x_values.LastElement() - x_values.FirstElement()) / double( x_values.GetSize() - 1)
                         : 0.0;
      m_BinSizeXYZ( 1) = y_values.GetSize() > size_t( 1)
                         ? ( y_values.LastElement() - y_values.FirstElement()) / double( y_values.GetSize() - 1)
                         : 0.0;
      m_BinSizeXYZ( 2) = z_values.GetSize() > size_t( 1)
                         ? ( z_values.LastElement() - z_values.FirstElement()) / double( z_values.GetSize() - 1)
                         : 0.0;
      m_LowerUpperBoundariesX.First() = x_values( 0) - m_BinSizeXYZ( 0) / 2.0;
      m_LowerUpperBoundariesY.First() = y_values( 0) - m_BinSizeXYZ( 1) / 2.0;
      m_LowerUpperBoundariesZ.First() = z_values( 0) - m_BinSizeXYZ( 2) / 2.0;
      m_LowerUpperBoundariesX.Second() = x_values.LastElement() + m_BinSizeXYZ( 0) / 2.0;
      m_LowerUpperBoundariesY.Second() = y_values.LastElement() + m_BinSizeXYZ( 1) / 2.0;
      m_LowerUpperBoundariesZ.Second() = z_values.LastElement() + m_BinSizeXYZ( 2) / 2.0;
      if( values.GetSize())
      {
        m_Histogram = Tensor< double>( x_values.GetSize(), y_values.GetSize(), z_values.GetSize(), &values( 0));
      }
      else
      {
        m_Histogram = Tensor< double>();
      }
      SetupBinning();

      //return
      return ISTREAM;
    }

    //! @brief setup bins according to the boundaries
    void Histogram3D::SetupBinning()
    {
      m_BinningX = linal::Vector< double>( m_Histogram.NumberLayers(), double( 0.0));
      m_BinningY = linal::Vector< double>( m_Histogram.GetNumberRows(), double( 0.0));
      m_BinningZ = linal::Vector< double>( m_Histogram.GetNumberCols(), double( 0.0));

      //binning in X direction
      double i( 0.5);
      for( double *ptr( m_BinningX.Begin()), *ptr_end( m_BinningX.End()); ptr != ptr_end; ++ptr, ++i)
      {
        *ptr = m_LowerUpperBoundariesX.First() + i * m_BinSizeXYZ.First();
      }

      //binning in Y direction
      i = 0.5;
      for( double *ptr( m_BinningY.Begin()), *ptr_end( m_BinningY.End()); ptr != ptr_end; ++ptr, ++i)
      {
        *ptr = m_LowerUpperBoundariesY.First() + i * m_BinSizeXYZ.Second();
      }

      //binning in Z direction
      i = 0.5;
      for( double *ptr( m_BinningZ.Begin()), *ptr_end( m_BinningZ.End()); ptr != ptr_end; ++ptr, ++i)
      {
        *ptr = m_LowerUpperBoundariesZ.First() + i * m_BinSizeXYZ.Third();
      }
    }

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_histogram.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_function_interface.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    //class Histogram

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Histogram::s_Instance
    (
      GetObjectInstances().AddInstance( new Histogram())
    );

    //! @brief string to indicate the lower boundary in output
    //! @return string to indicate the lower boundary in output
    const std::string &Histogram::GetLowerBoundaryString()
    {
      static std::string s_lower_boundary_string( "...<");
      return s_lower_boundary_string;
    }

    //! @brief string to indicate a bin
    //! @return string to indicate a bin
    const std::string &Histogram::GetBinString()
    {
      static std::string s_bin_string( "<..>");
      return s_bin_string;
    }

    //! @brief string to indicate the upper boundary in output
    //! @return string to indicate the upper boundary in output
    const std::string &Histogram::GetUpperBoundaryString()
    {
      static std::string s_upper_boundary_string( ">...");
      return s_upper_boundary_string;
    }

    //! @brief flag for setting the minimum value of the histogram
    //! @return flag for setting the minimum value of the histogram
    const util::ShPtr< command::FlagInterface> Histogram::GetFlagMin()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "histogram_min",
          "The minimum value to be binned by the histogram.",
          command::Parameter( "minimum", "double which is the minimum value to be binned by the histogram", "0")
        )
      );

      return s_flag;
    }

    //! @brief flag for setting the bin size of the histogram
    //! @return flag for setting the bin size of the histogram
    const util::ShPtr< command::FlagInterface> Histogram::GetFlagBinSize()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "histogram_binsize",
          "The size of the bins in the histogram.",
          command::Parameter( "binsize", "double which is the size of the bins in the histogram", "1")
        )
      );

      return s_flag;
    }

    //! @brief flag for setting the number of bins of the histogram
    //! @return flag for setting the number of bins of the histogram
    const util::ShPtr< command::FlagInterface> Histogram::GetFlagNumberOfBins()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "histogram_numbins",
          "The number of the bins in the histogram.",
          command::Parameter( "number of bins", "size_t which is the number of bins in the histogram", "10")
        )
      );

      return s_flag;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor - initializes all to 0
    Histogram::Histogram() :
      m_LowerUpperBoundaries( 0.0, 0.0),
      m_BinSize( 0.0),
      m_Histogram( 0),
      m_LowerUpperBoundariesCounts( 0.0, 0.0)
    {
    }

    //! @brief construct histogram from starting value MIN, BIN_SIZE and NUMBER_OF_BINS
    //! @param MIN the minimal value representing the left boundary
    //! @param BIN_SIZE the size (width) of one bin
    //! @param NUMBER_OF_BINS the number of bin in the histogram
    Histogram::Histogram
    (
      const double MIN,
      const double BIN_SIZE,
      const size_t NUMBER_OF_BINS
    ) :
      m_LowerUpperBoundaries( MIN, MIN + BIN_SIZE * NUMBER_OF_BINS),
      m_BinSize( BIN_SIZE),
      m_Histogram( NUMBER_OF_BINS),
      m_LowerUpperBoundariesCounts( 0, 0)
    {
    }

    //! @brief construct a histogram given a map of data with bins being the key and counts being the value
    //! @param DATA histogram data with key being the bin and the value being the counts for that bin
    Histogram::Histogram( const storage::Map< double, double> &DATA) :
      m_LowerUpperBoundaries( 0.0, 0.0),
      m_BinSize( 0.0),
      m_Histogram( DATA.GetSize()),
      m_LowerUpperBoundariesCounts( 0.0, 0.0)
    {
      // make sure the data is not empty
      if( DATA.IsEmpty())
      {
        return;
      }

      // get the information about the histogram contained in DATA
      const double smallest_bin( DATA.Begin()->first);
      const double largest_bin( DATA.ReverseBegin()->first);
      const size_t num_bins( DATA.GetSize());
      const double bin_size( num_bins > 1 ? ( largest_bin - smallest_bin) / double( num_bins - 1) : 1.0);

      // set the lower and upper values which are going to be one half the bin size above and below the smallest and
      // largest bins in DATA
      m_LowerUpperBoundaries.First() = smallest_bin - bin_size / 2.0;
      m_LowerUpperBoundaries.Second() = largest_bin + bin_size / 2.0;

      BCL_MessageDbg( "lower upper boundaries " + util::Format()( m_LowerUpperBoundaries));

      m_BinSize = bin_size;

      double *ptr( m_Histogram.Begin());

      // fill m_Histogram
      for
      (
        storage::Map< double, double>::const_iterator itr( DATA.Begin()), itr_end( DATA.End()); itr != itr_end;
        ++itr, ++ptr
      )
      {
        ( *ptr) = itr->second;
      }

      BCL_MessageDbg( "histogram " + util::Format()( *this));
    }

    //! copy constructor
    Histogram *Histogram::Clone() const
    {
      return new Histogram( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Histogram::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! get number of bins
    size_t Histogram::GetNumberOfBins() const
    {
      return m_Histogram.GetSize();
    }

    //! get binsize
    double const &Histogram::GetBinSize() const
    {
      return m_BinSize;
    }

    //! GetBoundaries
    const storage::VectorND< 2, double> &Histogram::GetBoundaries() const
    {
      return m_LowerUpperBoundaries;
    }

    //! GetBoundaries counts
    const storage::VectorND< 2, double> &Histogram::GetBoundariesCounts() const
    {
      return m_LowerUpperBoundariesCounts;
    }

    //! GetHistogram
    const linal::Vector< double> &Histogram::GetHistogram() const
    {
      return m_Histogram;
    }

    //! @brief a vector of the bin coordinates of the histogram
    //! @return the center of each bin in a vector
    linal::Vector< double> Histogram::GetBinning() const
    {
      return linal::FillVector< double>( m_Histogram.GetSize(), m_LowerUpperBoundaries.First() + 0.5 * m_BinSize, m_BinSize);
    }

    //! get sum of all counts
    double Histogram::GetSumOfAllCounts() const
    {
      return m_Histogram.Sum() + m_LowerUpperBoundariesCounts.First() + m_LowerUpperBoundariesCounts.Second();
    }

    //! This functions returns the sum of counts between [MIN, MAX)
    double Histogram::GetCountsInBetween( const double MIN, const double MAX) const
    {
      double sum( 0);

      // add counts below lower limit
      if( MIN < m_LowerUpperBoundaries.First())
      {
        sum += m_LowerUpperBoundariesCounts.First();
      }

      // add counts above higher limit
      if( MAX >= m_LowerUpperBoundaries.Second())
      {
        sum += m_LowerUpperBoundariesCounts.Second();
      }

      // add counts within the given interval
      for
      (
        const double
          *ptr( m_Histogram.Begin() + size_t( ( MIN - m_LowerUpperBoundaries.First()) / m_BinSize)),
          *ptr_end( m_Histogram.Begin() + size_t( ( MAX - m_LowerUpperBoundaries.First()) / m_BinSize));
        ptr != ptr_end;
        ++ptr
      )
      {
        sum += *ptr;
      }

      return sum;
    }

    //! @brief set the count for a bin
    //! @param INDEX index of the bin
    //! @param COUNT the new count
    void Histogram::SetCount( const size_t INDEX, const double COUNT)
    {
      m_Histogram( INDEX) = COUNT;
    }

  ////////////////
  // operations //
  ////////////////

    //! calculate the densitymap from the data and parameters
    Histogram &Histogram::CalculateHistogram( const storage::Vector< double> &VALUES_VECTOR)
    {
      //reset all counts to zero
      Reset();

      //count
      for
      (
        std::vector< double>::const_iterator itr( VALUES_VECTOR.Begin()), itr_end( VALUES_VECTOR.End());
        itr != itr_end; ++itr
      )
      {
        PushBack( *itr);
      }
      return *this;
    }

    //! pushback a value to its position in the histogram
    //! @param VALUE the value to be considered to derive the bin
    //! @param WEIGHT the increment for the bin that belongs to VALUE - 1 by default
    void Histogram::PushBack( const double VALUE, const double WEIGHT)
    {
      //checks that value and weight is defined
      if( !util::IsDefined( VALUE) || !util::IsDefined( WEIGHT))
      {
        return;
      }

      // lower as lower boundary
      if( VALUE < m_LowerUpperBoundaries.First())
      {
        m_LowerUpperBoundariesCounts.First() += WEIGHT;
        return;
      }
      else if( VALUE > m_LowerUpperBoundaries.Second())
      {
        m_LowerUpperBoundariesCounts.Second() += WEIGHT;
        return;
      }

      // determine the nominal bin for this value
      const size_t nominal_bin( ( VALUE - m_LowerUpperBoundaries.First()) / m_BinSize);

      // higher then upper boundary; add to outer counts bin
      if( nominal_bin >= m_Histogram.GetSize())
      {
        if( nominal_bin)
        {
          m_Histogram( nominal_bin - 1) += WEIGHT;
        }
        else
        {
          m_LowerUpperBoundariesCounts.First() += WEIGHT / 2.0;
          m_LowerUpperBoundariesCounts.Second() += WEIGHT / 2.0;
        }
      }
      else
      {
        // within one of the bins
        m_Histogram( nominal_bin) += WEIGHT;
      }
    }

    //! @brief combine this with a given histogram by adding up all counts
    //! all parameters have to be identical for this operation to work
    //! @param HISTOGRAM histogram to be added to this one
    //! @return true if it was successful - that is, when boundaries and bin sizes do match
    bool Histogram::Combine( const Histogram &HISTOGRAM)
    {
      // allow combining, if this is empty
      if( IsEmpty())
      {
        *this = HISTOGRAM;
        return true;
      }

      // check that the parameters agree
      if
      (
           m_LowerUpperBoundaries.First() != HISTOGRAM.m_LowerUpperBoundaries.First()
        || m_LowerUpperBoundaries.Second() != HISTOGRAM.m_LowerUpperBoundaries.Second()
        || m_BinSize != HISTOGRAM.m_BinSize
      )
      {
        BCL_MessageCrt
        (
          "combining histograms with different parameters:\n"
          "left boundary:  " + util::Format()( m_LowerUpperBoundaries.First()) + " != " + util::Format()( HISTOGRAM.m_LowerUpperBoundaries.First()) + "\n" +
          "right boundary: " + util::Format()( m_LowerUpperBoundaries.Second()) + " != " + util::Format()( HISTOGRAM.m_LowerUpperBoundaries.Second()) + "\n" +
          "bin size:       " + util::Format()( m_BinSize) + " != " + util::Format()( HISTOGRAM.m_BinSize)
        );
        return false;
      }

      // add the boundary counts and the histogram
      m_LowerUpperBoundariesCounts.First()  += HISTOGRAM.m_LowerUpperBoundariesCounts.First();
      m_LowerUpperBoundariesCounts.Second() += HISTOGRAM.m_LowerUpperBoundariesCounts.Second();
      m_Histogram                           += HISTOGRAM.m_Histogram;

      // end
      return true;
    }

    //! @brief add a pseudo-count to histogram and boundaries
    //! @param VALUE pseudo-count to be added
    void Histogram::AddPseudoCount( const double VALUE)
    {
      // add pseudo count to boundaries
      m_LowerUpperBoundariesCounts.First() += VALUE;
      m_LowerUpperBoundariesCounts.Second() += VALUE;
      // add pseudo count to boundaries
      m_Histogram += VALUE;
    }

    //! @brief reset all counts
    void Histogram::Reset()
    {
      m_Histogram = 0.0;
      m_LowerUpperBoundariesCounts = storage::VectorND< 2, double>( 0, 0);
    }

    //! @brief checks if there is any count
    //! @return
    bool Histogram::IsEmpty() const
    {
      return GetSumOfAllCounts() == double( 0.0);
    }

    //! @brief determine index of last bin with count == 0 from the back
    //! @param COUNT_THRESHOLD if count > COUNT_THRESHOLD, it is considered to contain information - default = 0.0
    //! @return index of last bin with counts in it - number of bins, if right boundary has counts, 0 if nothing has count
    size_t Histogram::GetIndexOfLastInformationContainingBin( const double COUNT_THRESHOLD) const
    {
      // if the boundary count is filled, the last index with information is the last bin
      if( m_LowerUpperBoundariesCounts.Second() > COUNT_THRESHOLD)
      {
        return GetNumberOfBins();
      }

      size_t last_index_information( 0);
      for( size_t i( 0), number_bins( GetNumberOfBins()); i < number_bins; ++i)
      {
        if( m_Histogram( i) > COUNT_THRESHOLD)
        {
          last_index_information = i;
        }
      }

      // end
      return last_index_information;
    }

    //! @brief determine index of first bin with count != 0 from the front
    //! @param COUNT_THRESHOLD if count > COUNT_THRESHOLD, it is considered to contain information - default = 0.0
    //! @return index of first bin with counts in it - number of bins, 0 if all have counts, number of bins, if all have no counts
    size_t Histogram::GetIndexOfFirstInformationContainingBin( const double COUNT_THRESHOLD) const
    {
      // if the boundary count is filled, the last index with information is the last bin
      if( m_LowerUpperBoundariesCounts.First() > COUNT_THRESHOLD)
      {
        return 0;
      }

      for( size_t i( 0), number_bins( GetNumberOfBins()); i < number_bins; ++i)
      {
        if( m_Histogram( i) > COUNT_THRESHOLD)
        {
          return i;
        }
      }

      // end
      return GetNumberOfBins();
    }

    //! @brief remove the bins before the given index
    //! counts are added to left boundary
    //! @param FIRST_INDEX index of bin before which histogram is cut
    void Histogram::RemoveBinsBeforeIndex( const size_t FIRST_INDEX)
    {
      // assert that the index is smaller than the number of bins
      BCL_Assert
      (
        FIRST_INDEX <= GetNumberOfBins(),
        "Cannot remove bins before index, " + util::Format()( FIRST_INDEX) +
        " since it is larger than the number of bins"
      );

      // add the counts of the to be removed bins to the boundary
      for( size_t i( 0); i != FIRST_INDEX; ++i)
      {
        m_LowerUpperBoundariesCounts.First() += m_Histogram( i);
      }

      // remove the bins
      m_Histogram = linal::Vector< double>( m_Histogram.Begin() + FIRST_INDEX, m_Histogram.End());

      // set new lower boundary
      m_LowerUpperBoundaries.First() = m_LowerUpperBoundaries.First() + double( FIRST_INDEX) * m_BinSize;
    }

    //! @brief reset the bins before the given index
    //! counts in bins and at boundary are set to 0.0
    //! @param FIRST_INDEX index of bin before which histogram is reset
    void Histogram::ResetBinsBeforeIndex( const size_t FIRST_INDEX)
    {
      // assert that the index is smaller than the number of bins
      BCL_Assert
      (
        FIRST_INDEX <= GetNumberOfBins(),
        "Cannot reset bins before index, " + util::Format()( FIRST_INDEX) +
        " since it is larger than the number of bins"
      );

      // reset counts of the bins up to first index
      for( size_t i( 0); i != FIRST_INDEX; ++i)
      {
        m_Histogram( i) = 0.0;
      }

      // set new lower boundary count
      m_LowerUpperBoundariesCounts.First() = 0.0;
    }

    //! @brief remove the bins after the given index
    //! will add the counts from the removed bins to the boundary counts
    void Histogram::RemoveBinsAfterIndex( const size_t LAST_INDEX)
    {
      // do nothing if last index is beyond the number of bins
      if( LAST_INDEX >= GetNumberOfBins())
      {
        return;
      }

      // add the counts of the to be removed bins to the boundary
      for( size_t i( LAST_INDEX + 1), number_bins( GetNumberOfBins()); i != number_bins; ++i)
      {
        m_LowerUpperBoundariesCounts.Second() += m_Histogram( i);
      }

      // remove the bins
      m_Histogram = linal::Vector< double>( LAST_INDEX + 1, m_Histogram.Begin());

      // set new upper boundary
      m_LowerUpperBoundaries.Second() = m_LowerUpperBoundaries.First() + m_Histogram.GetSize() * m_BinSize;
    }

    //! @brief reset the bins after the given index
    //! counts in bins and at boundary are set to 0.0
    //! @param LAST_INDEX index of bin after which histogram is reset
    void Histogram::ResetBinsAfterIndex( const size_t LAST_INDEX)
    {
      // do nothing if last index is beyond the number of bins
      if( LAST_INDEX >= GetNumberOfBins())
      {
        return;
      }

      // reset the counts of the bins to the boundary
      for( size_t i( LAST_INDEX + 1), number_bins( GetNumberOfBins()); i != number_bins; ++i)
      {
        m_Histogram( i) = 0.0;
      }

      // set new upper boundary count
      m_LowerUpperBoundariesCounts.Second() = 0.0;
    }

    //! @brief normalize the counts in the histogram
    //! each bin is divided by the total number of counts (even the boundary counts)
    void Histogram::Normalize()
    {
      // total sum
      double total( m_Histogram.Sum());
      total += m_LowerUpperBoundariesCounts.First();
      total += m_LowerUpperBoundariesCounts.Second();

      // normalize only if total is larger than 0
      if( total == double( 0))
      {
        return;
      }

      // normalize
      m_LowerUpperBoundariesCounts.First() /= total;
      m_LowerUpperBoundariesCounts.Second() /= total;
      m_Histogram /= total;
    }

    //! @brief calculate the mean
    //! @return weighted mean of x
    double Histogram::CalculateMean() const
    {
      double weighted_mean( 0);
      double total( 0);

      double current_bin_center( m_LowerUpperBoundaries.First() + 0.5 * m_BinSize);
      // iterate over counts
      for
      (
        const double *count( m_Histogram.Begin()), *count_end( m_Histogram.End());
        count != count_end;
        ++count, current_bin_center += m_BinSize
      )
      {
        weighted_mean += current_bin_center * ( *count);
        total += *count;
      }

      // end
      return weighted_mean / total;
    }

    //! @brief calculate the standard deviation
    //! @return standard deviation
    double Histogram::CalculateSD() const
    {
      double weighted_square_norm( 0);
      double weighted_mean( 0);
      double total( 0);

      double current_bin_center( m_LowerUpperBoundaries.First() + 0.5 * m_BinSize);
      // iterate over counts
      for
      (
        const double *count( m_Histogram.Begin()), *count_end( m_Histogram.End());
        count != count_end;
        ++count, current_bin_center += m_BinSize
      )
      {
        weighted_mean += current_bin_center * ( *count);
        total += *count;
        weighted_square_norm += Sqr( current_bin_center) * ( *count);
      }

      // calculate mean
      weighted_mean /= total;
      // calculate square norm
      weighted_square_norm /= total;

      // calculate standard deviation and return
      return Sqrt( Absolute( weighted_square_norm - Sqr( weighted_mean)));
    }

    //! @brief extends the histogram in the lower or upper direction
    //! @param NUM_BINS_LOWER the number of bins in the lower direction that the histogram should be extended
    //! @param LOWER_BIN_VALUES value that will be assigned to all of the newly extended bins in the lower direction
    //! @param FLOOR the desired lowest value for the lower boundary. bins won't be extended lower than this
    //! @param NUM_BINS_UPPER the number of bins in the upper direction that the histogram should be extended
    //! @param UPPER_BIN_VALUES value that will be assigned to all of the newly extended bins in the upper direction
    //! @param CEILING the desired highest value for the upper boundary. bins won't be extended higher than this
    //! @return bool indicating of extension was successful. Might fail if there are boundary counts since extension
    //!         might put the boundary counts into a bin but don't know from boundary counts alone.
    bool Histogram::ExtendBoundaries
    (
      const size_t NUM_BINS_LOWER,
      const double LOWER_BIN_VALUES,
      const double FLOOR,
      const size_t NUM_BINS_UPPER,
      const double UPPER_BIN_VALUES,
      const double CEILING
    )
    {
      // can't work if there are counts in the upper or lower boundaries
      if( m_LowerUpperBoundariesCounts.First() != 0 || m_LowerUpperBoundariesCounts.Second() != 0)
      {
        BCL_MessageStd( "can't extend boundaries of histogram with boundary counts");
        return false;
      }

      // make sure the floor and ceiling are outside of the boundaries
      if( FLOOR > m_LowerUpperBoundaries.First() || CEILING < m_LowerUpperBoundaries.Second())
      {
        BCL_MessageStd( "floor or ceiling within boundaries of histogram");
        return false;
      }

      // check how many bins in lower and upper direction the histogram can be extended
      const size_t lower_bin_room( Absolute( m_LowerUpperBoundaries.First() - FLOOR)    / m_BinSize);
      const size_t upper_bin_room( Absolute( CEILING - m_LowerUpperBoundaries.Second()) / m_BinSize);

      // set the number of bins to extend in lower direction
      const size_t lower_additional_bins( std::min( lower_bin_room, NUM_BINS_LOWER));

      // set the number of bins to extend in upper direction
      const size_t upper_additional_bins( std::min( upper_bin_room, NUM_BINS_UPPER));

      // set the upper and lower boundaries
      m_LowerUpperBoundaries.First() = m_LowerUpperBoundaries.First() - lower_additional_bins * m_BinSize;
      m_LowerUpperBoundaries.Second() = m_LowerUpperBoundaries.Second() + upper_additional_bins * m_BinSize;

      // add the new lower bins
      linal::Vector< double> new_values( lower_additional_bins + m_Histogram.GetSize() + upper_additional_bins);
      new_values.ReplaceElements( 0, linal::Vector< double>( lower_additional_bins, LOWER_BIN_VALUES));
      new_values.ReplaceElements( lower_additional_bins, m_Histogram);
      new_values.ReplaceElements
      (
        lower_additional_bins + m_Histogram.GetSize(),
        linal::Vector< double>( upper_additional_bins, UPPER_BIN_VALUES)
      );

      // set m_Histogram
      m_Histogram = new_values;

      return true;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write into std::ostream x and y values vertically
    //! @param OSTREAM output stream
    //! @param FORMAT_BINNING format for the bin middle coordinate
    //! @param FORMAT_VALUES format for the written values (counts)
    //! @param INDENT the indentation
    //! @return the ostream written to
    std::ostream &Histogram::WriteHorizontally
    (
      std::ostream &OSTREAM,
      const util::Format &FORMAT_BINNING,
      const util::Format &FORMAT_VALUES,
      const size_t INDENT
    ) const
    {
      //write header line
      io::Serialize::InsertIndent( OSTREAM, INDENT);
      OSTREAM << "\t\t" << GetLowerBoundaryString() << '\t';
      for( size_t i( 0); i < m_Histogram.GetSize(); ++i)
      {
        OSTREAM << GetBinString() << '\t';
      }
      OSTREAM << GetUpperBoundaryString() << '\n';
      io::Serialize::InsertIndent( OSTREAM, INDENT);
      OSTREAM << "center\t\t";

      //write x values
      //write lower boundary
      OSTREAM << FORMAT_BINNING( m_LowerUpperBoundaries.First()) << '\t';

      //write middle point of each bin
      for( double i( 0.5); i < m_Histogram.GetSize(); ++i)
      {
        OSTREAM << FORMAT_BINNING( m_LowerUpperBoundaries.First() + i * m_BinSize) << '\t';
      }

      //write upper boundary
      OSTREAM << FORMAT_BINNING( m_LowerUpperBoundaries.Second()) << '\n';
      io::Serialize::InsertIndent( OSTREAM, INDENT);
      OSTREAM << "counts\t\t";

      //write all counts
      OSTREAM << FORMAT_VALUES( m_LowerUpperBoundariesCounts.First()) << '\t';
      for( const double *ptr( m_Histogram.Begin()), *ptr_end( m_Histogram.End()); ptr != ptr_end; ++ptr)
      {
        OSTREAM << FORMAT_VALUES( *ptr) << '\t';
      }
      OSTREAM << FORMAT_VALUES( m_LowerUpperBoundariesCounts.Second()) << '\n';

      return OSTREAM;
    }

    //! write into std::ostream x and y values vertically
    std::ostream &Histogram::WriteVertically( std::ostream &OSTREAM, const util::Format &FORMAT_BINNING, const util::Format &FORMAT_VALUES) const
    {
      //write header line
      OSTREAM << "center\tcounts" << '\n';

      //write lower boundary
      OSTREAM << GetLowerBoundaryString() << "\t\t" << FORMAT_BINNING( m_LowerUpperBoundaries.First()) << '\t' << FORMAT_VALUES( m_LowerUpperBoundariesCounts.First()) << '\n';

      //write middle point of each bin and the according count
      double i( 0.5);
      for( const double *ptr( m_Histogram.Begin()), *ptr_end( m_Histogram.End()); ptr != ptr_end; ++ptr, ++i)
      {
        OSTREAM << GetBinString() << "\t\t"
               << FORMAT_BINNING( m_LowerUpperBoundaries.First() + i * m_BinSize) << '\t'
               << FORMAT_VALUES( *ptr) << '\n';
      }

      //write upper boundary
      OSTREAM << GetUpperBoundaryString() << "\t\t" << FORMAT_BINNING( m_LowerUpperBoundaries.Second()) << '\t' << FORMAT_VALUES( m_LowerUpperBoundariesCounts.Second()) << '\n';

      return OSTREAM;
    }

    //! read horizontally written Histogram from std::istream
    std::istream &Histogram::Read( std::istream &ISTREAM)
    {
      //return
      return ReadHorizontally( ISTREAM);
    }

    //! @brief write the data of the histogram in gnuplot format
    //! @param OSTREAM the stream to write gnuplot script to
    //! @param BINNING the data of the histgram that is going to be printed in gnuplot format
    //! @param COUNTS the counts that correspond to the bins
    //! @return the stream it was written to
    std::ostream &Histogram::WriteGnuPlotHeatMapFormatted
    (
      std::ostream &OSTREAM, const linal::Vector< double> &BINNING, const linal::Vector< double> &COUNTS
    )
    {
      for
      (
        const double *x( BINNING.Begin()), *x_end( BINNING.End()), *y( COUNTS.Begin()), *y_end( COUNTS.End());
        x != x_end && y != y_end;
        ++x, ++y
      )
      {
        OSTREAM << '\n';
        for( size_t i( 0); i < 2; ++i)
        {
          OSTREAM << i << '\t' << *x << '\t' << *y << '\n';
        }
      }

      // indicate end of data for heatmap
      OSTREAM << "e\n";

      // end
      return OSTREAM;
    }

    //! @brief write the data of the histogram in gnuplot format for making a line plot
    //! @param OSTREAM the stream to write gnuplot script to
    //! @param BINNING the data of the histgram that is going to be printed in gnuplot format
    //! @param COUNTS the counts that correspond to the bins
    //! @return the stream it was written to
    std::ostream &Histogram::WriteGnuPlotLinePlotFormatted
    (
      std::ostream &OSTREAM, const linal::Vector< double> &BINNING, const linal::Vector< double> &COUNTS
    )
    {
      // write out the data to file
      for
      (
        const double *x( BINNING.Begin()), *x_end( BINNING.End()), *y( COUNTS.Begin()), *y_end( COUNTS.End());
        x != x_end && y != y_end;
        ++x, ++y
      )
      {
        OSTREAM << *x << '\t' << *y << '\n';
      }

      // indicate end of data for gnuplot
      OSTREAM << "e\n";

      // end
      return OSTREAM;
    }

    //! @brief generate a linear gnuplot
    //! @param OSTREAM the stream to write gnuplot script to
    //! @param TITLE title for the linear map
    //! @return the stream it was written to
    std::ostream &Histogram::WriteLinearGnuplot( std::ostream &OSTREAM, const std::string &TITLE) const
    {
      // write comments
      OSTREAM << "# BCL generated linear plot from histogram\n";

      OSTREAM << "set terminal png transparent enhanced # size 2160,640 \n";
      OSTREAM << "set output \"" << TITLE << ".png\"\n";
      OSTREAM << "set encoding iso\n";
      OSTREAM << "set view map\n";
      OSTREAM << "set title \"" << TITLE << "\"\n";
      OSTREAM << "unset key\n\n";

      OSTREAM << "set xlabel \"x\"\n";
      OSTREAM << "set xrange [" << m_LowerUpperBoundaries.First() << ":" << m_LowerUpperBoundaries.Second() << "]\n";
      OSTREAM << "set autoscale y\n";
      OSTREAM << "plot '-' using 1:2 with lines\n";
      OSTREAM << "# number x values " << m_Histogram.GetSize() << " binning " << m_BinSize << '\n';
      OSTREAM << '\n';

      // write histogram data and return OSTREAM
      return WriteGnuPlotLinePlotFormatted( OSTREAM, GetBinning(), m_Histogram);
    }

    //! @brief read horizontally written Histogram from std::istream
    //! @param ISTREAM inout stream to read from
    //! @return the input stream read from
    std::istream &Histogram::ReadHorizontally( std::istream &ISTREAM)
    {
      // check whether file is valid Histogram file
      std::string identify;
      ISTREAM >> identify;
      BCL_Assert( identify == GetLowerBoundaryString(), "BCL_HISTOGRAM is not written horizontally");

      size_t number_of_bins( 0);

      // count the number of bins
      do
      {
        ISTREAM >> identify;
        ++number_of_bins;
      } while( !ISTREAM.eof() && identify != "center");

      // decrement for the boundary bins
      number_of_bins -= 2;

      //read lower boundary
      ISTREAM >> m_LowerUpperBoundaries.First();
      do
      {
        ISTREAM >> identify;
        if( util::IsNumerical( identify))
        {
          m_LowerUpperBoundaries.Second() = util::ConvertStringToNumericalValue< double>( identify);
        }
      } while( !ISTREAM.eof() && identify != "counts");

      //set histogram to new size
      m_Histogram = linal::Vector< int>( number_of_bins);

      //read count of lower boundary
      ISTREAM >> m_LowerUpperBoundariesCounts.First();

      for( double *ptr( m_Histogram.Begin()), *ptr_end( m_Histogram.End()); ptr != ptr_end; ++ptr)
      {
        ISTREAM >> *ptr;
      }

      // read count of upper boundary
      ISTREAM >> m_LowerUpperBoundariesCounts.Second();

      // calculate bin size
      m_BinSize = ( m_LowerUpperBoundaries.Second() - m_LowerUpperBoundaries.First()) / number_of_bins;

      //return
      return ISTREAM;
    }

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_kernel_function.h"

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
    const util::SiPtr< const util::ObjectInterface> KernelFunction::s_Instance
    (
      GetObjectInstances().AddInstance( new KernelFunction())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    KernelFunction::KernelFunction()
    {
    }

    //! @brief construct from properties
    KernelFunction::KernelFunction
    (
      const Range< double> &KERNEL
    ) :
      m_Kernel( KERNEL)
    {
    }

    //! @brief Clone function
    //! @return pointer to new KernelFunction
    KernelFunction *KernelFunction::Clone() const
    {
      return new KernelFunction( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &KernelFunction::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief access to range of kernel
    //! @return the range of the kernel
    Range< double> KernelFunction::GetKernel() const
    {
      return m_Kernel;
    }

    //! @brief access to range of kernel
    //! @param KERNEL the new range of the kernel
    void KernelFunction::SetKernel( const Range< double> &KERNEL)
    {
      m_Kernel = KERNEL;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator f( x) = y = 1 if x is in the kernel range
    //! @param ARGUMENT the x
    //! @return f( x)
    double KernelFunction::operator()( const double &ARGUMENT) const
    {
      return ( m_Kernel.IsWithin( ARGUMENT)) ? 1.0 : 0.0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &KernelFunction::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &KernelFunction::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_linear_least_squares.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_inversion_moore_penrose.h"
#include "storage/bcl_storage_pair.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> LinearLeastSquares::s_Instance
    (
      GetObjectInstances().AddInstance( new LinearLeastSquares())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LinearLeastSquares::LinearLeastSquares()
    {
    }

    //! @brief Clone is the virtual copy constructor
    LinearLeastSquares *LinearLeastSquares::Clone() const
    {
      return new LinearLeastSquares( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &LinearLeastSquares::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief SolutionAndChiSquared finds the solutions as well as the chi-squared
    //! chi-squared is the sum of residuals squared divided by the sum of of the y's squared
    storage::Pair< linal::Vector< double>, double> LinearLeastSquares::SolutionAndChiSquared
    (
      const linal::MatrixConstInterface< double> &X,
      const linal::VectorConstInterface< double> &Y
    )
    {
      // make a reference to X or X_Transposed, whichever has more rows
      // this way, the user does not have to remember whether X needs to be tall or wide, the function works either way
      if( X.GetNumberCols() > X.GetNumberRows())
      {
        return SolutionAndChiSquared( X.Transposed(), Y);
      }

      // Get the solutions
      linal::MatrixInversionMoorePenrose< double> inverter( X);
      const linal::Vector< double> solutions( inverter.Solve( Y));

      // calculate the sum of the solutions squared
      const double sum_solutions_squared( std::inner_product( Y.Begin(), Y.End(), Y.Begin(), 0.0));

      // calculate the sum of the residuals.
      // a residual is the experimental value minus the predicted value squared
      // the predicted value is the row of the matrix times the solutions
      double sum_residuals( 0.0);
      for( size_t index( 0), y_size( Y.GetSize()); index < y_size; index++)
      {
        // calculate the value given by the linear-least squares fit to the data
        const double linear_estimate( std::inner_product( solutions.Begin(), solutions.End(), X[ index], 0.0));

        // add the difference between the estimate and the actual value, squared, to the residual sum
        sum_residuals += std::pow( Y( index) - linear_estimate, 2.0);
      }

      return storage::Pair< linal::Vector< double>, double>
             (
               solutions,
               sum_residuals / sum_solutions_squared
             );
    }

    //! @brief SolutionAndChiSquared finds the solutions as well as the chi-squared
    //! chi-squared is the sum of residuals squared divided by the sum of of the y's squared
    storage::Pair< linal::Vector< float>, float> LinearLeastSquares::SolutionAndChiSquared
    (
      const linal::MatrixConstInterface< float> &X,
      const linal::VectorConstInterface< float> &Y
    )
    {
      // make a reference to X or X_Transposed, whichever has more rows
      // this way, the user does not have to remember whether X needs to be tall or wide, the function works either way
      if( X.GetNumberCols() > X.GetNumberRows())
      {
        return SolutionAndChiSquared( X.Transposed(), Y);
      }

      // Get the solutions
      linal::MatrixInversionMoorePenrose< float> inverter( X);
      const linal::Vector< float> solutions( inverter.Solve( Y));

      // calculate the sum of the solutions squared
      const float sum_solutions_squared( std::inner_product( Y.Begin(), Y.End(), Y.Begin(), 0.0));

      // calculate the sum of the residuals.
      // a residual is the experimental value minus the predicted value squared
      // the predicted value is the row of the matrix times the solutions
      float sum_residuals( 0.0);
      for( size_t index( 0), y_size( Y.GetSize()); index < y_size; index++)
      {
        // calculate the value given by the linear-least squares fit to the data
        const float linear_estimate( std::inner_product( solutions.Begin(), solutions.End(), X[ index], 0.0));

        // add the difference between the estimate and the actual value, squared, to the residual sum
        sum_residuals += std::pow( Y( index) - linear_estimate, float( 2.0));
      }

      return storage::Pair< linal::Vector< float>, float>
             (
               solutions,
               sum_residuals / sum_solutions_squared
             );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LinearLeastSquares::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &LinearLeastSquares::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  } // namespace math
} // namespace bcl

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
#include "math/bcl_math_log_likelihood.h"

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
    const util::SiPtr< const util::ObjectInterface> LogLikelihood::s_Instance
    (
      GetObjectInstances().AddInstance( new LogLikelihood())
    );

    //! an arbitary offset to the score, that shifts the confidence level, above which scores are returned negative
    //! in Tyka , DiMaio et.al. 2009 it was set to 0.5, which gives negative values for z-score > -1.0
    const double LogLikelihood::s_ArbitaryScoreConfidenceOffset( 0.5);

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LogLikelihood::LogLikelihood() :
      m_ZScore()
    {
    }

    //! @brief construct with Mean and Standard Deviation
    //! @param MEAN mean value of the considered distribution
    //! @param SIGMA standard deviation of the considered distribution
    LogLikelihood::LogLikelihood
    (
      const double MEAN,
      const double SIGMA
    ) :
      m_ZScore( MEAN, SIGMA)
    {
    }

    //! @brief Clone function
    //! @return pointer to new LogLikelihood
    LogLikelihood *LogLikelihood::Clone() const
    {
      return new LogLikelihood( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LogLikelihood::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &LogLikelihood::GetScheme() const
    {
      static const std::string s_scheme( "log(0.5*(1-erf((X-MEAN)/SD)))");
      return s_scheme;
    }

    //! @brief access to Mean
    //! @return the mean used
    double LogLikelihood::GetMean() const
    {
      return m_ZScore.GetMean();
    };

    //! @brief access to Standard Deviation
    //! @return the sigma used
    double LogLikelihood::GetSigma() const
    {
      return m_ZScore.GetSigma();
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculates the score of the value in the normal distribution with the given mean and sigma
    //! @param ARGUMENT value which z-score shall be calculated from
    //! @return score of the given value
    double LogLikelihood::operator()( const double &ARGUMENT) const
    {
      // make sure, x never reaches 1 to avoid -inf as a result
      const double x( erf( m_ZScore( ARGUMENT)) - 0.000001);
      return log( s_ArbitaryScoreConfidenceOffset * ( 1.0 - x));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LogLikelihood::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_ZScore, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &LogLikelihood::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_ZScore, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_minus_equals.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API MinusEquals< double>;
    template class BCL_API MinusEquals< float>;

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_mod_equals.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API ModEquals< double>;
    template class BCL_API ModEquals< float>;

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_mutate_vector.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_mutate_result.h"
#include "random/bcl_random_uniform_distribution.h"

// external includes - sorted alphabetically
#include <set>

namespace bcl
{
  namespace math
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateVector::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateVector())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateVector::MutateVector() :
      m_Constant( util::GetUndefinedSize_t())
    {
    }

    //! @brief constructor from mutate range, number params mutated and a boolean for using percentages
    //! @param NUMBER_PARAMS Total number of parameters
    //! @param MUTATE_RANGE the range allowed for all parameters
    //! @param NUMBER_PARAMS_MUTATED the number of params mutated at each call to the operator function
    //! @param USE_PERCENTAGES boolean to indicate whether the mutate ranges are percentages
    //! @param KEEP_POSITIVE keep the values as positive numbers
    //! @param CONSTANT index to keep constant. By default, no such index is defined
    MutateVector::MutateVector
    (
      const size_t NUMBER_PARAMS,
      const double MUTATE_RANGE,
      const size_t NUMBER_PARAMS_MUTATED,
      const bool USE_PERCENTAGE,
      const bool KEEP_POSITIVE,
      const size_t &CONSTANT
    ) :
      m_Range( NUMBER_PARAMS, MUTATE_RANGE),
      m_NumberParamsMutated( NUMBER_PARAMS_MUTATED),
      m_UsePercentage( USE_PERCENTAGE),
      m_KeepPositive( KEEP_POSITIVE),
      m_Constant( CONSTANT)
    {
      BCL_Assert
      (
        NUMBER_PARAMS >= NUMBER_PARAMS_MUTATED,
        "Request to mutate " + util::Format()( NUMBER_PARAMS_MUTATED) + " parameters but only " +
        util::Format()( NUMBER_PARAMS) + " parameters exist"
      );

      // make sure the given range is non-zero
      BCL_Assert( MUTATE_RANGE != double( 0.0), "given range is 0.0");
    }

    //! @brief constructor from mutate range, number params mutated and a boolean for using percentages
    //! @param MUTATE_RANGES a vector of acceptable ranges for each parameter
    //! @param NUMBER_PARAMS_MUTATED the number of params mutated at each call to the operator function
    //! @param USE_PERCENTAGE boolean to indicate whether the mutate ranges are percentages
    //! @param KEEP_POSITIVE keep the values as positive numbers
    //! @param CONSTANT index to keep constant. By default, no such index is defined
    MutateVector::MutateVector
    (
      const linal::Vector< double> &MUTATE_RANGES,
      const size_t NUMBER_PARAMS_MUTATED,
      const bool USE_PERCENTAGE,
      const bool KEEP_POSITIVE,
      const size_t &CONSTANT
    ) :
      m_Range( MUTATE_RANGES),
      m_NumberParamsMutated( NUMBER_PARAMS_MUTATED),
      m_UsePercentage( USE_PERCENTAGE),
      m_KeepPositive( KEEP_POSITIVE),
      m_Constant( CONSTANT)
    {
      BCL_Assert
      (
        MUTATE_RANGES.GetSize() >= NUMBER_PARAMS_MUTATED,
        "Request to mutate " + util::Format()( NUMBER_PARAMS_MUTATED) + " parameters but only " +
        util::Format()( MUTATE_RANGES.GetSize()) + " parameters exist"
      );

      size_t mutates_not_zero( 0);
      for( const double *ptr( m_Range.Begin()), *ptr_end( m_Range.End()); ptr != ptr_end; ++ptr)
      {
        if( *ptr != double( 0.0))
        {
          ++mutates_not_zero;
        }
      }

      // check that there are larger equal ranges not zero than parameters to be mutated
      BCL_Assert
      (
        mutates_not_zero >= NUMBER_PARAMS_MUTATED,
        "there are less ranges different from zero that parameters to be mutated: " +
        util::Format()( mutates_not_zero) + " < " + util::Format()( NUMBER_PARAMS_MUTATED)
      );
    }

    //! @brief Clone function
    //! @return pointer to new MutateVector
    MutateVector *MutateVector::Clone() const
    {
      return new MutateVector( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateVector::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator taking an VECTOR and returning a mutated Vector
    //! @param VECTOR Vector of parameters to be mutated
    //! @return MutateResult containing the mutated Vector
    MutateResult< linal::Vector< double> > MutateVector::operator()( const linal::Vector< double> &VECTOR) const
    {
      // indices that indicate which parameters to mutate
      std::set< size_t> indices;

      const size_t have_constant( m_Constant < VECTOR.GetSize() ? 1 : 0);
      if( m_NumberParamsMutated >= m_Range.GetSize() - have_constant)
      {
        // mutate all parameters being optimized
        for( size_t i( 0); i < m_Range.GetSize(); ++i)
        {
          // only insert index if param is actually to be mutated given by range
          if( m_Range( i) != double( 0.0) && i != m_Constant)
          {
            indices.insert( i);
          }
        }
      }
      else
      {
        // generate random indices to select which parameters to mutate
        while( indices.size() < m_NumberParamsMutated)
        {
          // select a value in the desired range; if there is a constant, the range is effectively one smaller
          // if the value ends up being >= the constant, we add 1 then to bring it back to the original range of
          // interest
          size_t chosen_value
          (
            random::GetGlobalRandom().Random< size_t>( 0, m_Range.GetSize() - 1 - have_constant)
          );
          if( chosen_value >= m_Constant)
          {
            ++chosen_value;
          }
          indices.insert( chosen_value);
        }
      }

      // create the vector of mutated parameters to return
      // initially a copy of original parameters - we overwrite the ones we mutate
      util::ShPtr< linal::Vector< double> > sp_mutated_vector( VECTOR.Clone());

      // loop over and mutate all parameters indicated by indices
      for
      (
        std::set< size_t>::const_iterator index_itr( indices.begin()), index_itr_end( indices.end());
        index_itr != index_itr_end;
        ++index_itr
      )
      {
        // get the index and look up the associated mutate range
        const size_t current_index( *index_itr);
        const double current_range( m_Range( current_index));

        if( m_KeepPositive)
        {
          // using random number generator and the mutate range, calculate the change to be applied
          double change( random::GetGlobalRandom().Random< double>( 0.0, 1.0));
          if( change < 0.45)
          {
            sp_mutated_vector->operator()( current_index) *= 0.75;
          }
          else
          {
            sp_mutated_vector->operator()( current_index) *= 1.333333;
          }
        }
        else
        {
          // using random number generator and the mutate range, calculate the change to be applied
          double change
          (
            random::GetGlobalRandom().Random< double>( -current_range, current_range)
          );

          // if percentages is being used then multiply the value with the percentage
          if( m_UsePercentage)
          {
            change *= sp_mutated_vector->operator()( current_index);
          }
          // change the value at the current_index correspondingly
          sp_mutated_vector->operator()( current_index) += change;
        }
      }

      // return the result
      return MutateResult< linal::Vector< double> >( sp_mutated_vector, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateVector::Read( std::istream &ISTREAM)
    {
      // read member variables
      io::Serialize::Read( m_Range              , ISTREAM);
      io::Serialize::Read( m_NumberParamsMutated, ISTREAM);
      io::Serialize::Read( m_UsePercentage      , ISTREAM);
      io::Serialize::Read( m_KeepPositive       , ISTREAM);
      io::Serialize::Read( m_Constant           , ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &MutateVector::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member variables
      io::Serialize::Write( m_Range              , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NumberParamsMutated, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_UsePercentage      , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_KeepPositive       , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Constant           , OSTREAM, INDENT) << '\n';

      // end
      return OSTREAM;
    }

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_piecewise_function.h"

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
    const util::SiPtr< const util::ObjectInterface> PiecewiseFunction::s_Instance
    (
      GetObjectInstances().AddInstance( new PiecewiseFunction())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PiecewiseFunction::PiecewiseFunction() :
      m_RangesFunctions()
    {
    }

    //! @brief construct from the ranges and their functions the piecewise function
    //! @param RANGES_FUNCTIONS List which gives the ranges and their functions as Pairs
    PiecewiseFunction::PiecewiseFunction
    (
      const storage::List
      <
        storage::Pair< Range< double>, util::ShPtr< FunctionInterfaceSerializable< double, double> > >
      > &RANGES_FUNCTIONS
    ) :
      m_RangesFunctions( RANGES_FUNCTIONS)
    {
      // check to make sure the ranges and functions define a valid piecewise function (i.e. no overlapping ranges
      // or no non existing functions)
      BCL_Assert
      (
        RangesFunctionsValid( m_RangesFunctions),
        "The ShPtr provided does not point to anything or the ranges overlap."
      );
    }

    //! @brief Clone function
    //! @return pointer to new class_name
    PiecewiseFunction *PiecewiseFunction::Clone() const
    {
      return new PiecewiseFunction( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PiecewiseFunction::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief GetRangesFunctions returns m_RangesFunctions
    //! @return returns List which is m_RangesFunctions
    const storage::List
    <
      storage::Pair< Range< double>, util::ShPtr< FunctionInterfaceSerializable< double, double> > >
    > &PiecewiseFunction::GetRangesFunctions() const
    {
      return m_RangesFunctions;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator() taking an x-value and returning the y-value based on "m_RangesFunctions"
    //! @param X_VALUE double which is the x value from which the y-value will be calculated
    //! @return returns a double which is the y-value based on X_VALUE and "m_RangesFunctions"
     double PiecewiseFunction::operator()( const double &X_VALUE) const
     {
       // iterate through the list of ranges and functions to see which range the X_VALUE lies within
       for
       (
         storage::List
         <
           storage::Pair< Range< double>, util::ShPtr< FunctionInterfaceSerializable< double, double> > >
         >::const_iterator iter( m_RangesFunctions.Begin()), iter_end( m_RangesFunctions.End());
         iter != iter_end; ++iter
       )
       {
         const Range< double> &range_a( iter->First());
         const util::ShPtr< FunctionInterfaceSerializable< double, double> > &function( iter->Second());

         // true if X_VALUE is within the range currently denoted by "iter"
         if( range_a.IsWithin( X_VALUE))
         {
           // execute the function to get the corresponding function value
           return function->operator()( X_VALUE);
         }
       }

       // if this point is reached then "X_VALUE" does not fall in any of the ranges so return undefined double
       return util::GetUndefined< double>();
     }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PiecewiseFunction::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_RangesFunctions, ISTREAM);

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PiecewiseFunction::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_RangesFunctions, OSTREAM, INDENT);

      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief determine whether the ranges overlap or if the ShPtr of the functions point to anything
    //! @param RANGES_FUNCTIONS_VALID gives a List with a pair of ranges and function interfaces to be analyzed
    //! @return bool which will tell whether or not the list is valid
    bool PiecewiseFunction::RangesFunctionsValid
    (
      const storage::List
      <
        storage::Pair< Range< double>, util::ShPtr< FunctionInterfaceSerializable< double, double> > >
      > &RANGES_FUNCTIONS_VALID
    )
    {
      // iterate through the list "RANGES_FUNCTIONS_VALID"
      for
      (
        storage::List
        <
          storage::Pair< Range< double>, util::ShPtr< FunctionInterfaceSerializable< double, double> > >
        >::const_iterator iter_a( RANGES_FUNCTIONS_VALID.Begin()), iter_end( RANGES_FUNCTIONS_VALID.End());
        iter_a != iter_end; ++iter_a
      )
      {
        const Range< double> &range_a( iter_a->First());
        const util::ShPtr< FunctionInterfaceSerializable< double, double> > &function( iter_a->Second());

        // for each function determine if it exists
        if( !function.IsDefined())
        {
          BCL_MessageCrt( "The function does not exist.")
          return false;
        }

        // iterate through the list "RANGES_FUNCTIONS_VALID" again and
        // for each range determine if it overlaps with any of the other ranges
        for
        (
          storage::List
          <
            storage::Pair< Range< double>, util::ShPtr< FunctionInterfaceSerializable< double, double> > >
          >::const_iterator iter_b( iter_a); iter_b != iter_end; ++iter_b
        )
        {
          const Range< double> &range_b( iter_b->First());

          // ignore the ranges that are equal to themselves
          if( iter_b == iter_a)
          {
            continue;
          }

          // true if "range_a" overlaps with "range_b" or vice versa
          if( range_a.DoesOverlap( range_b))
          {
            BCL_MessageCrt
            (
              "The following ranges overlap: "
                + range_a.GetString( util::Format()) + ", " + range_b.GetString( util::Format()) + "."
            );

            // return false indicating that the list of ranges and functions is not valid
            return false;
          }
        }
      }

      // return true indicating that none of the ranges overlap and all of the functions are defined
      return true;
    }

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_plus_equals.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API PlusEquals< double>;
    template class BCL_API PlusEquals< float>;

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_polynomial.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serialize.h"
#include "math/bcl_math_linear_least_squares.h"
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> Polynomial::s_Instance
    (
      util::Enumerated< Polynomial>::AddInstance( new Polynomial())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Polynomial::Polynomial()
    {
      UpdateScheme();
    }

    //! @brief virtual copy constructor
    //! @return pointer to new Polynomial
    Polynomial *Polynomial::Clone() const
    {
      return new Polynomial( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Polynomial::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief GetCoefficients returns the coefficients of the polynomial
    //! @return returns vector of coefficients
    const linal::Vector< double> &Polynomial::GetCoefficients() const
    {
      return m_Coefficients;
    }

    //! @brief SetCoefficients set the coefficients of the polynomial
    //! @param COEFFICIENTS the new coefficients to use
    void Polynomial::SetCoefficients( const linal::Vector< double> &COEFFICIENTS)
    {
      m_Coefficients = COEFFICIENTS;
      UpdateScheme();
    }

    //! @brief GetDegree get the highest power of x that this polynomial calculates
    //! @return highest power of x that this polynomial calculates
    size_t Polynomial::GetDegree() const
    {
      return m_Coefficients.IsEmpty() ? 0 : m_Coefficients.GetSize() - 1;
    }

    //! @brief get the polynomial's formula
    //! @return a string like "math::Polynomial(x) = 5.01 + 3.51*x + -0.42*x^2"
    const std::string &Polynomial::GetScheme() const
    {
      return m_Scheme;
    }

    //! @brief Get an iterator to the first coefficient
    //! @return an iterator to the first coefficient
    Polynomial::const_iterator Polynomial::Begin() const
    {
      return m_Coefficients.Begin();
    }

    //! @brief Get an iterator to the end of the coefficients
    //! @return an iterator to the end of the coefficients
    Polynomial::const_iterator Polynomial::End() const
    {
      return m_Coefficients.End();
    }

    //! @brief Get an iterator to the first coefficient
    //! @return an iterator to the first coefficient
    Polynomial::iterator Polynomial::Begin()
    {
      return m_Coefficients.Begin();
    }

    //! @brief Get an iterator to the end of the coefficients
    //! @return an iterator to the end of the coefficients
    Polynomial::iterator Polynomial::End()
    {
      return m_Coefficients.End();
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer Polynomial::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Represents a polynomial of the form f(x) = a0+a1*x+...+an*x^n");
      serializer.AddInitializer
      (
        "coefficients",
        "coefficients of the terms",
        io::Serialization::GetAgent( &m_Coefficients)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculates the derivative of the polynomial
    //! @return the derivative of the polynomial function
    Polynomial Polynomial::Derivative() const
    {
      if( GetDegree() == 0) // if polynomial is constant, derivative is zero, so just return the empty polynomial
      {
        return Polynomial();
      }

      // return a polynomial with 1
      linal::Vector< double> derivative_coefficients( GetDegree());

      // power law: d(an*x^n+an-1*x^(n-1)+...+a0) = n*an*x^(n-1)+(n-1)*an-1*x^(n-2)+...+a1
      for
      (
        size_t degree( 1), max_degree( GetDegree());
        degree <= max_degree;
        ++degree
      )
      {
        // d(m_Coefficients(degree)*x^degree)=degree*m_Coefficients(degree)*x^(degree-1)
        derivative_coefficients( degree - 1) = double( degree) * m_Coefficients( degree);
      }

      return MakeFromCoefficients( derivative_coefficients);
    }

    //! @brief create a polynomial with x ^ DEGREE that is the linear least squares fit to given coordinates
    //! @param X_COORDINATES the x coordinates to fit the polynomial to
    //! @param Y_COORDINATES the y coordinates to fit the polynomial to
    //! @param DEGREE the maximum degree that the polynomial should be of
    //! @return the resulting polynomial
    Polynomial Polynomial::MakeFromCoordinates
    (
      const linal::Vector< double> &X_COORDINATES,
      const linal::Vector< double> &Y_COORDINATES,
      const size_t &DEGREE
    )
    {
      BCL_Assert
      (
        X_COORDINATES.GetSize() == Y_COORDINATES.GetSize(),
        "given vectors of non-equal size"
      );

      const size_t number_coordinates( X_COORDINATES.GetSize());

      if( X_COORDINATES.IsEmpty()) // no coordinates to satisfy
      {
        return Polynomial(); // return a polynomial that always returns 0
      }
      else if( X_COORDINATES.GetSize() <= DEGREE) // under-determined system
      {
        // return a smaller polynomial that will pass through the points
        return MakeFromCoordinates( X_COORDINATES, Y_COORDINATES, number_coordinates - 1);
      }

      // given coordinates X(0...m), compute matrix
      // M =
      // {
      //   1,X(0),X(0)^2,X(0)^3,...
      //   1,X(1),X(1)^2,X(1)^3,...
      //   ........................
      //   1,X(m),X(m)^2,X(m)^3,...
      // }
      linal::Matrix< double> coordinate_powers( number_coordinates, DEGREE + 1);

      for( size_t coordinate_index( 0); coordinate_index < number_coordinates; ++coordinate_index)
      {
        double *const matrix_row( coordinate_powers[ coordinate_index]);

        double coordinate( X_COORDINATES( coordinate_index));
        double coordinate_to_degree_power( 1.0); // coordinate ^ degree
        for( size_t degree( 0); degree <= DEGREE; ++degree, coordinate_to_degree_power *= coordinate)
        {
          matrix_row[ degree] = coordinate_to_degree_power;
        }
      }

      // compute the least-squares approximation to the polynomial running through those points
      return MakeFromCoefficients
             (
               // just take the solutions, ignore the chi squared
               LinearLeastSquares::SolutionAndChiSquared( coordinate_powers, Y_COORDINATES).First()
             );
    }

    //! @brief create a polynomial with COEFFICIENTS
    //! @param COEFFICIENTS the coefficients of the new polynomial
    //! @return the resulting polynomial
    Polynomial Polynomial::MakeFromCoefficients( const linal::Vector< double> &COEFFICIENTS)
    {
      Polynomial poly;
      poly.SetCoefficients( COEFFICIENTS);

      return poly;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator() compute the polynomial using horner's method
    //! @param ARGUMENT x for the polynomial
    //! @return returns the result of evaluating the polynomial at x=ARGUMENT
    double Polynomial::operator ()( const double &ARGUMENT) const
    {
      if( m_Coefficients.IsEmpty()) // handle the trivial case of no coefficients
      {
        return 0.0;
      }

      double result( m_Coefficients.Last()); // start off with the highest-degree coefficient

      // evaluate polynomial using horner's method, which is faster and more numerically stable than computing powers
      // ((an*x+an-1)x+an-2)+...
      for
      (
        // set itr just before the last element, which we've already used
        const double *itr( m_Coefficients.End() - 2),
                     *itr_end( m_Coefficients.Begin() - 1); // reverse end
        itr != itr_end;
        --itr // iterate from highest degree to lowest degree
      )
      {
        result = result * ARGUMENT + *itr;
      }

      return result;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Polynomial::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Coefficients, ISTREAM);

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Polynomial::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Coefficients, OSTREAM, INDENT);

      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief update the scheme
    void Polynomial::UpdateScheme()
    {
      std::ostringstream formula;

      //write scheme
      formula << "F(x) = ";

      //write the formula
      if( !m_Coefficients.IsEmpty())
      {
        formula << m_Coefficients.First(); // constant term

        for( size_t power( 1), degree( GetDegree()); power <= degree; ++power)
        {
          formula << " + " << m_Coefficients( power) << "x"; // successive powers of x

          if( power > 1) // write x^power if the power > 1 (writing x^1 looks sloppy)
          {
            formula << "^" << power;
          }
        }
      }
      else
      {
        formula << "0.0";
      }

      m_Scheme = formula.str();
    }

  } // namespace math
} // namespace bcl

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
#include "math/bcl_math_power_equals.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API PowerEquals< double>;
    template class BCL_API PowerEquals< float>;

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_quaternion.h"

// includes from bcl - sorted alphabetically
#include "random/bcl_random_uniform_distribution.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Quaternion::s_Instance
    (
      GetObjectInstances().AddInstance( new Quaternion())
    );

    Quaternion Conjugate( const Quaternion &QUATERNION)
    {
      Quaternion quat;
      quat[ 0] = QUATERNION[ 0];
      for( size_t i( 1); i < Quaternion::s_Size; ++i)
      {
        quat[ i] = -QUATERNION[ i];
      }
      return quat;
    }

  ///////////////
  // operators //
  ///////////////

    Quaternion &Quaternion::operator *=( const double T)
    {
      std::transform( m_Values, m_Values + s_Size, m_Values, std::bind1st( std::multiplies< double>(), T));

      return *this;
    }

    Quaternion &Quaternion::operator =( const Quaternion &QUATERNION)
    {
      std::copy( QUATERNION.m_Values, QUATERNION.m_Values + s_Size, m_Values);

      return *this;
    }

    Quaternion &Quaternion::operator /=( const double T)
    {
      std::transform( m_Values, m_Values + s_Size, m_Values, std::bind1st( std::multiplies< double>(), double( 1) / T));

      return *this;
    }

  ////////////////
  // operations //
  ////////////////

    // J.J. Kuffner "Effective Sampling and Distance Metrices for 3D Rigid Body Path Planning" 2004 IEEE Int'l Conf. on Robotics and Automation (ICRA 2004)
    //! set random rotation
    Quaternion &Quaternion::SetRand()
    {
      double s( random::GetGlobalRandom().Random< double>( double( 1))),
             sigma1( Sqrt( s - 1)),
             sigma2( Sqrt( s)),
             theta1( 2 * g_Pi * random::GetGlobalRandom().Random< double>( double( 1))),
             theta2( 2 * g_Pi * random::GetGlobalRandom().Random< double>( double( 1)));

      m_Values[ 0] = std::cos( theta2) * sigma2;
      m_Values[ 1] = std::sin( theta1) * sigma1;
      m_Values[ 2] = std::cos( theta1) * sigma1;
      m_Values[ 3] = std::cos( theta2) * sigma2;

      return *this;
    }

    Quaternion &Quaternion::Equal( const double R, const double X, const double Y, const double Z)
    {
      m_Values[ 0] = R;
      m_Values[ 1] = X;
      m_Values[ 2] = Y;
      m_Values[ 3] = Z;

      return *this;
    }

    Quaternion &Quaternion::Equal( const double R, const linal::Vector3D &VECTOR3D)
    {
      m_Values[ 0] = R;
      std::copy( VECTOR3D.Begin(), VECTOR3D.End(), m_Values + 1);

      return *this;
    }

    double Quaternion::Norm() const
    {
      double tmp( 0);
      for( size_t i( 0); i < s_Size; ++i)
      {
        tmp += m_Values[ i] * m_Values[ i];
      }
      return Sqrt( tmp);
    }

  //////////////////////
  // input and output //
  //////////////////////

    std::istream &Quaternion::Read( std::istream &ISTREAM)
    {
      // read all values
      for( double *ptr( m_Values), *ptr_end( m_Values + s_Size); ptr != ptr_end; ++ptr)
      {
        ISTREAM >> *ptr;
      }

      return ISTREAM;
    }

    std::ostream &Quaternion::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // insert indent
      io::Serialize::InsertIndent( OSTREAM, INDENT);

      // write all values
      for( const double *ptr( m_Values), *ptr_end( m_Values + s_Size); ptr != ptr_end; ++ptr)
      {
        OSTREAM << *ptr << '\t';
      }

      // end
      return OSTREAM;
    }

    Quaternion &Quaternion::VectorToQuaternion( const linal::Vector3D &VECTOR3D)
    {
      m_Values[ 0] = 0;
      std::copy( VECTOR3D.Begin(), VECTOR3D.End(), m_Values + 1);

      return *this;
    }

    Quaternion operator *( const Quaternion &QUAT_1, const Quaternion &QUAT_2)
    {
      Quaternion quat;

      quat[ 0] = QUAT_1[ 0] * QUAT_2[ 0] - QUAT_1[ 1] * QUAT_2[ 1] - QUAT_1[ 2] * QUAT_2[ 2] - QUAT_1[ 3] * QUAT_2[ 3];
      quat[ 1] = QUAT_1[ 0] * QUAT_2[ 1] + QUAT_1[ 1] * QUAT_2[ 0] + QUAT_1[ 2] * QUAT_2[ 3] - QUAT_1[ 3] * QUAT_2[ 2];
      quat[ 2] = QUAT_1[ 0] * QUAT_2[ 2] + QUAT_1[ 2] * QUAT_2[ 0] + QUAT_1[ 3] * QUAT_2[ 1] - QUAT_1[ 1] * QUAT_2[ 3];
      quat[ 3] = QUAT_1[ 0] * QUAT_2[ 3] + QUAT_1[ 3] * QUAT_2[ 0] + QUAT_1[ 1] * QUAT_2[ 2] - QUAT_1[ 2] * QUAT_2[ 1];

      // end
      return quat;
    }

    Quaternion operator *( const linal::Vector3D &VECTOR3D, const Quaternion &QUATERNION)
    {
      Quaternion quat;
      quat[ 0] = 0;
      for( size_t i = 0; i < 3; ++i)
      {
        quat[ i + 1] = VECTOR3D( i);
      }

      quat = quat * QUATERNION;

      return quat;
    }

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_range.hpp"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

  //////////
  // data //
  //////////

    const std::string &RangeBorders::GetConditionLeftChars()
    {
      static const std::string s_condition_left_char( "[(");
      return s_condition_left_char;
    }
    const std::string &RangeBorders::GetConditionRightChars()
    {
      static const std::string s_condition_right_char( "])");
      return s_condition_right_char;
    }

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Range< double>;
    template class BCL_API Range< float>;
    template class BCL_API Range< int>;
    template class BCL_API Range< unsigned int>;
    template class BCL_API Range< unsigned long>;
    template class BCL_API Range< unsigned long long>;
    template class BCL_API Range< bool>;
    template class BCL_API Range< char>;

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_range_set.hpp"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API RangeSet< double>;
    template class BCL_API RangeSet< float>;
    template class BCL_API RangeSet< int>;
    template class BCL_API RangeSet< unsigned int>;
    template class BCL_API RangeSet< unsigned long>;
    template class BCL_API RangeSet< unsigned long long>;
    template class BCL_API RangeSet< bool>;
    template class BCL_API RangeSet< char>;

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_roc_curve.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_vector.h"
#include "math/bcl_math_binary_function_bind_second.h"
#include "math/bcl_math_comparisons.h"
#include "math/bcl_math_contingency_matrix_measures.h"
#include "math/bcl_math_cubic_spline_damped.h"
#include "math/bcl_math_linear_function.h"
#include "math/bcl_math_piecewise_function.h"
#include "math/bcl_math_running_average.h"
#include "math/bcl_math_spline_border_type.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_binary_function_stl_wrapper.h"
#include "util/bcl_util_logger_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from members
    //! @param FP # of false positives
    //! @param TP # of true positives
    //! @param CUTOFF actual cutoff value
    ROCCurve::Point::Point( const size_t &FP, const size_t &TP, const double &CUTOFF) :
      m_FalsePositives( FP),
      m_TruePositives( TP),
      m_Cutoff( CUTOFF)
    {
    }

    //! @brief Clone function
    //! @return pointer to new ROCCurve
    ROCCurve::Point *ROCCurve::Point::Clone() const
    {
      return new ROCCurve::Point( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ROCCurve::Point::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get number of true positives
    //! @return number of true positives = TP
    size_t ROCCurve::Point::GetNumberTruePositives() const
    {
      return m_TruePositives;
    }

    //! @brief get number of false positives
    //! @return number of false positives = FP
    size_t ROCCurve::Point::GetNumberFalsePositives() const
    {
      return m_FalsePositives;
    }

    //! @brief number predicted positives
    //! @return sum of TP and FP = P'
    size_t ROCCurve::Point::GetNumberPredictedPositives() const
    {
      return m_TruePositives + m_FalsePositives;
    }

    //! @brief get the cutoff
    //! @return the cutoff
    double ROCCurve::Point::GetCutoff() const
    {
      return m_Cutoff;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Construct a contingency matrix, given the final point in the ROC
    //! @param FINAL_POINT the final point in the ROC
    //! @return a contingency matrix for this object
    ContingencyMatrix ROCCurve::Point::GetContingencyMatrix( const Point &FINAL_POINT) const
    {
      return
        ContingencyMatrix
        (
          m_TruePositives,
          m_FalsePositives,
          FINAL_POINT.m_TruePositives - m_TruePositives,
          FINAL_POINT.m_FalsePositives - m_FalsePositives
        );
    }

    //! @brief add true positives to the point
    //! @param RESULT true for true positive, false for false positive
    //! @param CUTOFF the new cutoff
    void ROCCurve::Point::Update( const bool &RESULT, const double &CUTOFF)
    {
      // update the counts
      RESULT ? ++m_TruePositives : ++m_FalsePositives;
      // update the threshold
      m_Cutoff = CUTOFF;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ROCCurve::Point::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_FalsePositives, ISTREAM);
      io::Serialize::Read( m_TruePositives, ISTREAM);
      io::Serialize::Read( m_Cutoff, ISTREAM);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &ROCCurve::Point::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_FalsePositives, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_TruePositives, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Cutoff, OSTREAM, INDENT);
      return OSTREAM;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a list of pairs of prediction and test results not yet classified by a threshold
    //! @param TEST_RESULTS list of pairs of prediction ( double) and test outputs( double)
    //! @param THRESHOLD the threshold value to be used for classifying the outputs
    //! @param POSITIVES_ABOVE_THRESHOLD flag whether positives are above threshold
    //! ( classification means that results above the threshold are considered positive = 1)
    ROCCurve::ROCCurve
    (
      storage::List< storage::Pair< double, double> > &TEST_RESULTS,
      const double THRESHOLD,
      bool POSITIVES_ABOVE_THRESHOLD
    ) :
      m_SortedCounts()
    {
      // if positives are above threshold
      if( POSITIVES_ABOVE_THRESHOLD)
      {
        // sort in descending order
        TEST_RESULTS.Sort
        (
          storage::PairBinaryPredicateFirst< double, double>
          (
            util::BinaryFunctionSTLWrapper< std::greater< double> >()
          )
        );
      }
      else
      {
        // sort in ascending order
        TEST_RESULTS.Sort
        (
          storage::PairBinaryPredicateFirst< double, double>
          (
            util::BinaryFunctionSTLWrapper< std::less< double> >()
          )
        );
      }

      // initialize roc curve based on POSITIVES_ABOVE_THRESHOLD
      Initialize
      (
        POSITIVES_ABOVE_THRESHOLD
        ? ClassifyResults( TEST_RESULTS, Comparisons< double>::GetEnums().CreateUnaryPredicate( Comparisons< double>::GetEnums().e_Greater, THRESHOLD))
        : ClassifyResults( TEST_RESULTS, Comparisons< double>::GetEnums().CreateUnaryPredicate( Comparisons< double>::GetEnums().e_Less, THRESHOLD))
      );
    }

    //! @brief constructor from a list of pairs of prediction and test results not yet classified by a threshold
    //! @param TEST_RESULTS list of pairs of prediction ( double) and test outputs( double)
    //! @param THRESHOLD the threshold value to be used for classifying the outputs
    //! ( classification means that results above the threshold are considered positive = 1)
    ROCCurve::ROCCurve
    (
      const storage::List< storage::Pair< double, double> > &TEST_RESULTS,
      const double THRESHOLD
    ) :
      m_SortedCounts()
    {
      Initialize( ClassifyResults( TEST_RESULTS, Comparisons< double>::GetEnums().CreateUnaryPredicate( Comparisons< double>::GetEnums().e_Greater, THRESHOLD)));
    }

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ROCCurve::s_Instance
    (
      GetObjectInstances().AddInstance( new ROCCurve())
    );

  ////////////////
  // operations //
  ////////////////

    //! @brief initializes the roc curve with a list of pair of predicted values and expected outputs
    void ROCCurve::Initialize( const storage::List< storage::Pair< double, bool> > &RESULTS_CLASSIFIED)
    {
      // reset the counts list
      m_SortedCounts.Reset();

      // if there are no results, return
      if( RESULTS_CLASSIFIED.IsEmpty())
      {
        return;
      }

      // get an iterator range on the results
      storage::List< storage::Pair< double, bool> >::const_iterator result_itr( RESULTS_CLASSIFIED.Begin());

      // add the first result to the counts of true/false
      Point counts;
      counts.Update( result_itr->Second(), result_itr->First());
      ++result_itr;

      // insert the current counts
      m_SortedCounts.PushBack( counts);
      m_SortedCounts.AllocateMemory( RESULTS_CLASSIFIED.GetSize());

      // iterate over sorted results
      for
      (
        storage::List< storage::Pair< double, bool> >::const_iterator
          result_itr_prev( RESULTS_CLASSIFIED.Begin()),
          result_itr_end( RESULTS_CLASSIFIED.End());
        result_itr != result_itr_end;
        ++result_itr, ++result_itr_prev
      )
      {
        // check for new threshold value
        if( result_itr_prev->First() != result_itr->First())
        {
          // insert the current counts
          m_SortedCounts.PushBack( m_SortedCounts.LastElement());
        }

        // initialize the count that corresponds to the classification of the result
        m_SortedCounts.LastElement().Update( result_itr->Second(), result_itr->First());
      }
    }

    //! @brief returns the total number of false positives
    //! @return total number of false positives
    const size_t ROCCurve::GetNumberActualNegatives() const
    {
      BCL_Assert( !m_SortedCounts.IsEmpty(), "There are no results stored in this ROCCurve curve!");
      return m_SortedCounts.LastElement().GetNumberFalsePositives();
    }

    //! @brief returns the total number of true positives
    //! @return total number of true positives
    const size_t ROCCurve::GetNumberActualPositives() const
    {
      BCL_Assert( !m_SortedCounts.IsEmpty(), "There are no results stored in this ROCCurve curve!");
      return m_SortedCounts.LastElement().GetNumberTruePositives();
    }

    //! @brief calculates and returns the area under the curve
    //! @return curve integral
    double ROCCurve::Integral() const
    {
      return Integral( m_SortedCounts.End());
    }

    //! @brief calculates and returns the weighted area under the curve by weighting function
    //! @param WEIGHTING_FUNCTION FunctionInterface which weights the ROC curve and allows Integral calculation
    //! @return entire weighted curve integral
    double ROCCurve::Integral( const FunctionInterfaceSerializable< double, double> &WEIGHTING_FUNCTION, const bool ROC_X_AXIS_LOG10_SCALING) const
    {
      return Integral( m_SortedCounts.End(), WEIGHTING_FUNCTION, ROC_X_AXIS_LOG10_SCALING);
    }

    //! @brief calculates and returns the area under the curve for the given top fraction ( between 0 and 1)
    //! @param FRACTION top fraction ( between 0 and 1) to be considered when calculating integral
    //! @return curve integral
    double ROCCurve::Integral( const double FRACTION) const
    {
      return Integral( GetEndIteratorForInterval( FRACTION), LinearFunction( 0, 1));
    }

    //! @brief calculates and returns the area under the curve for the given top fraction ( between 0 and 1)
    //!        weighted by a given function
    //! @param FRACTION top fraction ( between 0 and 1) to be considered when calculating integral
    //! @param WEIGHTING_FUNCTION FunctionInterface which weights the ROC curve and allows Integral calculation
    //!        accordingly
    //! @return weighted curve integral
    double ROCCurve::Integral
    (
      const double FRACTION,
      const FunctionInterfaceSerializable< double, double> &WEIGHTING_FUNCTION,
      const bool ROC_X_AXIS_LOG10_SCALING
    ) const
    {
      return Integral( GetEndIteratorForInterval( FRACTION), WEIGHTING_FUNCTION, ROC_X_AXIS_LOG10_SCALING);
    }

    //! @brief calculates and returns the area under the curve for the given top fraction ( between 0 and 1),
    //!        the contingency matrix functions to plot on x- and y-axis can be specified
    //! @param ROC_X_AXIS function pointer to a member function in math::ContingencyMatrix
    //! @param ROC_Y_AXIS function pointer to a member function in math::ContingencyMatrix
    //! @param ROC_X_AXIS_LOG10_SCALING enables log10 scaling of x axis
    //! @return curve integral
    double ROCCurve::Integral
    (
      const function::MemberConst< ContingencyMatrix, double> &ROC_X_AXIS,
      const function::MemberConst< ContingencyMatrix, double> &ROC_Y_AXIS,
      const bool ROC_X_AXIS_LOG10_SCALING
    ) const
    {
      return Integral
      (
        m_SortedCounts.End(),
        LinearFunction( 0, 1),
        ROC_X_AXIS,
        ROC_Y_AXIS,
        Range< double>( 0, 1),
        ROC_X_AXIS_LOG10_SCALING
      );
    }

    //! @brief calculates and returns the area under the curve for the given top fraction ( between 0 and 1),
    //!        the contingency matrix functions to plot on x- and y-axis can be specified
    //! @param FRACTION top fraction ( between 0 and 1) to be considered when calculating integral
    //! @param ROC_X_AXIS function pointer to a member function in math::ContingencyMatrix
    //! @param ROC_Y_AXIS function pointer to a member function in math::ContingencyMatrix
    //! @param ROC_X_AXIS_LOG10_SCALING enables log10 scaling of x axis
    //! @return curve integral
    double ROCCurve::Integral
    (
      const double FRACTION,
      const function::MemberConst< ContingencyMatrix, double> &ROC_X_AXIS,
      const function::MemberConst< ContingencyMatrix, double> &ROC_Y_AXIS,
      const bool ROC_X_AXIS_LOG10_SCALING
    ) const
    {
      return Integral
      (
        GetEndIteratorForInterval( FRACTION),
        LinearFunction( 0, 1),
        ROC_X_AXIS,
        ROC_Y_AXIS,
        Range< double>( 0, 1),
        ROC_X_AXIS_LOG10_SCALING
      );
    }

    //! @brief calculates and returns the area under the curve for the given top fraction ( between 0 and 1),
    //!        the contingency matrix functions to plot on x- and y-axis can be specified
    //! @param ROC_X_AXIS function pointer to a member function in ContingencyMatrix
    //! @param ROC_Y_AXIS function pointer to a member function in ContingencyMatrix
    //! @param ROC_X_AXIS_FRACTION_CUTOFF x axis value cutoff rather then specifying a number of counts
    //! @param ROC_X_AXIS_LOG10_SCALING enables log10 scaling of x axis
    //! @return curve integral
    double ROCCurve::Integral
    (
      const function::MemberConst< ContingencyMatrix, double> &ROC_X_AXIS,
      const function::MemberConst< ContingencyMatrix, double> &ROC_Y_AXIS,
      const Range< double> &ROC_X_AXIS_FRACTION_CUTOFF,
      const bool ROC_X_AXIS_LOG10_SCALING
    ) const
    {
      return Integral
      (
        m_SortedCounts.End(),
        LinearFunction( 0, 1),
        ROC_X_AXIS,
        ROC_Y_AXIS,
        ROC_X_AXIS_FRACTION_CUTOFF,
        ROC_X_AXIS_LOG10_SCALING
      );
    }

    //! @brief calculates and returns the area under the curve for the given top fraction ( between 0 and 1)
    //!        weighted by a given function
    //!        the contingency matrix functions to plot on x- and y-axis can be specified
    //! @param FRACTION top fraction ( between 0 and 1) to be considered when calculating integral
    //! @param WEIGHTING_FUNCTION FunctionInterface which weights the ROC curve and allows Integral calculation
    //!        accordingly
    //! @param ROC_X_AXIS function pointer to a member function in math::ContingencyMatrix
    //! @param ROC_Y_AXIS function pointer to a member function in math::ContingencyMatrix
    //! @param ROC_X_AXIS_LOG10_SCALING enables log10 scaling of x axis
    //! @return weighted curve integral
    double ROCCurve::Integral
    (
      const double FRACTION,
      const FunctionInterfaceSerializable< double, double> &WEIGHTING_FUNCTION,
      const function::MemberConst< ContingencyMatrix, double> &ROC_X_AXIS,
      const function::MemberConst< ContingencyMatrix, double> &ROC_Y_AXIS,
      const bool ROC_X_AXIS_LOG10_SCALING
    ) const
    {
      return Integral
      (
        GetEndIteratorForInterval( FRACTION),
        WEIGHTING_FUNCTION,
        ROC_X_AXIS,
        ROC_Y_AXIS,
        Range< double>( 0, 1),
        ROC_X_AXIS_LOG10_SCALING
      );
    }

    //! @brief calculates and returns the area under the curve for the given top fraction ( between 0 and 1)
    //!        weighted by a given function
    //!        the contingency matrix functions to plot on x- and y-axis can be specified
    //! @param WEIGHTING_FUNCTION FunctionInterface which weights the ROC curve and allows Integral calculation
    //!        accordingly
    //! @param ROC_X_AXIS function pointer to a member function in ContingencyMatrix
    //! @param ROC_Y_AXIS function pointer to a member function in ContingencyMatrix
    //! @param ROC_X_AXIS_FRACTION_CUTOFF x axis cutoff range rather then specifying a number of counts
    //! @param ROC_X_AXIS_LOG10_SCALING enables log10 scaling of x axis
    //! @return weighted curve integral
    double ROCCurve::Integral
    (
      const FunctionInterfaceSerializable< double, double> &WEIGHTING_FUNCTION,
      const function::MemberConst< ContingencyMatrix, double> &ROC_X_AXIS,
      const function::MemberConst< ContingencyMatrix, double> &ROC_Y_AXIS,
      const Range< double> &ROC_X_AXIS_FRACTION_CUTOFF,
      const bool ROC_X_AXIS_LOG10_SCALING
    ) const
    {
      return Integral
      (
        m_SortedCounts.End(),
        WEIGHTING_FUNCTION,
        ROC_X_AXIS,
        ROC_Y_AXIS,
        ROC_X_AXIS_FRACTION_CUTOFF,
        ROC_X_AXIS_LOG10_SCALING
      );
    }

    //! @brief calculates ContingencyMatrix with true/false positive/negatives
    //! @param FRACTION consider the first fraction as being predicted positive, the remainder predicted negative
    ContingencyMatrix ROCCurve::ContingencyMatrixFraction( const double FRACTION) const
    {
      if( m_SortedCounts.IsEmpty())
      {
        return ContingencyMatrix();
      }

      // closed interval [0,1] for valid fraction
      static const Range< double> s_valid_fraction_range( 0.0, 1.0);

      // check that fraction is between 0 and 1
      BCL_Assert
      (
        s_valid_fraction_range.IsWithin( FRACTION),
        "given fraction not in interval " + s_valid_fraction_range.GetString() + ": " + util::Format()( FRACTION)
      );

      // compute the # of predictions desired
      const size_t desired_predictions( std::min( size_t( std::ceil( GetNumberResults() * FRACTION)), GetNumberResults()));

      // number of predictions in first fraction
      storage::Vector< Point>::const_iterator
        itr_fraction_last( m_SortedCounts.Begin()), itr_fraction_end( m_SortedCounts.Last());

      // walk in the list until the # of predicted values is correct
      while
      (
        itr_fraction_last != itr_fraction_end
        && itr_fraction_last->GetNumberPredictedPositives() < desired_predictions
      )
      {
        ++itr_fraction_last;
      }

      // if the total counts is exactly equal to the desired predictions, return the desired contingency matrix
      if( itr_fraction_last->GetNumberPredictedPositives() == desired_predictions)
      {
        return itr_fraction_last->GetContingencyMatrix( m_SortedCounts.LastElement());
      }

      // get the counts at the itr, which is the minimal # of counts greater than or equal to the desired # of predictions
      Point counts_gte_desired( *itr_fraction_last);

      // get the counts just before itr, which is the minimal # of counts less than the desired # of predictions
      Point counts_lt_desired;

      if( itr_fraction_last != m_SortedCounts.Begin())
      {
        // go to the previous value;
        counts_lt_desired = *--itr_fraction_last;
      }

      const size_t total_lt( counts_lt_desired.GetNumberPredictedPositives());
      const size_t total_gte( counts_gte_desired.GetNumberPredictedPositives());

      // compute true and false positives slopes
      const double true_positives_slope
      (
        double( counts_gte_desired.GetNumberTruePositives() - counts_lt_desired.GetNumberTruePositives())
        / double( total_gte - total_lt)
      );
      const double false_positives_slope
      (
        double( counts_gte_desired.GetNumberFalsePositives() - counts_lt_desired.GetNumberFalsePositives())
        / double( total_gte - total_lt)
      );

      // compute true positives and false positives using interpolation between the two values
      const size_t true_positives
      (
        counts_lt_desired.GetNumberTruePositives() + true_positives_slope * ( desired_predictions - total_lt)
      );
      const size_t false_positives
      (
        counts_lt_desired.GetNumberFalsePositives() + false_positives_slope * ( desired_predictions - total_lt)
      );

      // end
      return
        Point( false_positives, true_positives, counts_lt_desired.GetCutoff()).GetContingencyMatrix
        (
          m_SortedCounts.LastElement()
        );
    }

    //! @brief calculates ContingencyMatrix with true/false positive/negatives
    //! @param FRACTION consider the first fraction as being predicted positive, the remainder predicted negative
    double ROCCurve::CutoffFraction( const double FRACTION) const
    {
      if( m_SortedCounts.IsEmpty())
      {
        return 0.0;
      }

      /////////////
      if( FRACTION == 0.0)
      {
        return m_SortedCounts.FirstElement().GetCutoff();
      }

      if( FRACTION == 1.0)
      {
        return m_SortedCounts.LastElement().GetCutoff();
      }
      ////////////////

      // closed interval [0,1] for valid fraction
      static const Range< double> s_valid_fraction_range( 0.0, 1.0);

      // check that fraction is between 0 and 1
      BCL_Assert
      (
        s_valid_fraction_range.IsWithin( FRACTION),
        "given fraction not in interval " + s_valid_fraction_range.GetString() + ": " + util::Format()( FRACTION)
      );

      // compute the # of predictions desired
//      const size_t desired_predictions( std::min( size_t( std::ceil( GetNumberResults() * FRACTION)), GetNumberResults()));
      const float desired_predictions( float( GetNumberResults()) * FRACTION);

      // number of predictions in first fraction
      storage::Vector< Point>::const_iterator
        itr_fraction_last( m_SortedCounts.Begin()), itr_fraction_end( m_SortedCounts.Last());

      // walk in the list until the # of predicted values is correct
      while
      (
        itr_fraction_last != itr_fraction_end
        && float( itr_fraction_last->GetNumberPredictedPositives()) < desired_predictions
      )
      {
        ++itr_fraction_last;
      }

      // if the total counts is exactly equal to the desired predictions, return the desired contingency matrix
      if( itr_fraction_last->GetNumberPredictedPositives() == desired_predictions)
      {
        return itr_fraction_last->GetCutoff();
      }

      // get the counts at the itr, which is the minimal # of counts greater than or equal to the desired # of predictions
      Point counts_gte_desired( *itr_fraction_last);

      // get the counts just before itr, which is the minimal # of counts less than the desired # of predictions
      Point counts_lt_desired;

      if( itr_fraction_last != m_SortedCounts.Begin())
      {
        // go to the previous value;
        counts_lt_desired = *--itr_fraction_last;
      }
      else
      {
          return itr_fraction_last->GetCutoff();
      }
//
//
      // {
        //if( ++itr_fraction_last != m_SortedCounts.End())
       // {
         // counts_lt_desired = *itr_fraction_last;
       // }
      //}

      const float total_lt( float( counts_lt_desired.GetNumberPredictedPositives()));
      const float total_gte( float( counts_gte_desired.GetNumberPredictedPositives()));

      // compute true and false positives slopes
      const double cutoff_slope
      (
        double( counts_gte_desired.GetCutoff() - counts_lt_desired.GetCutoff())
        / double( total_gte - total_lt)
      );
      return cutoff_slope * ( desired_predictions - total_lt) + counts_lt_desired.GetCutoff();
    }

    //! @brief calculates and returns the area under the curve for the counts between provided regions
    //! @param END iterator to the end of the interval for which the integral is going to be calculated
    //! @return curve integral for values between specified regions
    double ROCCurve::Integral
    (
      const storage::Vector< Point>::const_iterator &END
    ) const
    {
      // return integral with no weighting by a function
      return Integral( END, LinearFunction( 0, 1));
    }

    //! @brief calculates and returns the area under the curve plotting False Positive Rate on the x axis
    //!        and True Positive Rate on the y axis for the counts between provided regions
    //!        weighted with a mathematical function
    //! @param END iterator to the end of the interval for which the integral is going to be calculated
    //! @param WEIGHTING_FUNCTION FunctionInterface which weights the ROC curve and allows Integral calculation
    //!        accordingly
    //! @return weighted curve integral for values between specified regions
    double ROCCurve::Integral
    (
      const storage::Vector< Point>::const_iterator &END,
      const FunctionInterfaceSerializable< double, double> &WEIGHTING_FUNCTION,
      const bool ROC_X_AXIS_LOG10_SCALING
    ) const
    {
      return Integral
      (
        END,
        WEIGHTING_FUNCTION,
        &ContingencyMatrix::GetFalsePositiveRate, // X-coordinate
        &ContingencyMatrix::GetTruePositiveRate,  // Y-coordinate
        Range< double>( 0, 1),
        ROC_X_AXIS_LOG10_SCALING
      );
    }

    //! @brief calculates and returns the area under the curve for the counts between provided regions
    //!        weighted with a mathematical function
    //! @param END iterator to the end of the interval for which the integral is going to be calculated
    //! @param WEIGHTING_FUNCTION FunctionInterface which weights the ROC curve and allows Integral calculation
    //!        accordingly
    //! @param ROC_X_AXIS function pointer to a member function in ContingencyMatrix
    //! @param ROC_Y_AXIS function pointer to a member function in ContingencyMatrix
    //! @param ROC_X_AXIS_RANGE range on x axis that specifies a begin and end fraction
    //!        rather then specifying a number of counts through param END
    //! @param ROC_X_AXIS_LOG10_SCALING enables log10 scaling of x axis
    //! @return weighted curve integral for values between specified regions,
    //!         the integral is normalized according to the specified range depending on the x axis scaling
    double ROCCurve::Integral
    (
      const storage::Vector< Point>::const_iterator &END,
      const FunctionInterfaceSerializable< double, double> &WEIGHTING_FUNCTION,
      const function::MemberConst< ContingencyMatrix, double> &ROC_X_AXIS,
      const function::MemberConst< ContingencyMatrix, double> &ROC_Y_AXIS,
      const Range< double> ROC_X_AXIS_RANGE,
      const bool ROC_X_AXIS_LOG10_SCALING
    ) const
    {
      BCL_Assert
      (
        ROC_X_AXIS_RANGE.GetMax() > 0 && ROC_X_AXIS_RANGE.GetMin() <= 1,
        "x axis cutoff has to be between 0 < x_axis_cutoff <= 1! current range is " + util::Format()( ROC_X_AXIS_RANGE)
      );

      storage::Vector< Point>::const_iterator counts_itr( m_SortedCounts.Begin());

      // if interval is empty return 0
      if( counts_itr == END)
      {
        BCL_MessageStd
        (
          "This interval asked for ROCCurve curve is empty, therefore returning 0 for integral"
        );
        return double( 0);
      }

      // total count of data points
      const size_t total_count( GetNumberResults() - 1);

      // initialize minimum of roc curve x-axis dependent on enabled log scaling
      const double roc_min
      (
        ROC_X_AXIS_LOG10_SCALING
        ? ROC_X_AXIS_RANGE.GetMin() == 0
          ? std::log10( 1 / double( total_count))
          : std::log10( ROC_X_AXIS_RANGE.GetMin())
        : ROC_X_AXIS_RANGE.GetMin()
      );

      // x axis step size = 1 / (N + P)
      const double step_size
      (
        ROC_X_AXIS_LOG10_SCALING
        ? ( std::log10( ROC_X_AXIS_RANGE.GetMax()) - roc_min) / double( total_count)
        : ( 1.0 - roc_min) / double( total_count)
      );

      // runnning average of integral segments
      RunningAverage< double> integral;

      // map for associating x-axis values to y-axis values
      storage::Map< double, RunningAverage< double> > xy_plot( ComputeXAxisYAxisMapping( ROC_X_AXIS, ROC_Y_AXIS));

      // get the actual end of the x-axis value based on the END iterator provided
      float end_fraction
      (
        END != m_SortedCounts.End()
        ? ROC_X_AXIS( END->GetContingencyMatrix( m_SortedCounts.LastElement()))
        : (
              ROC_X_AXIS_LOG10_SCALING
            ? Pow( 10.0, total_count * step_size + roc_min)
            : total_count * step_size + roc_min
          )
      );
      if( end_fraction > ROC_X_AXIS_RANGE.GetMax())
      {
        end_fraction = ROC_X_AXIS_RANGE.GetMax();
      }

      BCL_Assert( xy_plot.GetSize(), "Roc curve: x axis value map is empty!")

      if( xy_plot.GetSize() == size_t( 1))
      {
        return xy_plot.Begin()->second;
      }

      // iterate over every count and calculate the integral of every integration segment
      for( size_t distance_counter( 0); distance_counter < total_count; ++distance_counter)
      {
        // determine fraction of x axis
        double fraction( distance_counter * step_size);

        // if log scaling is activated
        if( ROC_X_AXIS_LOG10_SCALING)
        {
          fraction = Pow( 10.0, fraction + roc_min);
        }
        else
        {
          fraction += roc_min;
        }

        // stop if x axis range is reached
        if( fraction > end_fraction)
        {
          break;
        }
        else if( fraction < ROC_X_AXIS_RANGE.GetMin())
        {
          continue;
        }

        storage::Map< double, RunningAverage< double> >::const_iterator itr_xy_plot_lower, itr_xy_plot_upper;

        // determine boundaries for every integration segment
        itr_xy_plot_lower = itr_xy_plot_upper = xy_plot.UpperBound( fraction);

        if( itr_xy_plot_upper == xy_plot.End())
        {
          itr_xy_plot_lower = --itr_xy_plot_upper;
        }
        if( itr_xy_plot_lower != xy_plot.Begin())
        {
          --itr_xy_plot_lower;
        }

        if( itr_xy_plot_lower == itr_xy_plot_upper)
        {
          integral += itr_xy_plot_lower->second * WEIGHTING_FUNCTION( fraction);
          continue;
        }

        // final boundary values for iteration segment
        const double upper_x( itr_xy_plot_upper->first);
        const double upper_y( itr_xy_plot_upper->second);
        const double lower_x( itr_xy_plot_lower->first);
        const double lower_y( itr_xy_plot_lower->second);

        double integral_delta( 0);

        // final boundary values when log scaling is enabled
        const double slope( ( upper_y - lower_y) / ( upper_x - lower_x));
        const double y_icept( lower_y - lower_x * slope);

        integral_delta = fraction * slope + y_icept;

        // averaged integral
        integral += integral_delta * WEIGHTING_FUNCTION( fraction);
      }

      // if integral is zero, write out a message
      if( integral == 0)
      {
        BCL_MessageStd( "calculated integral has value 0 !");
        return integral;
      }

      BCL_MessageDbg( "calculated integral: " + util::Format()( integral));
      return integral;
    }

    //! @brief compute the x axis and y axis mapping for computing an interpolated integral under the curve
    //! @param ROC_X_AXIS function pointer to a member function in ContingencyMatrix
    //! @param ROC_Y_AXIS function pointer to a member function in ContingencyMatrix
    //! @return x axis and y axis mapping
    storage::Map< double, RunningAverage< double> > ROCCurve::ComputeXAxisYAxisMapping
    (
      const function::MemberConst< ContingencyMatrix, double> &ROC_X_AXIS,
      const function::MemberConst< ContingencyMatrix, double> &ROC_Y_AXIS
    ) const
    {
      storage::Map< double, RunningAverage< double> > xy_plot;

      // add data for the 0 point
      {
        ContingencyMatrix con_matrix
        (
          0,                          // true_positives
          0,                          // false_positives
          GetNumberActualPositives(), // false_negatives
          GetNumberActualNegatives()  // true_negatives
        );

        // only add the value if it is defined on both axis!
        const double x_axis_value( ROC_X_AXIS( con_matrix));
        const double y_axis_value( ROC_Y_AXIS( con_matrix));
        if( util::IsDefined( x_axis_value) && util::IsDefined( y_axis_value))
        {
          xy_plot[ x_axis_value] += y_axis_value;
        }
      }

      if( !m_SortedCounts.IsEmpty())
      {
        Point last_point( m_SortedCounts.LastElement());
        // iterate over rates
        for
        (
          storage::Vector< Point>::const_iterator
            counts_itr( m_SortedCounts.Begin()), counts_itr_end( m_SortedCounts.End());
          counts_itr != counts_itr_end;
          ++counts_itr
        )
        {
          ContingencyMatrix con_matrix( counts_itr->GetContingencyMatrix( last_point));
          xy_plot[ ROC_X_AXIS( con_matrix)] += ROC_Y_AXIS( con_matrix);
        }
      }

      return xy_plot;
    }

    //! @brief Convert into a map from threshold to contingency matrix
    storage::Map< double, ContingencyMatrix> ROCCurve::ToMap() const
    {
      storage::Map< double, ContingencyMatrix> matrices;

      Point last_point( m_SortedCounts.LastElement());
      // iterate over rates
      for
      (
        storage::Vector< Point>::const_iterator
          counts_itr( m_SortedCounts.Begin()), counts_itr_end( m_SortedCounts.End());
        counts_itr != counts_itr_end;
        ++counts_itr
      )
      {
        matrices[ counts_itr->GetCutoff()] = counts_itr->GetContingencyMatrix( last_point);
      }
      return matrices;
    }

    //! @brief get the optima for a particular contingency matrix measure
    //! @param MEASURE the measure to optimize
    //! @return an iterator to the maximum point for the particular measure and the maximum value
    std::pair< storage::Vector< ROCCurve::Point>::const_iterator, double> ROCCurve::GetMaxima
    (
      const util::FunctionInterfaceSerializable< ContingencyMatrix, double> &MEASURE
    ) const
    {
      Point last_point( m_SortedCounts.LastElement());
      storage::Vector< ROCCurve::Point>::const_iterator itr_best( m_SortedCounts.Begin());
      double best_val( util::GetUndefined< double>());
      // iterate over points
      for
      (
        storage::Vector< Point>::const_iterator
          counts_itr( m_SortedCounts.Begin()), counts_itr_end( m_SortedCounts.End());
        counts_itr != counts_itr_end;
        ++counts_itr
      )
      {
        const double local_metric_value( MEASURE( counts_itr->GetContingencyMatrix( last_point)));
        if
        (
          util::IsDefined( local_metric_value)
          && ( !util::IsDefined( best_val) || local_metric_value > best_val)
        )
        {
          best_val = local_metric_value;
          itr_best = counts_itr;
        }
      }
      return std::make_pair( itr_best, best_val);
    }

    //! @brief creates a thinned version of the ca1ounts so that only every Nth( PERIODICITY * n) element is kept
    //! @brief This function is to be used for plotting.
    //! @param PERIODICITY periodicity with which elements should be selected
    //! @return thinned ROC curve counts
    ROCCurve ROCCurve::GetThinnedRocCurvePeriodicity( const size_t PERIODICITY) const
    {
      // check that fraction is between 1 and total size
      BCL_Assert
      (
        PERIODICITY != 0 && PERIODICITY < GetNumberResults(),
        "The provided periodicity: " + util::Format()( PERIODICITY) +
        " cannot be 0 or larger than total number of results: " + util::Format()( GetNumberResults())
      );

      // initialize the new counts list
      storage::Vector< Point> new_counts;

      // initialize iterator to the beginning
      storage::Vector< Point>::const_iterator iterator( m_SortedCounts.Begin());

      // initialize end iterator
      const storage::Vector< Point>::const_iterator end_iterator( m_SortedCounts.End());

      // start the target # of predictions counter
      size_t desired_predictions( PERIODICITY);

      // iterate until you hit the end
      while( iterator != end_iterator)
      {
        // push back this element to new counts
        new_counts.PushBack( *iterator);

        // walk in the list until the # of predicted values is correct
        while( iterator != m_SortedCounts.End() && iterator->GetNumberPredictedPositives() < desired_predictions)
        {
          ++iterator;
        }
        desired_predictions += PERIODICITY;
      }

      // create a new ROCCurve to hold the counts
      ROCCurve roc_with_new_counts;

      // set the counts
      roc_with_new_counts.m_SortedCounts = new_counts;

      // return the roc curve with the new counts
      return roc_with_new_counts;
    }

    //! @brief creates a thinned version of the counts so that only TOTAL_NUMBER_ENTRIES elements are kept
    //! @brief This function is to be used for plotting.
    //! @param TOTAL_NUMBER_ENTRIES periodicity with which elements should be selected
    //! @return thinned ROC curve counts
    ROCCurve ROCCurve::GetThinnedRocCurveTotal( const size_t TOTAL_NUMBER_ENTRIES) const
    {
      return GetThinnedRocCurvePeriodicity( std::max( size_t( 1), m_SortedCounts.GetSize() / TOTAL_NUMBER_ENTRIES));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream from a formatted output
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ROCCurve::ReadPlottingTable( std::istream &ISTREAM)
    {
      // reset the counts
      m_SortedCounts.Reset();

      // read the identifier to make sure the file has correct header
      std::string identifier;
      ISTREAM >> identifier;
      BCL_Assert( identifier == "FALSE_POSITIVE", "The provided identifier is not correct" + identifier);

      // read the remaining of the header
      ISTREAM >> identifier;

      // read the number of examples
      size_t number_examples, example_ctr( 0);
      ISTREAM >> number_examples;

      // while number of examples is not reached
      while( example_ctr < number_examples && !ISTREAM.eof())
      {
        // read a new count pair and insert it into m_SortedCounts
        size_t first, second;
        ISTREAM >> first >> second;
        m_SortedCounts.PushBack( Point( first, second, double( example_ctr)));

        // increase the example ctr
        ++example_ctr;
      }

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream in a formatted output style the result between BEGIN and END
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &ROCCurve::WritePlottingTable
    (
      std::ostream &OSTREAM
    ) const
    {
      // print number of elements
      OSTREAM << "FALSE_POSITIVE\tTRUE_POSITIVE\t" << m_SortedCounts.GetSize() << '\n';

      // x axis step size = 1 / (N + P)
      const double step_size( 1 / double( GetNumberResults()));

      // total count of data points
      const double max_fraction( 1);

      ContingencyMatrix con_matrix;

      // iterate over counts
      for( double fraction( step_size); fraction <= max_fraction; fraction += step_size)
      {
        con_matrix = ContingencyMatrixFraction( fraction);

        // output
        OSTREAM << con_matrix.GetNumberFalsePositives() << '\t' << con_matrix.GetNumberTruePositives() << '\n';
      }

      // end
      return OSTREAM;
    }

    //! @brief write to std::ostream in a formatted output style the result between BEGIN and END
    //! @param OSTREAM output stream to write to
    //! @param FORMAT format the rate
    //! @return output stream which was written to
    std::ostream &ROCCurve::WriteRatePlottingTable
    (
      std::ostream &OSTREAM,
      const util::Format &FORMAT
    ) const
    {
      // print number of elements
      OSTREAM << "FALSE_POSITIVE_RATE\tTRUE_POSITIVE_RATE\t"
              << GetNumberActualNegatives() << '\t' << GetNumberActualPositives()
              << '\n';

      // end
      return WriteRatePlottingTableGeneric
      (
        OSTREAM,
        &ContingencyMatrix::GetFalsePositiveRate, // X-coordinate
        &ContingencyMatrix::GetTruePositiveRate,  // Y-coordinate
        Range< double>( 0.0, 1.0),                // X-axis range
        FORMAT
      );
    }

    //! @brief write to std::ostream in a formatted output style the result between BEGIN and END
    //! @param OSTREAM output stream to write to
    //! @param ROC_X_AXIS function pointer to a member function in ContingencyMatrix, specifies what is plotted
    //!        on x axis of roc curve
    //! @param ROC_Y_AXIS function pointer to a member function in ContingencyMatrix, specifies what is plotted
    //!        on y axis of roc curve
    //! @param ROC_X_AXIS_RANGE x axis value range specifying the x axis value when to begin and stop plotting
    //!        default is set to [ 0.0, 1.0], 0% to 100%
    //! @param FORMAT format the rate
    //! @param PLOTTING_TABLE_TITLE title that describes columns for x axis values and y axis values
    //! @return output stream which was written to
    std::ostream &ROCCurve::WriteRatePlottingTableGeneric
    (
      std::ostream &OSTREAM,
      const function::MemberConst< ContingencyMatrix, double> &ROC_X_AXIS,
      const function::MemberConst< ContingencyMatrix, double> &ROC_Y_AXIS,
      const Range< double> &ROC_X_AXIS_RANGE,
      const util::Format &FORMAT,
      const std::string &PLOTTING_TABLE_TITLE
    ) const
    {
      BCL_Assert
      (
        ROC_X_AXIS_RANGE.GetMin() >= 0 && ROC_X_AXIS_RANGE.GetMax() <= 1,
        "WriteRatePlottingTableGeneric: x axis cutoff range has to be between [ 0.0, 1.0] but is " + util::Format()( ROC_X_AXIS_RANGE)
      );

      if( !PLOTTING_TABLE_TITLE.empty())
      {
        // print title
        OSTREAM << PLOTTING_TABLE_TITLE << '\n';
      }

      // iterate over rates
      for
      (
        storage::Vector< Point>::const_iterator
          counts_itr( m_SortedCounts.Begin()), counts_itr_end( m_SortedCounts.End());
        counts_itr != counts_itr_end;
        ++counts_itr
      )
      {
        // update running contingency matrix at current rate
        ContingencyMatrix con_matrix( counts_itr->GetContingencyMatrix( m_SortedCounts.LastElement()));

        const double roc_x_axis( ROC_X_AXIS( con_matrix));

        // if x axis cutoff did not enter allowed range yet then continue integration
        if( roc_x_axis < ROC_X_AXIS_RANGE.GetMin())
        {
          continue;
        }

        // if x axis cutoff exceeding allowed maximum range then break integration
        if( roc_x_axis > ROC_X_AXIS_RANGE.GetMax())
        {
          break;
        }

        // output
        OSTREAM << FORMAT( roc_x_axis) << '\t' << FORMAT( ROC_Y_AXIS( con_matrix)) << '\n';
      }

      // end
      return OSTREAM;
    }

    //! @brief write all common contingency matrix measures to the file, return a vector of strings for each measure that was written
    //! @param DATA_FILENAME file name for the table
    //! @param MEASURES a list of measures of interest
    //! @return strings of all measures that were written
    storage::Vector< std::string> ROCCurve::WriteRatePlottingTableComplete
    (
      const std::string &DATA_FILENAME,
      const storage::Vector< std::string> &MEASURES
    ) const
    {
      io::OFStream table;
      io::File::MustOpenOFStream( table, DATA_FILENAME);
      // create all measures
      storage::Map< std::string, ContingencyMatrixMeasures> measure_map;

      // cutoff and local-PPV are special (since they are really unknown by the contingency matrix),
      // so get its alias and check whether it will be used
      size_t cutoff_index( util::GetUndefined< size_t>());
      size_t local_ppv_index( util::GetUndefined< size_t>());
      for( size_t measure( 0); measure < ContingencyMatrixMeasures::s_NumberMeasures; ++measure)
      {
        ContingencyMatrixMeasures this_measure( static_cast< ContingencyMatrixMeasures::Measure>( measure));
        measure_map[ this_measure.GetAlias()] = this_measure;
      }
      storage::Vector< std::string> measure_strings( MEASURES.IsEmpty() ? ContingencyMatrixMeasures::s_NumberMeasures : MEASURES.GetSize());
      storage::Vector< ContingencyMatrixMeasures> all_measures( MEASURES.IsEmpty() ? ContingencyMatrixMeasures::s_NumberMeasures : MEASURES.GetSize());
      if( MEASURES.IsEmpty())
      {
        for( size_t measure( 0); measure < ContingencyMatrixMeasures::s_NumberMeasures; ++measure)
        {
          all_measures( measure) = ContingencyMatrixMeasures( static_cast< ContingencyMatrixMeasures::Measure>( measure));
          measure_strings( measure) = all_measures( measure).GetAlias();
          table << measure_strings( measure) << '\t';
        }
        cutoff_index = size_t( ContingencyMatrixMeasures::e_Cutoff);
        local_ppv_index = size_t( ContingencyMatrixMeasures::e_LocalPPV);
      }
      else
      {
        measure_strings = MEASURES;
        for( size_t measure( 0), n_measures( MEASURES.GetSize()); measure < n_measures; ++measure)
        {
          all_measures( measure) = measure_map[ MEASURES( measure)];
          if( all_measures( measure).GetMeasure() == ContingencyMatrixMeasures::e_Cutoff)
          {
            cutoff_index = measure;
          }
          else if( all_measures( measure).GetMeasure() == ContingencyMatrixMeasures::e_LocalPPV)
          {
            local_ppv_index = measure;
          }
          table << measure_strings( measure) << '\t';
        }
      }
      table << '\n';

      PiecewiseFunction local_ppv;
      if( util::IsDefined( local_ppv_index))
      {
        local_ppv = GetLocalPPVCurve();
      }

      // iterate over rates
      util::Format format;
      for
      (
        storage::Vector< Point>::const_iterator
          counts_itr( m_SortedCounts.Begin()), counts_itr_end( m_SortedCounts.End());
        counts_itr != counts_itr_end;
        ++counts_itr
      )
      {
        // update running contingency matrix at current rate
        ContingencyMatrix con_matrix( counts_itr->GetContingencyMatrix( m_SortedCounts.LastElement()));

        size_t measure_id( 0);
        for
        (
          storage::Vector< ContingencyMatrixMeasures>::const_iterator
            itr_measure( all_measures.Begin()), itr_measure_end( all_measures.End());
          itr_measure != itr_measure_end;
          ++itr_measure, ++measure_id
        )
        {
          if( measure_id != cutoff_index && measure_id != local_ppv_index)
          {
            table << format( ( *itr_measure)( con_matrix)) << '\t';
          }
          else if( measure_id == cutoff_index)
          {
            table << format( counts_itr->GetCutoff()) << '\t';
          }
          else if( measure_id == local_ppv_index)
          {
            table << local_ppv( counts_itr->GetCutoff()) << '\t';
          }
        }
        table << '\n';
      }

      io::File::CloseClearFStream( table);

      // end
      return measure_strings;
    }

    //! @brief get the localized PPV curve
    //! The localized PPV-curve yields a more precise estimate of the local PPVs for a given range
    //! The standard ROC-curve answers the question: what is the PPV of all points predicted at or above a given
    //! cutoff (assuming +parity). The localized ppv-curve, by contrast, provides the PPV of values near a particular
    //! cutoff. This curve is produced by iteratively finding the optimal PPV interest, then removing all
    //! entries (P/N) at or above/below (depending on parity) from the ROC curve. PPV is essentially the only metric
    //! (other than the obscure false discovery rate = 1-PPV) that is independent of negatives and the total number of
    //! positives, and so is the only metric that can be computed locally in this manner.
    //! @return a piecewise function for the PPV curve, which is guaranteed to be monotonic
    PiecewiseFunction ROCCurve::GetLocalPPVCurve() const
    {
      ROCCurve copy( *this);
      const ContingencyMatrixMeasures ppv( ContingencyMatrixMeasures::e_PositivePredictiveValue);

      // true if positive cases should be predicted higher than negative cases
      const bool positives_higher( m_SortedCounts.FirstElement().GetCutoff() > m_SortedCounts.LastElement().GetCutoff());

      // keep track of all x, y points for construction of spline
      storage::Vector< double> xs, ys;

      // keep finding the point with the highest PPV value; then remove all predicted positives contributing to that ppv
      // after recording it in the xs (cutoff) and ys (ppv) vector.
      while( !copy.m_SortedCounts.IsEmpty())
      {
        std::pair< storage::Vector< ROCCurve::Point>::const_iterator, double> optima( copy.GetMaxima( ppv));
        storage::Vector< ROCCurve::Point>::const_iterator itr_optima( optima.first);
        ROCCurve::Point optimal_point( *itr_optima);

        // add this optima to the splined points
        xs.PushBack( copy.CutoffFraction( optimal_point.GetNumberPredictedPositives() / ( 2.0 * copy.GetNumberResults())));

//       if(xs.GetSize() >= size_t(2) && xs.LastElement() < xs(xs.GetSize()-2) && copy.m_SortedCounts.GetSize() <= size_t(5))
//       {
//        BCL_MessageStd(util::Format()(copy.m_SortedCounts));
//       }

        ys.PushBack( optima.second);

        // remove the elements that form this optimal PPV from the curve
        copy.m_SortedCounts.RemoveElements( 0, std::distance( copy.GetSortedCounts().Begin(), itr_optima) + 1);
        // subtract the predicted positives from the copy's points
        for
        (
          storage::Vector< ROCCurve::Point>::iterator
            itr_copy( copy.m_SortedCounts.Begin()), itr_copy_end( copy.m_SortedCounts.End());
          itr_copy != itr_copy_end;
          ++itr_copy
        )
        {
          *itr_copy =
            ROCCurve::Point
            (
              itr_copy->GetNumberFalsePositives() - optimal_point.GetNumberFalsePositives(),
              itr_copy->GetNumberTruePositives() - optimal_point.GetNumberTruePositives(),
              itr_copy->GetCutoff()
            );
        }
      }

      // spline requires sorted points; but if positives are higher, then x is a descending vector
      if( positives_higher)
      {
        std::reverse( xs.Begin(), xs.End());
        std::reverse( ys.Begin(), ys.End());
      }

      // create a damped spline to ensure monotonicity

      // no begin/end slope.
      // Assume that the best ppv on this data is really the best PPV.
      // In general there's no way of knowing what the actual maximum that the model will predict is or what the ppv is
      // at that value, without having access to training dataset predictions
      CubicSplineDamped csvd;
      csvd.Train
      (
        linal::Vector< double>( xs.Begin(), xs.End()),
        linal::Vector< double>( ys.Begin(), ys.End()),
        0.0,
        0.0
      );

      // create a trivial piecewise function. This allows for easy editing later on if we decide to use a different
      // function
      storage::List
      <
        storage::Pair< Range< double>, util::ShPtr< FunctionInterfaceSerializable< double, double> > >
      > functions;
      functions.PushBack
      (
        storage::Pair< Range< double>, util::ShPtr< FunctionInterfaceSerializable< double, double> > >
        (
          Range< double>
          (
            RangeBorders::e_LeftClosed,
            -std::numeric_limits< double>::max(),
            std::numeric_limits< double>::max(),
            RangeBorders::e_RightClosed
          ),
          util::ShPtr< FunctionInterfaceSerializable< double, double> >( csvd.Clone())
        )
      );

      return PiecewiseFunction( functions);
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ROCCurve::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SortedCounts, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &ROCCurve::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SortedCounts, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief returns iterator on the the end of the specified fraction of the list
    //! @param FRACTION top fraction ( between 0 and 1) to be used in consequent operations
    //! @return iterator to ending of the requested interval
    storage::Vector< ROCCurve::Point>::const_iterator
    ROCCurve::GetEndIteratorForInterval( const double FRACTION) const
    {
      // assert that the fraction is between 0 and including 1
      BCL_Assert
      (
        FRACTION > 0.0 && FRACTION <= 1.0,
        "The provided fraction should be between 0 and including 1 and not" + util::Format()( FRACTION)
      );

      // initialize an iterator
      storage::Vector< Point>::const_iterator iterator( m_SortedCounts.Begin());

      // compute the # of predictions desired
      const size_t desired_predictions( size_t( GetNumberResults() * FRACTION));

      // walk in the list until the # of predicted values is correct
      while
      (
        iterator != m_SortedCounts.End() &&
        iterator->GetNumberPredictedPositives() < desired_predictions
      )
      {
        ++iterator;
      }

      // return
      return iterator;
    }

    //! @brief static function to classify provided results and return the classified expected outputs
    //! @param UNCLASSIFIED_RESULTS results that are not yet classified
    //! @param UNARY_CLASSIFIER classifier to be used for classifying expected outputs
    //! @return classified results
    storage::List< storage::Pair< double, bool> > ROCCurve::ClassifyResults
    (
      const storage::List< storage::Pair< double, double> > &UNCLASSIFIED_RESULTS,
      const FunctionInterfaceSerializable< double, bool> &UNARY_CLASSIFIER
    )
    {
      // initialize the list to be returned
      storage::List< storage::Pair< double, bool> > classified_results;

      // iterate over provided results
      for
      (
        storage::List< storage::Pair< double, double> >::const_iterator
          result_itr( UNCLASSIFIED_RESULTS.Begin()), result_itr_end( UNCLASSIFIED_RESULTS.End());
        result_itr != result_itr_end; ++result_itr
      )
      {
        if( util::IsDefined( result_itr->Second()))
        {
          // push back the classified result
          classified_results.PushBack
          (
            storage::Pair< double, bool>( result_itr->First(), UNARY_CLASSIFIER( result_itr->Second()))
          );
        }
      }

      // end
      return classified_results;
    }

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_rotation_matrix_2d.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"
#include "random/bcl_random_uniform_distribution.h"

// external includes - sorted alphabetically
#include <math.h>

namespace bcl
{
  namespace math
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> RotationMatrix2D::s_Instance
    (
      GetObjectInstances().AddInstance( new RotationMatrix2D())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! construct empty RotationMatrix2D
    RotationMatrix2D::RotationMatrix2D()
    {
      // set the matrix up as the unitary matrix (no rotation)
      m_RotationMatrix2D( 0, 0) = m_RotationMatrix2D( 1, 1) = 1.0;
    }

    //! construct RotationMatrix2D from RotationAngle
    RotationMatrix2D::RotationMatrix2D( const double &ANGLE)
    {
      SetAngle( ANGLE);
    }

    //! construct RotationMatrix2D from a 2x2 matrix
    RotationMatrix2D::RotationMatrix2D( const linal::Matrix2x2< double> &ROTATION) :
      m_RotationMatrix2D( ROTATION)
    {
    }

    //! virtual copy constructor
    RotationMatrix2D *RotationMatrix2D::Clone() const
    {
      return new RotationMatrix2D( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RotationMatrix2D::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! return the matrix
    const linal::Matrix2x2< double> &RotationMatrix2D::GetMatrix() const
    {
      return m_RotationMatrix2D;
    }

    //! @brief return the angle that this matrix rotates by
    //! @return rotation angle that this matrix represents (in radians)
    double RotationMatrix2D::GetAngle() const
    {
      return acos( *m_RotationMatrix2D.Begin());
    }

    //! @brief set the angle that this matrix rotates by
    //! @param THETA rotation angle that this matrix represents (in radians)
    void RotationMatrix2D::SetAngle( const double &THETA)
    {
      SetFromCosSin( cos( THETA), sin( THETA));
    }

    //! @brief set given cos(theta) and sin(theta)
    //! @param COS cosine of theta
    //! @param SIN sine theta
    void RotationMatrix2D::SetFromCosSin( const double &COS, const double &SIN)
    {
      m_RotationMatrix2D( 0, 0) = m_RotationMatrix2D( 1, 1) = COS;
      m_RotationMatrix2D( 1, 0) = -( m_RotationMatrix2D( 0, 1) = -SIN);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief set random rotation from 0 - MAX_ROTATIONANGLE
    //! @param MAX_ROTATIONANGLE the maximum rotation angle
    RotationMatrix2D &RotationMatrix2D::SetRand( const double MAX_ROTATIONANGLE)
    {
      //random angle between 0 and 2pi
      SetAngle( MAX_ROTATIONANGLE * random::GetGlobalRandom().Double());
      return *this;
    }

    //! @brief set the matrix up to be a givens rotation matrix; which is a numerically stable solution used in tridiagonalization
    //! @param A, B the vector to rotate in the givens rotation matrix problem:
    //! | c -s || A | = | r |
    //! | s  c || B |   | 0 |
    //! @return R (equal to pythag(A,B))
    double RotationMatrix2D::MakeGivens( const double &A, const double &B)
    {
      const double abs_a( Absolute( A)), abs_b( Absolute( B));
      if( abs_a > abs_b)
      {
        if( B == 0.0)
        {
          SetFromCosSin( A > 0.0 ? 1.0 : -1.0, 0.0);
          return abs_a;
        }
        // |A| > |B|, neither B nor A are 0
        const double t( B / A);
        const double u( A > 0.0 ? std::sqrt( 1.0 + t * t) : -std::sqrt( 1.0 + t * t));
        const double c( 1.0 / u);
        SetFromCosSin( c, -c * t);
        return A * u;
      }
      else if( A == 0.0)
      {
        BCL_Assert( A != B, "Givens rotation undefined if A and B are both 0!");
        SetFromCosSin( 0.0, B > 0.0 ? -1.0 : 1.0);
        return abs_b;
      }

      // |A| <= |B|, neither B nor A are 0
      const double t( A / B);
      const double u( B > 0.0 ? std::sqrt( 1.0 + t * t) : -std::sqrt( 1.0 + t * t));
      const double s( -1.0 / u);
      SetFromCosSin( -s * t, s);
      return B * u;
    }

  ///////////////
  // operators //
  ///////////////

    //! operator *= RotationMatrix2D
    RotationMatrix2D &RotationMatrix2D::operator *=( const RotationMatrix2D &MATRIX)
    {
      m_RotationMatrix2D *= ( MATRIX.m_RotationMatrix2D);
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write RotationMatrix2D to std::ostream using the given util::Format
    std::ostream &RotationMatrix2D::Write( std::ostream &OSTREAM, const size_t INDENT, const util::Format &FORMAT) const
    {
      io::Serialize::Write( m_RotationMatrix2D, OSTREAM, INDENT, FORMAT) << '\n';
      return OSTREAM;
    }

    //! write RotationMatrix2D to std::ostream using the given util::Format
    std::ostream &RotationMatrix2D::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_RotationMatrix2D, OSTREAM, INDENT) << '\n';
      return OSTREAM;
    }

    //! read RotationMatrix2D from io::IFStream
    std::istream &RotationMatrix2D::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_RotationMatrix2D, ISTREAM);
      return ISTREAM;
    }

  } // namespace math
} // namespace bcl

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
#include "math/bcl_math_rotation_matrix_3d.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "random/bcl_random_uniform_distribution.h"
#include "util/bcl_util_logger_interface.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> RotationMatrix3D::s_Instance
    (
      GetObjectInstances().AddInstance( new RotationMatrix3D())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! construct empty RotationMatrix3D
    RotationMatrix3D::RotationMatrix3D()
    {
      // set the matrix up as the unitary matrix (no rotation)
      m_RotationMatrix3D( 0, 0) = m_RotationMatrix3D( 1, 1) = m_RotationMatrix3D( 2, 2) = 1.0;
    }

    //! construct RotationMatrix3D from RotationAxis and RotationAngle
    RotationMatrix3D::RotationMatrix3D( const coord::Axis &AXIS, const double &ANGLE)
    {
      Init( AXIS, ANGLE);
    }

    //! construct RotationMatrix3D from any linal::Vector3D around which shall be rotated and the Angle.
    RotationMatrix3D::RotationMatrix3D( const linal::Vector3D &AXIS, const double ANGLE)
    {
      if( !ANGLE)
      {
        // angle is 0, ignore the axis (which does not matter in this case)
        // set the matrix up as the unitary matrix (no rotation)
        m_RotationMatrix3D( 0, 0) = m_RotationMatrix3D( 1, 1) = m_RotationMatrix3D( 2, 2) = 1.0;
      }
      else
      {
        // normalize the provided axis
        linal::Vector3D norm_axis( AXIS);
        norm_axis.Normalize();

        const double rcos( cos( ANGLE));

        // minor improvements in performance and numeric stability can be had by precomputing
        // norm(AXIS) * (1-cos(ANGLE)) and
        // norm(AXIS) * sin(ANGLE)
        linal::Vector3D norm_axis_remain( norm_axis);
        norm_axis_remain *= double( 1.0) - rcos;

        linal::Vector3D norm_axis_rsin( norm_axis);
        norm_axis_rsin *= sin( ANGLE);

        //first col
        m_RotationMatrix3D( 0, 0) =                rcos + norm_axis.X() * norm_axis_remain.X();
        m_RotationMatrix3D( 1, 0) =  norm_axis_rsin.Z() + norm_axis.Y() * norm_axis_remain.X();
        m_RotationMatrix3D( 2, 0) = -norm_axis_rsin.Y() + norm_axis.Z() * norm_axis_remain.X();

        //second col
        m_RotationMatrix3D( 0, 1) = -norm_axis_rsin.Z() + norm_axis.X() * norm_axis_remain.Y();
        m_RotationMatrix3D( 1, 1) =                rcos + norm_axis.Y() * norm_axis_remain.Y();
        m_RotationMatrix3D( 2, 1) =  norm_axis_rsin.X() + norm_axis.Z() * norm_axis_remain.Y();

        //third col
        m_RotationMatrix3D( 0, 2) =  norm_axis_rsin.Y() + norm_axis.X() * norm_axis_remain.Z();
        m_RotationMatrix3D( 1, 2) = -norm_axis_rsin.X() + norm_axis.Y() * norm_axis_remain.Z();
        m_RotationMatrix3D( 2, 2) =                rcos + norm_axis.Z() * norm_axis_remain.Z();
      }
    }

    //! construct RotationMatrix3D from three RotationAxis and three RotationAngles
    RotationMatrix3D::RotationMatrix3D( const coord::Axis AXES[ 3], const double ANGLE[ 3])
    {
      Init( AXES[ 0], ANGLE[ 0]);
      RotationMatrix3D b( AXES[ 1], ANGLE[ 1]), c( AXES[ 2], ANGLE[ 2]);
      operator *= ( b);
      operator *= ( c);
    }

    // assume rotation of PHI around x-axis, THETA around y-axis and PSI around z-axis
    //! construct RotationMatrix3D from three Euler angles
    RotationMatrix3D::RotationMatrix3D( const double PHI, const double THETA, const double PSI)
    {
      // pre calculate some values for efficiency
      const double c_x( cos( PHI));         // A
      const double s_x( sin( PHI));         // B
      const double c_y( cos( THETA));       // C
      const double s_y( sin( THETA));       // D
      const double c_z( cos( PSI));         // E
      const double s_z( sin( PSI));         // F
      const double c_x_s_y( c_x * s_y);     // AD
      const double s_x_s_y( s_x * s_y);     // BD

      m_RotationMatrix3D( 0, 0) =  c_y * c_z ;
      m_RotationMatrix3D( 0, 1) = -c_y * s_z;
      m_RotationMatrix3D( 0, 2) =  s_y;

      m_RotationMatrix3D( 1, 0) =  s_x_s_y * c_z + c_x * s_z;
      m_RotationMatrix3D( 1, 1) = -s_x_s_y * s_z + c_x * c_z;
      m_RotationMatrix3D( 1, 2) = -s_x * c_y;

      m_RotationMatrix3D( 2, 0) = -c_x_s_y * c_z + s_x * s_z;
      m_RotationMatrix3D( 2, 1) =  c_x_s_y * s_z + s_x * c_z;
      m_RotationMatrix3D( 2, 2) =  c_x * c_y;
    }

    //! construct RotationMatrix3D from a 3x3 matrix
    RotationMatrix3D::RotationMatrix3D( const linal::Matrix3x3< double> &ROTATION) :
      m_RotationMatrix3D( ROTATION)
    {
    }

    //! virtual copy constructor
    RotationMatrix3D *RotationMatrix3D::Clone() const
    {
      return new RotationMatrix3D( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RotationMatrix3D::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! return the matrix
    const linal::Matrix3x3< double> &RotationMatrix3D::GetMatrix() const
    {
      return m_RotationMatrix3D;
    }

    //! @brief return the X Y or Z axis in respect to their reference coordinate system x, y and z (euclidean unit vectors)
    //! @param AXIS axis that will be returned
    //! @return axis in coordinates of reference coordinate system
    linal::Vector3D RotationMatrix3D::GetAxis( const coord::Axis &AXIS) const
    {
      // axis is stored in the columns
      return linal::Vector3D( m_RotationMatrix3D( 0, AXIS), m_RotationMatrix3D( 1, AXIS), m_RotationMatrix3D( 2, AXIS));
    }

    //! @brief return the axis of rotation; result is a unit vector
    //! @return axis of rotation; if this rotation matrix is applied to it, the vector will be unchanged
    linal::Vector3D RotationMatrix3D::GetRotationAxis() const
    {
      linal::Vector3D ret
      (
        m_RotationMatrix3D( 2, 1) - m_RotationMatrix3D( 1, 2),
        m_RotationMatrix3D( 0, 2) - m_RotationMatrix3D( 2, 0),
        m_RotationMatrix3D( 1, 0) - m_RotationMatrix3D( 0, 1)
      );
      const double norm( ret.Norm());
      if( norm == 0.0)
      {
        ret += 1.0 / Sqrt( 3.0);
      }
      else
      {
        ret /= norm;
      }
      return ret;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief return the euler angles for that rotation matrix
    //! http://www.gamedev.net/reference/articles/article1691.asp
    //! @return linal::Vector3D with 1st to 3rd element being alpha, beta, gamma - rotations around x y and z
    linal::Vector3D RotationMatrix3D::EulerAnglesXYZ() const
    {
      linal::Vector3D euler_angles_x_y_z( 0.0, 0.0, 0.0);

      euler_angles_x_y_z( 1) = asin( m_RotationMatrix3D( 0,2));
      const double cos_b     = cos( euler_angles_x_y_z( 1))   ;
      // check for no gimbal lock
      if( Absolute( cos_b) > 0.005)
      {
        // x axis angle
        euler_angles_x_y_z( 0) = atan2( -m_RotationMatrix3D( 1, 2) / cos_b, m_RotationMatrix3D( 2, 2) / cos_b);
        // z axis angle
        euler_angles_x_y_z( 2) = atan2( -m_RotationMatrix3D( 0, 1) / cos_b, m_RotationMatrix3D( 0, 0) / cos_b);
      }
      // there was a gimbal lock
      else
      {
        // x axis angle
        euler_angles_x_y_z( 0) = 0;
        // z axis angle
        euler_angles_x_y_z( 2) = atan2( m_RotationMatrix3D( 1, 0), m_RotationMatrix3D( 1, 1));
      }

      // add 2 pi to all negative angles
      for( double *ptr( euler_angles_x_y_z.Begin()), *ptr_end( euler_angles_x_y_z.End()); ptr != ptr_end; ++ptr)
      {
        if( *ptr < double( 0.0))
        {
          *ptr += 2 * g_Pi;
        }
      }

      // end
      return euler_angles_x_y_z;
    }

    //! @brief return the Euler angles according to z-x'-z'' convention
    //! @return Euler angles according to z-x'-z'' convention
    storage::VectorND< 3, double> RotationMatrix3D::EulerAnglesZXZ() const
    {
      const double beta( std::acos( m_RotationMatrix3D( 2, 2)));
      const double gamma
      (
        std::acos
        (
          m_RotationMatrix3D( 2, 0) / beta == 0.0
          ? std::numeric_limits< double>::epsilon()
          : std::sin( beta)
        )
      );
      const double alpha
      (
        std::asin
        (
          m_RotationMatrix3D( 0, 2) / beta == 0.0
          ? std::numeric_limits< double>::epsilon()
          : std::sin( beta)
        )
      );
      const storage::VectorND< 3, double> euler_angles( alpha, beta, gamma);

      return euler_angles;
    }

    //http://journals.iucr.org/j/issues/2002/05/00/vi0166/index.html
    //! calculate effective rotation angle
    double RotationMatrix3D::EffectiveRotationAngle() const
    {
      //make sure that initial argument for std::acos is between [-1,1]
      return std::acos( std::max( std::min( ( m_RotationMatrix3D.Trace() - 1.0) / 2.0, 1.0), -1.0));
    }

    // J.J. Kuffner "Effective Sampling and Distance Metrices for 3D Rigid Body Path Planning" 2004 IEEE Int'l Conf. on Robotics and Automation (ICRA 2004)
    //! set random rotation from 0 - Pi
    RotationMatrix3D &RotationMatrix3D::SetRand( const double MAX_ROTATIONANGLE)
    {
      do
      {
        //random angle between -pi and pi
        const double a1( MAX_ROTATIONANGLE * ( 2.0 * random::GetGlobalRandom().Double() - 1));
        const double a3( MAX_ROTATIONANGLE * ( 2.0 * random::GetGlobalRandom().Double() - 1));
        //random cos angle between -pi/2 and 3pi/2
        double a2
        (
          std::acos
          (
            random::GetGlobalRandom().Random< double>( double( -1.0), double( 1.0))
          ) + g_Pi * 0.5
        );
        if( random::GetGlobalRandom().Boolean())
        {
          a2 += ( a2 < g_Pi ? g_Pi : -g_Pi);
        }
        a2 *= MAX_ROTATIONANGLE / g_Pi;

        ( *this = RotationMatrix3D( a1, a2, a3));
      } while( this->EffectiveRotationAngle() > MAX_ROTATIONANGLE);
      return *this;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief construct from alpha, beta and gamma angles according to z-x-z convention
    //! @param ALPHA rotation angle around LCS z-axis
    //! @param BETA rotation angle around LCS x'-axis
    //! @param GAMMA rotation angle around LCS z''-axis
    RotationMatrix3D RotationMatrix3D::CreateZXZ( const double ALPHA, const double BETA, const double GAMMA)
    {
      // 3x3 matrix with 0.0 as elements
      linal::Matrix3x3< double> zxz_rotation;

      zxz_rotation( 0, 0) =  cos( ALPHA) * cos( GAMMA) - sin( ALPHA) * cos( BETA) * sin( GAMMA);
      zxz_rotation( 0, 1) = -cos( ALPHA) * sin( GAMMA) - sin( ALPHA) * cos( BETA) * cos( GAMMA);
      zxz_rotation( 0, 2) =  sin(  BETA) * sin( ALPHA);

      zxz_rotation( 1, 0) =  sin( ALPHA) * cos( GAMMA) + cos( ALPHA) * cos( BETA) * sin( GAMMA);
      zxz_rotation( 1, 1) = -sin( ALPHA) * sin( GAMMA) + cos( ALPHA) * cos( BETA) * cos( GAMMA);
      zxz_rotation( 1, 2) = -sin(  BETA) * cos( ALPHA);

      zxz_rotation( 2, 0) =  sin( BETA) * sin( GAMMA);
      zxz_rotation( 2, 1) =  sin( BETA) * cos( GAMMA);
      zxz_rotation( 2, 2) =  cos( BETA);

      return RotationMatrix3D( zxz_rotation);
    }

    //! private helper function to set rotation matrix data from RotationAxis and RotationAngles
    void RotationMatrix3D::Init( const coord::Axis &AXIS, const double &ANGLE)
    {
      m_RotationMatrix3D( AXIS, AXIS) = 1;
      m_RotationMatrix3D( ( AXIS + 1) % 3, ( AXIS + 1) % 3) =  cos( ANGLE);
      m_RotationMatrix3D( ( AXIS + 2) % 3, ( AXIS + 2) % 3) =  cos( ANGLE);
      m_RotationMatrix3D( ( AXIS + 1) % 3, ( AXIS + 2) % 3) = -sin( ANGLE);
      m_RotationMatrix3D( ( AXIS + 2) % 3, ( AXIS + 1) % 3) =  sin( ANGLE);
    }

  ///////////////
  // operators //
  ///////////////

    //! operator *= RotationMatrix3D
    RotationMatrix3D &RotationMatrix3D::operator *=( const RotationMatrix3D &MATRIX)
    {
      m_RotationMatrix3D *= ( MATRIX.m_RotationMatrix3D);
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write RotationMatrix3D to std::ostream using the given util::Format
    std::ostream &RotationMatrix3D::Write( std::ostream &OSTREAM, const size_t INDENT, const util::Format &FORMAT) const
    {
      io::Serialize::Write( m_RotationMatrix3D, OSTREAM, INDENT, FORMAT) << '\n';
      return OSTREAM;
    }

    //! write RotationMatrix3D to std::ostream using the given util::Format
    std::ostream &RotationMatrix3D::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_RotationMatrix3D, OSTREAM, INDENT) << '\n';
      return OSTREAM;
    }

    //! read RotationMatrix3D from io::IFStream
    std::istream &RotationMatrix3D::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_RotationMatrix3D, ISTREAM);
      return ISTREAM;
    }

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_template_instantiations.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    template class BCL_API FunctionInterfaceSerializable< linal::Vector< double>, double>;

    template class BCL_API FunctionInterfaceSerializable< assemble::ProteinModel, double>;

    template class BCL_API BinaryFunctionInterface< assemble::SSEGeometryInterface, assemble::SSE, double>;

    template class BCL_API FunctionInterfaceSerializable< double, double>;

    template class BCL_API FunctionInterfaceSerializable< storage::VectorND< 2, linal::Vector< double> >, storage::VectorND< 2, linal::Vector< double> > >;

    template class BCL_API MutateCombine< assemble::SSEPool>;

    template class BCL_API MutateDecisionNode< assemble::ProteinModel>;

    template class BCL_API MutateDecisionNode< assemble::SSEPool>;

    template class BCL_API MutateInterface< assemble::Domain>;

    template class BCL_API MutateInterface< assemble::ProteinModel>;

    template class BCL_API MutateInterface< TransformationMatrix3D>;

    template class BCL_API SumFunctionMixin< score::ProteinModel>;

    template class BCL_API BinarySumFunction< assemble::SSEPool, biol::Membrane, double, double>;

    template class BCL_API MutateMoveWrapper< assemble::SSE>;

    template class BCL_API MutatePerturbation< assemble::ProteinModel>;

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_tensor.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    template class BCL_API Tensor< double>;

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_times_equals.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API TimesEquals< double>;
    template class BCL_API TimesEquals< float>;

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_transformation_matrix_3d.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"
#include "linal/bcl_linal_matrix_inversion_gauss_jordan.h"
#include "linal/bcl_linal_matrix_operations.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_rotation_matrix_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> TransformationMatrix3D::s_Instance
    (
      GetObjectInstances().AddInstance( new TransformationMatrix3D())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! construct TransformationMatrix3D and set to unitmatrix (empty constructor)
    TransformationMatrix3D::TransformationMatrix3D() :
      m_TransformationMatrix3D( 4, 4, 0.0)
    {
      SetUnitFromZero();
    }

    //! copy constructor
    TransformationMatrix3D::TransformationMatrix3D( const TransformationMatrix3D &MATRIX) :
      m_TransformationMatrix3D( MATRIX.m_TransformationMatrix3D)
    {
    }

    //! @brief move constructor
    TransformationMatrix3D::TransformationMatrix3D( TransformationMatrix3D && MATRIX) :
      m_TransformationMatrix3D( std::move( MATRIX.m_TransformationMatrix3D))
    {
    }

    //! construct TransformationMatrix3D from 4linal::Vector3D ( x, y and z axis and the origin)
    TransformationMatrix3D::TransformationMatrix3D
    (
      const linal::Vector3D &XAXIS,
      const linal::Vector3D &YAXIS,
      const linal::Vector3D &ZAXIS,
      const linal::Vector3D &ORIGIN
    ) :
      m_TransformationMatrix3D( 4, 4)
    {
      double *trans( m_TransformationMatrix3D.Begin());
      for( const double *x_ptr = XAXIS.Begin(); x_ptr != XAXIS.End(); ++x_ptr)
      {
        ( *( trans++)) = ( *x_ptr);
      }

      trans++;
      for( const double *y_ptr = YAXIS.Begin(); y_ptr != YAXIS.End(); ++y_ptr)
      {
        ( *( trans++)) = ( *y_ptr);
      }

      trans++;
      for( const double *z_ptr = ZAXIS.Begin(); z_ptr != ZAXIS.End(); ++z_ptr)
      {
        ( *( trans++)) = ( *z_ptr);
      }

      trans++;
      for( const double *o_ptr = ORIGIN.Begin(); o_ptr != ORIGIN.End(); ++o_ptr)
      {
        ( *( trans++)) = ( *o_ptr);
      }

      //corner is 1
      ( *trans) = 1;
    }

    //! construct TransformationMatrix3D from MatrixND
    TransformationMatrix3D::TransformationMatrix3D( const RotationMatrix3D &MATRIX) :
      m_TransformationMatrix3D( 4, 4, 0.0)
    {
      SetUnitFromZero();
      SetRotation( MATRIX);
    }

    //! construct TransformationMatrix3D from Matrix< double> 4x4
    TransformationMatrix3D::TransformationMatrix3D( const linal::Matrix< double> &MATRIX) :
      m_TransformationMatrix3D( MATRIX)
    {
      BCL_Assert( MATRIX.IsSquare() && MATRIX.GetNumberRows() == size_t( 4), "Argument has to be 4x4 matrix");
    }

    //! @brief construct TransformationMatrix3D from translation
    //! @param TRANSLATION the translation to perform
    TransformationMatrix3D::TransformationMatrix3D( const linal::Vector3D &TRANSLATION) :
      m_TransformationMatrix3D( 4, 4, 0.0)
    {
      SetUnitFromZero();
      m_TransformationMatrix3D.ReplaceRow( 3, TRANSLATION);
    }

    //! construct TransformationMatrix3D from translation
    TransformationMatrix3D::TransformationMatrix3D( const double &X, const double &Y, const double &Z) :
      m_TransformationMatrix3D( 4, 4, 0.0)
    {
      SetUnitFromZero();
      m_TransformationMatrix3D( 3, 0) = X;
      m_TransformationMatrix3D( 3, 1) = Y;
      m_TransformationMatrix3D( 3, 2) = Z;
    }

    //! construct TransformationMatrix3D from coord::Axis and RotationAngle
    TransformationMatrix3D::TransformationMatrix3D( const coord::Axis &AXIS, const double &ANGLE) :
      m_TransformationMatrix3D( 4, 4, 0.0)
    {
      SetUnitFromZero();
      RotationMatrix3D a( AXIS, ANGLE);
      SetRotation( a.m_RotationMatrix3D);
    }

    //! construct TransformationMatrix3D from three coord::Axis and three RotationAngles
    TransformationMatrix3D::TransformationMatrix3D( const coord::Axis AXES[ 3], const double ANGLE[ 3]) :
      m_TransformationMatrix3D( 4, 4, 0.0)
    {
      SetUnitFromZero();
      SetRotation( RotationMatrix3D( AXES, ANGLE));
    }

    //! construct from vector with six components - three rotations and three translations
    //! @param VECTOR elements 1-3 rotation x, y, z; elements 4-6 translations x, y, z
    TransformationMatrix3D::TransformationMatrix3D( const linal::Vector< double> &VECTOR) :
      m_TransformationMatrix3D( 4, 4, 0.0)
    {
      BCL_Assert( VECTOR.GetSize() == 6, "argument vector needs to have 6 elements");
      SetUnitFromZero();
      SetRotation( RotationMatrix3D( &coord::GetAxes().e_X, VECTOR.Begin()));
      std::copy( VECTOR.Begin() + 3, VECTOR.End(), m_TransformationMatrix3D.Begin() + 12);
    }

    //! @brief constructor from a definition state
    //! @param DEFINITION_STATE defined or undefined
    TransformationMatrix3D::TransformationMatrix3D( const util::UndefinedObject DEFINITION_STATE) :
      m_TransformationMatrix3D( 4, 4, 0.0)
    {
      SetUnitFromZero();
      // set the scaling factor to undefined
      m_TransformationMatrix3D( 3, 3) = util::GetUndefinedDouble();
    }

    //! @brief constructor taking the member variable parameters
    //! @param ROTATION_AXIS the axis around which the rotation will take place
    //! @param ROTATION_ORIGIN the origin point around which the rotation will occur
    //! @param ROTATION the amount of rotation
    TransformationMatrix3D::TransformationMatrix3D
    (
      const coord::LineSegment3D &ROTATION_AXIS, const linal::Vector3D &ROTATION_ORIGIN, const double ROTATION
    ) :
      m_TransformationMatrix3D( 4, 4, 0.0)
    {
      SetUnit();
      m_TransformationMatrix3D.ReplaceRow( 3, -ROTATION_ORIGIN);

      // add the desired rotation according to "m_RotationAxis" and "m_Rotation" to "transform"
      this->operator ()( RotationMatrix3D( ROTATION_AXIS.GetDirection(), ROTATION));

      // add the translation of moving back to the original position
      this->operator ()( ROTATION_ORIGIN);
    }

    //! @brief constructor taking two lines, computes the transformation A->B
    //! @param A the starting line
    //! @param B the ending line
    TransformationMatrix3D::TransformationMatrix3D( const coord::LineSegment3D &A, const coord::LineSegment3D &B) :
      m_TransformationMatrix3D( 4, 4, 0.0)
    {
      SetUnit();
      // get the angle and axis between the lines

      // get projected angle between lines A_A-A_B and B_A-B_B
      const double angle( linal::ProjAngle( A.GetStartPoint(), A.GetEndPoint(), B.GetStartPoint(), B.GetEndPoint()));

      // get rotation axis perpendicular to the plane containing two lines
      linal::Vector3D rotation_axis
      (
        linal::CrossProduct( A.GetDirection(), B.GetDirection())
      );

      // translate the midpoint on the line between the two points to the origin
      const linal::Vector3D mid_point_b( ( B.GetStartPoint() + B.GetEndPoint()) / 2.0);
      ( *this)( -mid_point_b);

      // handle collinear case
      // Choose a random axis that is orthogonal to them
      while( rotation_axis.SquareNorm() == 0.0)
      {
        rotation_axis.SetRandomTranslation( 1.0);
        rotation_axis = linal::CrossProduct( rotation_axis, A.GetDirection());
      }

      RotationMatrix3D rotation( rotation_axis, angle);
      ( *this)( rotation);

      // translate the midpoint back to the midpoint of a
      const linal::Vector3D mid_point_a( ( A.GetStartPoint() + A.GetEndPoint()) / 2.0);
      ( *this)( mid_point_a);
    }

    //! copy constructor
    TransformationMatrix3D *TransformationMatrix3D::Clone() const
    {
      return new TransformationMatrix3D( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &TransformationMatrix3D::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! access and modify the matrix
    void TransformationMatrix3D::SetMatrix( const linal::Matrix< double> &MATRIX)
    {
      BCL_Assert( MATRIX.IsSquare() && MATRIX.GetNumberRows() == size_t( 4), "Argument has to be 4x4 matrix");
      m_TransformationMatrix3D = MATRIX;
    }

    //! read only access
    const linal::Matrix< double> &TransformationMatrix3D::GetMatrix() const
    {
      return m_TransformationMatrix3D;
    }

  ///////////////
  // operators //
  ///////////////

    //! operator = TransformationsMatrix3D
    TransformationMatrix3D &TransformationMatrix3D::operator =( const TransformationMatrix3D &MATRIX)
    {
      m_TransformationMatrix3D = MATRIX.m_TransformationMatrix3D;
      return *this;
    }

    //! operator = TransformationsMatrix3D
    TransformationMatrix3D &TransformationMatrix3D::operator =( TransformationMatrix3D && MATRIX)
    {
      m_TransformationMatrix3D = std::move( MATRIX.m_TransformationMatrix3D);
      return *this;
    }

    //! operator *= TransformationMatrix3D
    TransformationMatrix3D &TransformationMatrix3D::operator *=( const TransformationMatrix3D &MATRIX)
    {
      m_TransformationMatrix3D = m_TransformationMatrix3D * MATRIX.m_TransformationMatrix3D;
      return *this;
    }

    //! operator += TRANSFORMATIONMATRIX3D
    TransformationMatrix3D &TransformationMatrix3D::operator +=( const TransformationMatrix3D &MATRIX)
    {
      m_TransformationMatrix3D += MATRIX.m_TransformationMatrix3D;
      return *this;
    }

    //! operator /= SCALAR
    TransformationMatrix3D &TransformationMatrix3D::operator /=( const double &SCALAR)
    {
      m_TransformationMatrix3D /= SCALAR;
      return *this;
    }

    //! operator *= SCALAR
    TransformationMatrix3D &TransformationMatrix3D::operator *=( const double &SCALAR)
    {
      m_TransformationMatrix3D *= SCALAR;
      return *this;
    }

    //! apply translation
    TransformationMatrix3D &TransformationMatrix3D::operator()( const TransformationMatrix3D &MATRIX)
    {
      return operator *=( TransformationMatrix3D( MATRIX));
    }

    //!  apply translation
    TransformationMatrix3D &TransformationMatrix3D::operator()( const double &X, const double &Y, const double &Z)
    {
      return operator *=( TransformationMatrix3D( X, Y, Z));
    }

    //!  apply translation
    TransformationMatrix3D &TransformationMatrix3D::operator()( const linal::Vector3D &TRANSLATION)
    {
      return operator *=( TransformationMatrix3D( TRANSLATION));
    }

    //! apply rotation
    TransformationMatrix3D &TransformationMatrix3D::operator()( const RotationMatrix3D &ROTATION)
    {
      return operator *=( TransformationMatrix3D( ROTATION));
    }

    //!  apply coord::Axis and RotationAngle
    TransformationMatrix3D &TransformationMatrix3D::operator()( const coord::Axis &AXIS, const double &ANGLE)
    {
      return operator *=( TransformationMatrix3D( AXIS, ANGLE));
    }

    //!  apply coord::Axis and three RotationAngles
    TransformationMatrix3D &TransformationMatrix3D::operator()( const coord::Axis AXES[ 3], const double ANGLE[ 3])
    {
      return operator *=( TransformationMatrix3D( AXES, ANGLE));
    }

  ////////////////
  // operations //
  ////////////////

    //! Set Unit Matrix
    TransformationMatrix3D &TransformationMatrix3D::SetUnit()
    {
      m_TransformationMatrix3D.SetZero();
      SetUnitFromZero();
      return *this;
    }

    //! Inverse
    TransformationMatrix3D &TransformationMatrix3D::Invert()
    {
      linal::MatrixInversionGaussJordan< double> inverter( m_TransformationMatrix3D);
      BCL_Assert( inverter.IsDefined(), "Could not compute the inverse of " + util::Format()( m_TransformationMatrix3D));
      m_TransformationMatrix3D = inverter.ComputeInverse();
      return *this;
    }

    //! return rotationmatrix
    RotationMatrix3D TransformationMatrix3D::GetRotation() const
    {
      linal::Matrix3x3< double> rotation;
      rotation.ReplaceRow( 0, linal::VectorConstReference< double>( 3, m_TransformationMatrix3D[ 0]));
      rotation.ReplaceRow( 1, linal::VectorConstReference< double>( 3, m_TransformationMatrix3D[ 1]));
      rotation.ReplaceRow( 2, linal::VectorConstReference< double>( 3, m_TransformationMatrix3D[ 2]));
      return RotationMatrix3D( rotation);
    }

    //! return Translation
    linal::Vector3D TransformationMatrix3D::GetTranslation() const
    {
      return linal::Vector3D( m_TransformationMatrix3D[ 3]);
    }

    linal::Vector3D TransformationMatrix3D::GetAxis( const coord::Axis &AXIS) const
    {
      return linal::Vector3D( m_TransformationMatrix3D[ AXIS]);
    }

    linal::Vector3D TransformationMatrix3D::GetOrigin() const
    {
      return linal::Vector3D( m_TransformationMatrix3D[ 3]);
    }

    //! @brief set the translation for this transformation matrix
    //! @param TRANSLATION translation vector
    void TransformationMatrix3D::SetTranslation( const linal::Vector3D &TRANSLATION)
    {
      m_TransformationMatrix3D.ReplaceRow( 3, TRANSLATION);
    }

    //! @brief set the rotation matrix for this class
    //! @param ROTATION the rotation matrix to use
    void TransformationMatrix3D::SetRotation( const RotationMatrix3D &ROTATION)
    {
      m_TransformationMatrix3D.ReplaceRow( 0, ROTATION.GetMatrix().GetRow( 0));
      m_TransformationMatrix3D.ReplaceRow( 1, ROTATION.GetMatrix().GetRow( 1));
      m_TransformationMatrix3D.ReplaceRow( 2, ROTATION.GetMatrix().GetRow( 2));
    }

    //! modification of the axis
    void TransformationMatrix3D::SetAxis( const coord::Axis &AXIS, const linal::Vector3D &NEW_AXIS)
    {
      for( size_t i = 0; i < 3; ++i)
      {
        m_TransformationMatrix3D( AXIS, i) = NEW_AXIS( i);
      }
    }

    void TransformationMatrix3D::SetOrigin( const linal::Vector3D &AXIS)
    {
      for( size_t i = 0; i < 3; ++i)
      {
        m_TransformationMatrix3D( 3, i) = AXIS( i);
      }
    }

    //! @brief returns whether this transformation is defined or not by looking at the scaling factor
    //! @return whether this transformation is defined or not
    bool TransformationMatrix3D::IsDefined() const
    {
      // return whether scaling factor is defined
      return util::IsDefined( m_TransformationMatrix3D( 3, 3));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write TransformationMatrix3D to std::ostream using the given util::Format
    std::ostream &TransformationMatrix3D::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_TransformationMatrix3D, OSTREAM, INDENT) << '\n';
      return OSTREAM;
    }

    //! write TransformationMatrix3D to std::ostream using the given util::Format
    std::ostream &TransformationMatrix3D::Write( std::ostream &OSTREAM, const size_t INDENT, const util::Format &FORMAT) const
    {
      io::Serialize::Write( m_TransformationMatrix3D, OSTREAM, INDENT, FORMAT) << '\n';
      return OSTREAM;
    }

    //! read TransformationMatrix3D from io::IFStream
    std::istream &TransformationMatrix3D::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_TransformationMatrix3D, ISTREAM);
      return ISTREAM;
    }

    //! Set Unitary matrix, assuming that the matrix is currently composed of zeros
    void TransformationMatrix3D::SetUnitFromZero()
    {
      m_TransformationMatrix3D( 0, 0) = 1.0;
      m_TransformationMatrix3D( 1, 1) = 1.0;
      m_TransformationMatrix3D( 2, 2) = 1.0;
      m_TransformationMatrix3D( 3, 3) = 1.0;
    }

    //! boolean operator TRANSFORMATIONMATRIX3D_A == TRANSFORMATIONMATRIX3D_B
    bool operator ==
    (
      const TransformationMatrix3D &TRANSFORMATIONMATRIX3D_A,
      const TransformationMatrix3D &TRANSFORMATIONMATRIX3D_B
    )
    {
      return
        std::equal
        (
          TRANSFORMATIONMATRIX3D_A.GetMatrix().Begin(),
          TRANSFORMATIONMATRIX3D_A.GetMatrix().End(),
          TRANSFORMATIONMATRIX3D_B.GetMatrix().Begin()
        );
    }

    //! boolean operator TRANSFORMATIONMATRIX3D_A != TRANSFORMATIONMATRIX3D_B
    bool operator !=
    (
      const TransformationMatrix3D &TRANSFORMATIONMATRIX3D_A,
      const TransformationMatrix3D &TRANSFORMATIONMATRIX3D_B
    )
    {
      return !( TRANSFORMATIONMATRIX3D_A == TRANSFORMATIONMATRIX3D_B);
    }

    //! return Inverse of TransformationMatrix3D
    TransformationMatrix3D Inverse( const TransformationMatrix3D &MATRIX)
    {
      return TransformationMatrix3D( MATRIX).Invert();
    }

    //! check if two transformation matrices are similar within a given rotational and translational tolerance
    bool SimilarWithinTolerance
    (
      const TransformationMatrix3D &TRANSFORMATIONMATRIX3D_A,
      const TransformationMatrix3D &TRANSFORMATIONMATRIX3D_B,
      const double TRANSLATION_TOLERANCE,
      const double ROTATION_TOLERANCE_RAD
    )
    {
      TransformationMatrix3D diff( TRANSFORMATIONMATRIX3D_A);
      diff( Inverse( TRANSFORMATIONMATRIX3D_B));

      // check for similar Translation and similar rotation
      return (
                 ( diff.GetTranslation().Norm() < TRANSLATION_TOLERANCE)
              && ( diff.GetRotation().EffectiveRotationAngle() < ROTATION_TOLERANCE_RAD)
             );
    }

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_tricubic_spline.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"
#include "math/bcl_math_bicubic_spline.h"
#include "math/bcl_math_cubic_spline.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> TricubicSpline::s_Instance
    (
      GetObjectInstances().AddInstance( new TricubicSpline())
    );

  ////////////////
  // operations //
  ////////////////

    //! train TricubicSpline
    void TricubicSpline::Train
    (
      const SplineBorderType BORDER[ 3], const double START[ 3], const double DELTA[ 3], const Tensor< double> &RESULTS,
      const bool LINCONT[ 3], const storage::Pair< double, double> FIRSTBE[ 3]
    )
    {
      //check, if the points are given in positive direction
      BCL_Assert( DELTA[ coord::GetAxes().e_X]>0 , "deltax <= 0 not supported");
      BCL_Assert( DELTA[ coord::GetAxes().e_Y]>0 , "deltay <= 0 not supported");
      BCL_Assert( DELTA[ coord::GetAxes().e_Z]>0 , "deltaz <= 0 not supported");

      //determine values for all three dimensions
      const int dimz( RESULTS.GetNumberCols());
      const int dimy( RESULTS.GetNumberRows());
      const int dimx( RESULTS.NumberLayers());

      //assigning values
      m_Border[ coord::GetAxes().e_X]  = BORDER[ coord::GetAxes().e_X];
      m_Border[ coord::GetAxes().e_Y]  = BORDER[ coord::GetAxes().e_Y];
      m_Border[ coord::GetAxes().e_Z]  = BORDER[ coord::GetAxes().e_Z];
      m_Start[ coord::GetAxes().e_X]   = START[ coord::GetAxes().e_X];
      m_Start[ coord::GetAxes().e_Y]   = START[ coord::GetAxes().e_Y];
      m_Start[ coord::GetAxes().e_Z]   = START[ coord::GetAxes().e_Z];
      m_Delta[ coord::GetAxes().e_X]   = DELTA[ coord::GetAxes().e_X];
      m_Delta[ coord::GetAxes().e_Y]   = DELTA[ coord::GetAxes().e_Y];
      m_Delta[ coord::GetAxes().e_Z]   = DELTA[ coord::GetAxes().e_Z];

      BCL_Assert
      (
        m_Border[0] != e_NotAKnot &&
        m_Border[1] != e_NotAKnot &&
        m_Border[2] != e_NotAKnot ,
        "The specified boundary condition is not yet supported"
      );

      m_Values   = RESULTS;
      m_Dsecox   = RESULTS;
      m_Dsecoy   = RESULTS;
      m_Dsecoz   = RESULTS;
      m_Dsecoxy  = RESULTS;
      m_Dsecoxz  = RESULTS;
      m_Dsecoyz  = RESULTS;
      m_Dsecoxyz = RESULTS;

      m_LinCont[ coord::GetAxes().e_X] = LINCONT[ coord::GetAxes().e_X];
      m_LinCont[ coord::GetAxes().e_Y] = LINCONT[ coord::GetAxes().e_Y];
      m_LinCont[ coord::GetAxes().e_Z] = LINCONT[ coord::GetAxes().e_Z];
      m_FirstBe[ coord::GetAxes().e_X] = FIRSTBE[ coord::GetAxes().e_X];
      m_FirstBe[ coord::GetAxes().e_Y] = FIRSTBE[ coord::GetAxes().e_Y];
      m_FirstBe[ coord::GetAxes().e_Z] = FIRSTBE[ coord::GetAxes().e_Z];

      //train seven times for fxx, fyy, fzz, fxxyy, fxxzz, fyyzz, fxxyyzz
      //reduction to Spline2D by training only 2D-layers at the same time
      for( int layer( 0); layer < dimx; ++layer)
      {
        linal::Matrix< double> values( dimy, dimz);
        for( int row( 0); row < dimy; ++row)
        {
          for( int col( 0); col < dimz; ++col)
          {
            values( row, col) = m_Values( layer, row, col);
          }
        }
        BicubicSpline bs;
        SplineBorderType border[2]            = { BORDER[ coord::GetAxes().e_Y],   BORDER[ coord::GetAxes().e_Z]};
        double start[2]                       = { START[ coord::GetAxes().e_Y],    START[ coord::GetAxes().e_Z]};
        double delta[2]                       = { DELTA[ coord::GetAxes().e_Y],    DELTA[ coord::GetAxes().e_Z]};
        bool lin_cont[2]                      = { LINCONT[ coord::GetAxes().e_Y],  LINCONT[ coord::GetAxes().e_Z]};
        storage::Pair< double, double> firstbe[ 2] = { FIRSTBE[ coord::GetAxes().e_Y], FIRSTBE[ coord::GetAxes().e_Z]};

        bs.Train( border, start, delta, values, lin_cont, firstbe);
        m_Dsecoy.ReplaceLayer(  layer, bs.GetDsecox());
        m_Dsecoz.ReplaceLayer(  layer, bs.GetDsecoy());
        m_Dsecoyz.ReplaceLayer( layer, bs.GetDsecoxy());
      }

      for( int row( 0); row < dimy; ++row)
      {
        for( int col( 0); col < dimz; ++col)
        {
          linal::Vector< double> values( dimx), dsecoz( dimx), dsecoy( dimx), dsecoyz( dimx);
          for( int layer( 0); layer < dimx; ++layer)
          {
            values(   layer) = m_Values(  layer, row, col);
            dsecoy(   layer) = m_Dsecoy(  layer, row, col);
            dsecoz(   layer) = m_Dsecoz(  layer, row, col);
            dsecoyz(  layer) = m_Dsecoyz( layer, row, col);
          }
          CubicSpline cs, csz, csy, csyz;
          cs.Train(   BORDER[ coord::GetAxes().e_X], START[ coord::GetAxes().e_X], DELTA[ coord::GetAxes().e_X], values , FIRSTBE[ coord::GetAxes().e_X]);
          csy.Train(  BORDER[ coord::GetAxes().e_X], START[ coord::GetAxes().e_X], DELTA[ coord::GetAxes().e_X], dsecoy , FIRSTBE[ coord::GetAxes().e_X]);
          csz.Train(  BORDER[ coord::GetAxes().e_X], START[ coord::GetAxes().e_X], DELTA[ coord::GetAxes().e_X], dsecoz , FIRSTBE[ coord::GetAxes().e_X]);
          csyz.Train( BORDER[ coord::GetAxes().e_X], START[ coord::GetAxes().e_X], DELTA[ coord::GetAxes().e_X], dsecoyz, FIRSTBE[ coord::GetAxes().e_X]);
          for( int layer( 0); layer < dimx; ++layer)
          {
            m_Dsecox(   layer, row, col) = cs.GetDsecox()(    layer);
            m_Dsecoxy(  layer, row, col) = csy.GetDsecox()(   layer);
            m_Dsecoxz(  layer, row, col) = csz.GetDsecox()(   layer);
            m_Dsecoxyz( layer, row, col) = csyz.GetDsecox()(  layer);
          }
        }
      }
      return;
    }

    //! return value at certain (x, y, z)
    double TricubicSpline::F( const double x, const double y, const double z) const
    {
      const int dimx( m_Values.NumberLayers());
      const int dimy( m_Values.GetNumberRows());
      const int dimz( m_Values.GetNumberCols());

      //check if argument is in range for non-periodic splines
      if
      (
            ( m_Border[ coord::GetAxes().e_X] != e_Periodic)
        && ( ( x < m_Start[ coord::GetAxes().e_X]) || ( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X] < x))
      )
      {
        BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");
        BCL_MessageDbg( "argument out of range, using linear continuation");
        if( x < m_Start[ coord::GetAxes().e_X])
        {
          return F( m_Start[ coord::GetAxes().e_X], y, z) + ( x - m_Start[ coord::GetAxes().e_X]) * dFdx( m_Start[ coord::GetAxes().e_X], y, z);
        }
        if( x > m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X])
        {
          return F( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X], y, z)
                + ( x - m_Start[ coord::GetAxes().e_X] - ( dimx - 1) * m_Delta[ coord::GetAxes().e_X])
                * dFdx( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X], y, z);
        }
      }

      //check if argument y is in range for non-periodic splines
      if
      (
            ( m_Border[ coord::GetAxes().e_Y] != e_Periodic)
        && ( ( y < m_Start[ coord::GetAxes().e_Y]) || ( m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y] < y))
      )
      {
        BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");
        BCL_MessageDbg( "argument out of range, using linear continuation");
        if( y < m_Start[ coord::GetAxes().e_Y])
        {
          return F( x, m_Start[ coord::GetAxes().e_Y], z) + ( y - m_Start[ coord::GetAxes().e_Y]) * dFdy( x, m_Start[ coord::GetAxes().e_Y], z);
        }
        if( y > m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y])
        {
          return F( x, m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y], z)
                + ( y - m_Start[ coord::GetAxes().e_Y] - ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y])
                * dFdy( x, m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y], z);
        }
      }

      //check if argument z is in range for non-periodic splines
      if
      (
            ( m_Border[ coord::GetAxes().e_Z] != e_Periodic)
        && ( ( z < m_Start[ coord::GetAxes().e_Z]) || ( m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z] < z))
      )
      {
        BCL_Assert( m_LinCont[ coord::GetAxes().e_Z], "argument out of range for non-periodic spline!");
        BCL_MessageDbg( "argument out of range, using linear continuation");
        if( z < m_Start[ coord::GetAxes().e_Z])
        {
          return F( x, y, m_Start[ coord::GetAxes().e_Z]) + ( z - m_Start[ coord::GetAxes().e_Z]) * dFdz( x, y, m_Start[ coord::GetAxes().e_Z]);
        }
        if( z > m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z])
        {
          return F( x, y, m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z])
                + ( z - m_Start[ coord::GetAxes().e_Z] - ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z])
                * dFdz( x, y, m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z]);
        }
      }

      //determine i with m_Start[ coord::GetAxes().e_X]+(i-1)*m_Delta[ coord::GetAxes().e_X] < x < m_Start[ coord::GetAxes().e_X]+i*m_Delta[ coord::GetAxes().e_X] for the correct supporting points
      int i( int( floor( ( x - m_Start[ coord::GetAxes().e_X]) / m_Delta[ coord::GetAxes().e_X])));
      while( m_Start[ coord::GetAxes().e_X] + i * m_Delta[ coord::GetAxes().e_X] < x) i++;

      //determine j with m_Start[ coord::GetAxes().e_Y]+(i-1)*m_Delta[ coord::GetAxes().e_Y] < y < m_Start[ coord::GetAxes().e_Y]+i*m_Delta[ coord::GetAxes().e_Y] for the correct supporting points
      int j( int( floor( ( y - m_Start[ coord::GetAxes().e_Y]) / m_Delta[ coord::GetAxes().e_Y])));
      while( m_Start[ coord::GetAxes().e_Y] + j * m_Delta[ coord::GetAxes().e_Y] < y) j++;

      // determine k with m_Start[ coord::GetAxes().e_Z]+(k-1)*m_Deltaz < z < m_Start[ coord::GetAxes().e_Z]+k*m_Deltaz for the correct supporting points
      int k( int( floor( ( z - m_Start[ coord::GetAxes().e_Z]) / m_Delta[ coord::GetAxes().e_Z])));
      while( m_Start[ coord::GetAxes().e_Z] + k * m_Delta[ coord::GetAxes().e_Z] < z) k++;

      // this formula was derived from the Spline2D formula
      // the idea is to 'combine every part of the 2D formula
      // with dzm, dzp, dz3m, dz3p'

      const double delta_aktx( x - m_Start[ coord::GetAxes().e_X] - (i-1) * m_Delta[ coord::GetAxes().e_X]);
      const double delta_akty( y - m_Start[ coord::GetAxes().e_Y] - (j-1) * m_Delta[ coord::GetAxes().e_Y]);
      const double delta_aktz( z - m_Start[ coord::GetAxes().e_Z] - (k-1) * m_Delta[ coord::GetAxes().e_Z]);

      const double dxp( delta_aktx / m_Delta[ coord::GetAxes().e_X]);
      const double dxm( 1 - dxp);
      const double dx3p( ( dxp*dxp*dxp - dxp) * Sqr( m_Delta[ coord::GetAxes().e_X]) / 6);
      const double dx3m( ( dxm*dxm*dxm - dxm) * Sqr( m_Delta[ coord::GetAxes().e_X]) / 6);

      const double dyp( delta_akty / m_Delta[ coord::GetAxes().e_Y]);
      const double dym( 1 - dyp);
      const double dy3p( ( dyp*dyp*dyp - dyp) * Sqr( m_Delta[ coord::GetAxes().e_Y]) / 6);
      const double dy3m( ( dym*dym*dym - dym) * Sqr( m_Delta[ coord::GetAxes().e_Y]) / 6);

      const double dzp( delta_aktz / m_Delta[ coord::GetAxes().e_Z]);
      const double dzm( 1 - dzp);
      const double dz3p( ( dzp*dzp*dzp - dzp) * Sqr( m_Delta[ coord::GetAxes().e_Z]) / 6);
      const double dz3m( ( dzm*dzm*dzm - dzm) * Sqr( m_Delta[ coord::GetAxes().e_Z]) / 6);

      //generate positive values to prevent some problems with the indices
      while( i<1) i += dimx;
      while( j<1) j += dimy;
      while( k<1) k += dimz;

      return
        dzm
        *(dxm*(dym*m_Values( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Values( (i-1)%dimx  , j%dimy, (k-1)%dimz))
        + dxp*(dym*m_Values(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Values(i%dimx      , j%dimy, (k-1)%dimz))
        +dx3m*(dym*m_Dsecox( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecox( (i-1)%dimx  , j%dimy, (k-1)%dimz))
        +dx3p*(dym*m_Dsecox(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecox(i%dimx      , j%dimy, (k-1)%dimz))
        + dxm*(dy3m*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoy( (i-1)%dimx , j%dimy, (k-1)%dimz))
        + dxp*(dy3m*m_Dsecoy(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoy(i%dimx     , j%dimy, (k-1)%dimz))
        +dx3m*(dy3m*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxy( (i-1)%dimx, j%dimy, (k-1)%dimz))
        +dx3p*(dy3m*m_Dsecoxy(i%dimx    , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxy(i%dimx    , j%dimy, (k-1)%dimz)))

        +dzp
        *(dxm*(dym*m_Values( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Values( (i-1)%dimx  , j%dimy, k%dimz))
        + dxp*(dym*m_Values(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Values(i%dimx      , j%dimy, k%dimz))
        +dx3m*(dym*m_Dsecox( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecox( (i-1)%dimx  , j%dimy, k%dimz))
        +dx3p*(dym*m_Dsecox(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecox(i%dimx      , j%dimy, k%dimz))
        + dxm*(dy3m*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoy( (i-1)%dimx , j%dimy, k%dimz))
        + dxp*(dy3m*m_Dsecoy(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoy(i%dimx     , j%dimy, k%dimz))
        +dx3m*(dy3m*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxy( (i-1)%dimx, j%dimy, k%dimz))
        +dx3p*(dy3m*m_Dsecoxy(i%dimx    , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxy(i%dimx    , j%dimy, k%dimz)))

        +dz3m
        *(dxm*(dym*m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoz( (i-1)%dimx  , j%dimy, (k-1)%dimz))
        + dxp*(dym*m_Dsecoz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoz(i%dimx      , j%dimy, (k-1)%dimz))
        +dx3m*(dym*m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoxz( (i-1)%dimx  , j%dimy, (k-1)%dimz))
        +dx3p*(dym*m_Dsecoxz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoxz(i%dimx      , j%dimy, (k-1)%dimz))
        + dxm*(dy3m*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoyz( (i-1)%dimx , j%dimy, (k-1)%dimz))
        + dxp*(dy3m*m_Dsecoyz(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoyz(i%dimx     , j%dimy, (k-1)%dimz))
        +dx3m*(dy3m*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxyz( (i-1)%dimx, j%dimy, (k-1)%dimz))
        +dx3p*(dy3m*m_Dsecoxyz(i%dimx    , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxyz(i%dimx    , j%dimy, (k-1)%dimz)))

        +dz3p
        *(dxm*(dym*m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecoz( (i-1)%dimx  , j%dimy, k%dimz))
        + dxp*(dym*m_Dsecoz(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecoz(i%dimx      , j%dimy, k%dimz))
        +dx3m*(dym*m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecoxz( (i-1)%dimx  , j%dimy, k%dimz))
        +dx3p*(dym*m_Dsecoxz(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecoxz(i%dimx      , j%dimy, k%dimz))
        + dxm*(dy3m*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoyz( (i-1)%dimx , j%dimy, k%dimz))
        + dxp*(dy3m*m_Dsecoyz(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoyz(i%dimx     , j%dimy, k%dimz))
        +dx3m*(dy3m*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxyz( (i-1)%dimx, j%dimy, k%dimz))
        +dx3p*(dy3m*m_Dsecoxyz(i%dimx    , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxyz(i%dimx    , j%dimy, k%dimz)))
      ;
    }

    //! return partial derivative at certain (x, y, z) for x
    double TricubicSpline::dFdx( const double x, const double y, const double z) const
    {
      const int dimx( m_Values.NumberLayers());
      const int dimy( m_Values.GetNumberRows());
      const int dimz( m_Values.GetNumberCols());

      //check if argument x is in range for non-periodic splines
      if
      (
            ( m_Border[ coord::GetAxes().e_X] != e_Periodic)
        && ( ( x < m_Start[ coord::GetAxes().e_X]) || ( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X] < x))
      )
      {
        BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");
        BCL_MessageDbg( "argument out of range, using linear continuation");
        if( x < m_Start[ coord::GetAxes().e_X])
        {
          return dFdx( m_Start[ coord::GetAxes().e_X], y, z);
        }
        if( x > m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X])
        {
          return dFdx( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X], y, z);
        }
      }

      //check if argument y is in range for non-periodic splines
      if
      (
            ( m_Border[ coord::GetAxes().e_Y] != e_Periodic)
        && ( ( y < m_Start[ coord::GetAxes().e_Y] || m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y] < y))
      )
      {
        BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");
        BCL_MessageDbg( "argument out of range, using linear continuation");
        if( y < m_Start[ coord::GetAxes().e_Y])
        {
          return dFdx( x, m_Start[ coord::GetAxes().e_Y], z);
        }
        if( y > m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y])
        {
          return dFdx( x, m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y], z);
        }
      }

      //check if argument z is in range for non-periodic splines
      if
      (
            ( m_Border[ coord::GetAxes().e_Z] != e_Periodic)
        && ( ( z < m_Start[ coord::GetAxes().e_Z] || m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z] < z)))
      {
        BCL_Assert( m_LinCont[ coord::GetAxes().e_Z], "argument out of range for non-periodic spline!");
        BCL_MessageDbg( "argument out of range, using linear continuation");
        if( z < m_Start[ coord::GetAxes().e_Z])
        {
          return dFdx( x, y, m_Start[ coord::GetAxes().e_Z]);
        }
        if( z > m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z])
        {
          return dFdx( x, y, m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z]);
        }
      }

      //determine i with m_Start+(i-1)*m_Delta < x < m_Start+i*m_Delta for the correct supporting points
      int i( int( floor( ( x - m_Start[ coord::GetAxes().e_X]) / m_Delta[ coord::GetAxes().e_X])));
      while( m_Start[ coord::GetAxes().e_X] + i * m_Delta[ coord::GetAxes().e_X] < x) i++;

      //the same for j
      int j( int( floor( ( y - m_Start[ coord::GetAxes().e_Y]) / m_Delta[ coord::GetAxes().e_Y])));
      while( m_Start[ coord::GetAxes().e_Y] + j * m_Delta[ coord::GetAxes().e_Y] < y) j++;

      //the same for k
      int k( int( floor( ( z - m_Start[ coord::GetAxes().e_Z]) / m_Delta[ coord::GetAxes().e_Z])));
      while( m_Start[ coord::GetAxes().e_Z] + k * m_Delta[ coord::GetAxes().e_Z] < z) k++;

      const double delta_aktx( x - m_Start[ coord::GetAxes().e_X] - ( i - 1) * m_Delta[ coord::GetAxes().e_X]);
      const double delta_akty( y - m_Start[ coord::GetAxes().e_Y] - ( j - 1) * m_Delta[ coord::GetAxes().e_Y]);
      double delta_aktz( z - m_Start[ coord::GetAxes().e_Z] - ( k - 1) * m_Delta[ coord::GetAxes().e_Z]);

      const double dxp( delta_aktx / m_Delta[ coord::GetAxes().e_X]);
      const double dxm( 1 - dxp);

      const double dyp( delta_akty / m_Delta[ coord::GetAxes().e_Y]);
      const double dym( 1 - dyp);
      const double dy3p( ( dyp * dyp * dyp - dyp) * Sqr( m_Delta[ coord::GetAxes().e_Y]) / 6);
      const double dy3m( ( dym * dym * dym - dym) * Sqr( m_Delta[ coord::GetAxes().e_Y]) / 6);

      const double dzp( delta_aktz / m_Delta[ coord::GetAxes().e_Z]);
      const double dzm( 1 - dzp);
      const double dz3p( ( dzp * dzp * dzp - dzp) * Sqr( m_Delta[ coord::GetAxes().e_Z]) / 6);
      const double dz3m( ( dzm * dzm * dzm - dzm) * Sqr( m_Delta[ coord::GetAxes().e_Z]) / 6);

      //generate positive values to prevent some problems with the indizes
      while( i<1) i += dimx;
      while( j<1) j += dimy;
      while( k<1) k += dimz;

      return
        dzm
        * (
          -(dym*m_Values( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Values( (i-1)%dimx  , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dym*m_Values(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Values(i%dimx      , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecox( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecox( (i-1)%dimx  , j%dimy, (k-1)%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecox(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecox(i%dimx      , j%dimy, (k-1)%dimz))
          -(dy3m*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoy( (i-1)%dimx , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dy3m*m_Dsecoy(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoy(i%dimx     , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxy( (i-1)%dimx, j%dimy, (k-1)%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxy(i%dimx    , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxy(i%dimx    , j%dimy, (k-1)%dimz))
        )
        +dzp
        *(
          -(dym*m_Values( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Values( (i-1)%dimx  , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dym*m_Values(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Values(i%dimx      , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecox( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecox( (i-1)%dimx  , j%dimy, k%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecox(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecox(i%dimx      , j%dimy, k%dimz))
          -(dy3m*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoy( (i-1)%dimx , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dy3m*m_Dsecoy(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoy(i%dimx     , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxy( (i-1)%dimx, j%dimy, k%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxy(i%dimx    , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxy(i%dimx    , j%dimy, k%dimz))
        )
        +dz3m
        *(
          -(dym*m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoz( (i-1)%dimx  , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dym*m_Dsecoz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoz(i%dimx      , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoxz( (i-1)%dimx  , j%dimy, (k-1)%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecoxz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoxz(i%dimx      , j%dimy, (k-1)%dimz))
          -(dy3m*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoyz( (i-1)%dimx , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dy3m*m_Dsecoyz(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoyz(i%dimx     , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxyz( (i-1)%dimx, j%dimy, (k-1)%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxyz(i%dimx    , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxyz(i%dimx    , j%dimy, (k-1)%dimz))
        )
        +dz3p
        *(
          -(dym*m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecoz( (i-1)%dimx  , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dym*m_Dsecoz(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecoz(i%dimx      , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecoxz( (i-1)%dimx  , j%dimy, k%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecoxz(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecoxz(i%dimx      , j%dimy, k%dimz))
          -(dy3m*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoyz( (i-1)%dimx , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dy3m*m_Dsecoyz(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoyz(i%dimx     , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxyz( (i-1)%dimx, j%dimy, k%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxyz(i%dimx    , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxyz(i%dimx    , j%dimy, k%dimz))
        )
      ;
    }

    //! return partial derivative at certain (x, y, z) for y
    double TricubicSpline::dFdy( const double x, const double y, const double z) const
    {
      const int dimx( m_Values.NumberLayers());
      const int dimy( m_Values.GetNumberRows());
      const int dimz( m_Values.GetNumberCols());

      //check if argument x is in range for non-periodic splines
      if
      (
            ( m_Border[ coord::GetAxes().e_X] != e_Periodic)
        && ( ( x < m_Start[ coord::GetAxes().e_X]) || ( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X] < x)))
      {
        BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");
        BCL_MessageDbg( "argument out of range, using linear continuation");
        if( x < m_Start[ coord::GetAxes().e_X])
        {
          return dFdy( m_Start[ coord::GetAxes().e_X], y, z);
        }
        if( x > m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X])
        {
          return dFdy( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X], y, z);
        }
      }

      //check if argument y is in range for non-periodic splines
      if
      (
            ( m_Border[ coord::GetAxes().e_Y] != e_Periodic)
        && ( ( y < m_Start[ coord::GetAxes().e_Y] || m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y] < y))
      )
      {
        BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");
        BCL_MessageDbg( "argument out of range, using linear continuation");
        if( y < m_Start[ coord::GetAxes().e_Y])
        {
          return dFdy( x, m_Start[ coord::GetAxes().e_Y], z);
        }
        if( y > m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y])
        {
          return dFdy( x, m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y], z);
        }
      }

      //check if argument z is in range for non-periodic splines
      if
      (
            ( m_Border[ coord::GetAxes().e_Z] != e_Periodic)
        && ( ( z < m_Start[ coord::GetAxes().e_Z] || m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z] < z))
      )
      {
        BCL_Assert( m_LinCont[ coord::GetAxes().e_Z], "argument out of range for non-periodic spline!");
        BCL_MessageDbg( "argument out of range, using linear continuation");
        if( z < m_Start[ coord::GetAxes().e_Z])
        {
          return dFdy( x, y, m_Start[ coord::GetAxes().e_Z]);
        }
        if( z > m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z])
        {
          return dFdy( x, y, m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z]);
        }
      }

      //determine i with m_Start[ coord::GetAxes().e_X]+(i-1)*m_Delta[ coord::GetAxes().e_X] < x < m_Start[ coord::GetAxes().e_X]+i*m_Delta[ coord::GetAxes().e_X]
      //for the correct supporting points
      int i( int( floor( ( x - m_Start[ coord::GetAxes().e_X]) / m_Delta[ coord::GetAxes().e_X])));
      while( m_Start[ coord::GetAxes().e_X] + i * m_Delta[ coord::GetAxes().e_X] < x) i++;

      //determine j with m_Start[ coord::GetAxes().e_Y]+(j-1)*m_Delta[ coord::GetAxes().e_Y] < y < m_Start[ coord::GetAxes().e_Y]+j*m_Delta[ coord::GetAxes().e_Y]
      //for the correct supporting points
      int j( int( floor( ( y - m_Start[ coord::GetAxes().e_Y]) / m_Delta[ coord::GetAxes().e_Y])));
      while( m_Start[ coord::GetAxes().e_Y] + j * m_Delta[ coord::GetAxes().e_Y] < y) j++;

      //determine k with m_Start[ coord::GetAxes().e_Z]+(k-1)*m_Delta[ coord::GetAxes().e_Z] < z < m_Start[ coord::GetAxes().e_Z]+k*m_Delta[ coord::GetAxes().e_Z]
      //for the correct supporting points
      int k( int( floor( ( z - m_Start[ coord::GetAxes().e_Z]) / m_Delta[ coord::GetAxes().e_Z])));
      while( m_Start[ coord::GetAxes().e_Z] + k * m_Delta[ coord::GetAxes().e_Z] < z) k++;

      //generate some auxiliary variables
      const double delta_aktx( x - m_Start[ coord::GetAxes().e_X] - ( i - 1) * m_Delta[ coord::GetAxes().e_X]);
      const double delta_akty( y - m_Start[ coord::GetAxes().e_Y] - ( j - 1) * m_Delta[ coord::GetAxes().e_Y]);
      const double delta_aktz( z - m_Start[ coord::GetAxes().e_Z] - ( k - 1) * m_Delta[ coord::GetAxes().e_Z]);

      //see F(x, y, z) for a short explanation of the values
      const double dxp( delta_aktx / m_Delta[ coord::GetAxes().e_X]);
      const double dxm( 1 - dxp);
      const double dx3p( ( dxp * dxp * dxp - dxp) * Sqr( m_Delta[ coord::GetAxes().e_X]) / 6);
      const double dx3m( ( dxm * dxm * dxm - dxm) * Sqr( m_Delta[ coord::GetAxes().e_X]) / 6);

      const double dyp( delta_akty / m_Delta[ coord::GetAxes().e_Y]);
      const double dym( 1 - dyp);

      const double dzp( delta_aktz / m_Delta[ coord::GetAxes().e_Z]);
      const double dzm( 1 - dzp);
      const double dz3p( ( dzp * dzp * dzp - dzp) * Sqr( m_Delta[ coord::GetAxes().e_Z]) / 6);
      const double dz3m( ( dzm * dzm * dzm - dzm) * Sqr( m_Delta[ coord::GetAxes().e_Z]) / 6);

      //generate positive values to prevent some problems with the indizes
      while( i <= 0) i += dimx;
      while( j <= 0) j += dimy;
      while( k <= 0) k += dimz;

      return
        dzm
        * (
            dxm*(-m_Values( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+m_Values( (i-1)%dimx  , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxp*(-m_Values(i%dimx      , (j-1)%dimy, (k-1)%dimz)+m_Values(i%dimx      , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3m*(-m_Dsecox( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+m_Dsecox( (i-1)%dimx  , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3p*(-m_Dsecox(i%dimx      , (j-1)%dimy, (k-1)%dimz)+m_Dsecox(i%dimx      , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxm*(-(3*dym*dym-1)*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoy( (i-1)%dimx , j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          + dxp*(-(3*dym*dym-1)*m_Dsecoy(i%dimx     , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoy(i%dimx     , j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3m*(-(3*dym*dym-1)*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoxy( (i-1)%dimx, j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3p*(-(3*dym*dym-1)*m_Dsecoxy(i%dimx    , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoxy(i%dimx    , j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
        )
        +dzp
        *(
            dxm*(-m_Values( (i-1)%dimx  , (j-1)%dimy, k%dimz)+m_Values( (i-1)%dimx  , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxp*(-m_Values(i%dimx      , (j-1)%dimy, k%dimz)+m_Values(i%dimx      , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3m*(-m_Dsecox( (i-1)%dimx  , (j-1)%dimy, k%dimz)+m_Dsecox( (i-1)%dimx  , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3p*(-m_Dsecox(i%dimx      , (j-1)%dimy, k%dimz)+m_Dsecox(i%dimx      , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxm*(-(3*dym*dym-1)*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoy( (i-1)%dimx , j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          + dxp*(-(3*dym*dym-1)*m_Dsecoy(i%dimx     , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoy(i%dimx     , j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3m*(-(3*dym*dym-1)*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoxy( (i-1)%dimx, j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3p*(-(3*dym*dym-1)*m_Dsecoxy(i%dimx    , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoxy(i%dimx    , j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
        )
        +dz3m
        *(
            dxm*(-m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+m_Dsecoz( (i-1)%dimx  , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxp*(-m_Dsecoz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+m_Dsecoz(i%dimx      , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3m*(-m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+m_Dsecoxz( (i-1)%dimx  , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3p*(-m_Dsecoxz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+m_Dsecoxz(i%dimx      , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxm*(-(3*dym*dym-1)*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoyz( (i-1)%dimx , j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          + dxp*(-(3*dym*dym-1)*m_Dsecoyz(i%dimx     , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoyz(i%dimx     , j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3m*(-(3*dym*dym-1)*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoxyz( (i-1)%dimx, j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3p*(-(3*dym*dym-1)*m_Dsecoxyz(i%dimx    , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoxyz(i%dimx    , j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
        )
        +dz3p
        *(
            dxm*(-m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+m_Dsecoz( (i-1)%dimx  , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxp*(-m_Dsecoz(i%dimx      , (j-1)%dimy, k%dimz)+m_Dsecoz(i%dimx      , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3m*(-m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+m_Dsecoxz( (i-1)%dimx  , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3p*(-m_Dsecoxz(i%dimx      , (j-1)%dimy, k%dimz)+m_Dsecoxz(i%dimx      , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxm*(-(3*dym*dym-1)*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoyz( (i-1)%dimx , j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          + dxp*(-(3*dym*dym-1)*m_Dsecoyz(i%dimx     , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoyz(i%dimx     , j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3m*(-(3*dym*dym-1)*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoxyz( (i-1)%dimx, j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3p*(-(3*dym*dym-1)*m_Dsecoxyz(i%dimx    , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoxyz(i%dimx    , j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
        )
      ;
    }

    //! return partial derivative at certain (x, y, z) for z
    double TricubicSpline::dFdz( const double x, const double y, const double z) const
    {
      const int dimx( m_Values.NumberLayers());
      const int dimy( m_Values.GetNumberRows());
      const int dimz( m_Values.GetNumberCols());

      //check if argument x is in range for non-periodic splines
      if
      (
            ( m_Border[ coord::GetAxes().e_X] != e_Periodic)
        && ( ( x < m_Start[ coord::GetAxes().e_X]) || ( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X] < x)))
      {
        BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");
        BCL_MessageDbg( "argument out of range, using linear continuation");
        if( x < m_Start[ coord::GetAxes().e_X])
        {
          return dFdz( m_Start[ coord::GetAxes().e_X], y, z);
        }
        if( x > m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X])
        {
          return dFdz( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X], y, z);
        }
      }

      //check if argument y is in range for non-periodic splines
      if
      (
            ( m_Border[ coord::GetAxes().e_Y] != e_Periodic)
        && ( ( y < m_Start[ coord::GetAxes().e_Y] || m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y] < y))
      )
      {
        BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");
        BCL_MessageDbg( "argument out of range, using linear continuation");
        if( y < m_Start[ coord::GetAxes().e_Y])
        {
          return dFdz( x, m_Start[ coord::GetAxes().e_Y], z);
        }
        if( y > m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y])
        {
          return dFdz( x, m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y], z);
        }
      }

      //check if argument z is in range for non-periodic splines
      if
      (
            ( m_Border[ coord::GetAxes().e_Z] != e_Periodic)
        && ( ( z < m_Start[ coord::GetAxes().e_Z]) || ( m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z] < z))
      )
      {
        BCL_Assert( m_LinCont[ coord::GetAxes().e_Z], "argument out of range for non-periodic spline!");
        BCL_MessageDbg( "argument out of range, using linear continuation");
        if( z < m_Start[ coord::GetAxes().e_Z])
        {
          return dFdz( x, y, m_Start[ coord::GetAxes().e_Z]);
        }
        if( z > m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z])
        {
          return dFdz( x, y, m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z]);
        }
      }

      //determine i with m_Start+(i-1)*m_Delta < x < m_Start+i*m_Delta for the correct supporting points
      int i( int( floor( ( x - m_Start[ coord::GetAxes().e_X]) / m_Delta[ coord::GetAxes().e_X])));
      while( m_Start[ coord::GetAxes().e_X] + i * m_Delta[ coord::GetAxes().e_X] < x) i++;

      //the same for j and y
      int j( int( floor( ( y - m_Start[ coord::GetAxes().e_Y]) / m_Delta[ coord::GetAxes().e_Y])));
      while( m_Start[ coord::GetAxes().e_Y] + j * m_Delta[ coord::GetAxes().e_Y] < y) j++;

      //the same for k and z
      int k( int( floor( ( z - m_Start[ coord::GetAxes().e_Z]) / m_Delta[ coord::GetAxes().e_Z])));
      while( m_Start[ coord::GetAxes().e_Z] + k * m_Delta[ coord::GetAxes().e_Z] < z) k++;

      //generate some auxiliary variables
      const double delta_aktx( x-m_Start[ coord::GetAxes().e_X]-(i-1)*m_Delta[ coord::GetAxes().e_X]);
      const double delta_akty( y-m_Start[ coord::GetAxes().e_Y]-(j-1)*m_Delta[ coord::GetAxes().e_Y]);
      const double delta_aktz( z-m_Start[ coord::GetAxes().e_Z]-(k-1)*m_Delta[ coord::GetAxes().e_Z]);

      const double dxp( delta_aktx / m_Delta[ coord::GetAxes().e_X]);
      const double dxm( 1 - dxp);
      const double dx3p( ( dxp * dxp * dxp - dxp) * Sqr( m_Delta[ coord::GetAxes().e_X]) / 6);
      const double dx3m( ( dxm * dxm * dxm - dxm) * Sqr( m_Delta[ coord::GetAxes().e_X]) / 6);

      const double dyp( delta_akty / m_Delta[ coord::GetAxes().e_Y]);
      const double dym( 1 - dyp);
      const double dy3p( ( dyp * dyp * dyp - dyp) * Sqr( m_Delta[ coord::GetAxes().e_Y]) / 6);
      const double dy3m( ( dym * dym * dym - dym) * Sqr( m_Delta[ coord::GetAxes().e_Y]) / 6);

      const double dzp( delta_aktz / m_Delta[ coord::GetAxes().e_Z]);
      const double dzm( 1 - dzp);

      //generate positive values to prevent some problems with the indizes
      while( i <= 0) i += dimx;
      while( j <= 0) j += dimy;
      while( k <= 0) k += dimz;

      return
        -(dxm*(dym*m_Values( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Values( (i-1)%dimx  , j%dimy, (k-1)%dimz))
        + dxp*(dym*m_Values(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Values(i%dimx      , j%dimy, (k-1)%dimz))
        +dx3m*(dym*m_Dsecox( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecox( (i-1)%dimx  , j%dimy, (k-1)%dimz))
        +dx3p*(dym*m_Dsecox(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecox(i%dimx      , j%dimy, (k-1)%dimz))
        + dxm*(dy3m*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoy( (i-1)%dimx , j%dimy, (k-1)%dimz))
        + dxp*(dy3m*m_Dsecoy(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoy(i%dimx     , j%dimy, (k-1)%dimz))
        +dx3m*(dy3m*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxy( (i-1)%dimx, j%dimy, (k-1)%dimz))
        +dx3p*(dy3m*m_Dsecoxy(i%dimx    , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxy(i%dimx    , j%dimy, (k-1)%dimz)))
        / m_Delta[ coord::GetAxes().e_Z]

        +(dxm*(dym*m_Values( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Values( (i-1)%dimx  , j%dimy, k%dimz))
        + dxp*(dym*m_Values(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Values(i%dimx      , j%dimy, k%dimz))
        +dx3m*(dym*m_Dsecox( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecox( (i-1)%dimx  , j%dimy, k%dimz))
        +dx3p*(dym*m_Dsecox(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecox(i%dimx      , j%dimy, k%dimz))
        + dxm*(dy3m*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoy( (i-1)%dimx , j%dimy, k%dimz))
        + dxp*(dy3m*m_Dsecoy(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoy(i%dimx     , j%dimy, k%dimz))
        +dx3m*(dy3m*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxy( (i-1)%dimx, j%dimy, k%dimz))
        +dx3p*(dy3m*m_Dsecoxy(i%dimx    , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxy(i%dimx    , j%dimy, k%dimz)))
        / m_Delta[ coord::GetAxes().e_Z]

        -(3*dzm*dzm-1)*m_Delta[ coord::GetAxes().e_Z]/ 6
        *(dxm*(dym*m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoz( (i-1)%dimx  , j%dimy, (k-1)%dimz))
        + dxp*(dym*m_Dsecoz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoz(i%dimx      , j%dimy, (k-1)%dimz))
        +dx3m*(dym*m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoxz( (i-1)%dimx  , j%dimy, (k-1)%dimz))
        +dx3p*(dym*m_Dsecoxz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoxz(i%dimx      , j%dimy, (k-1)%dimz))
        + dxm*(dy3m*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoyz( (i-1)%dimx , j%dimy, (k-1)%dimz))
        + dxp*(dy3m*m_Dsecoyz(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoyz(i%dimx     , j%dimy, (k-1)%dimz))
        +dx3m*(dy3m*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxyz( (i-1)%dimx, j%dimy, (k-1)%dimz))
        +dx3p*(dy3m*m_Dsecoxyz(i%dimx    , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxyz(i%dimx    , j%dimy, (k-1)%dimz)))

        +(3*dzp*dzp-1)*m_Delta[ coord::GetAxes().e_Z]/ 6
        *(dxm*(dym*m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecoz( (i-1)%dimx  , j%dimy, k%dimz))
        + dxp*(dym*m_Dsecoz(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecoz(i%dimx      , j%dimy, k%dimz))
        +dx3m*(dym*m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecoxz( (i-1)%dimx  , j%dimy, k%dimz))
        +dx3p*(dym*m_Dsecoxz(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecoxz(i%dimx      , j%dimy, k%dimz))
        + dxm*(dy3m*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoyz( (i-1)%dimx , j%dimy, k%dimz))
        + dxp*(dy3m*m_Dsecoyz(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoyz(i%dimx     , j%dimy, k%dimz))
        +dx3m*(dy3m*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxyz( (i-1)%dimx, j%dimy, k%dimz))
        +dx3p*(dy3m*m_Dsecoxyz(i%dimx    , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxyz(i%dimx    , j%dimy, k%dimz)))
      ;
    }

    //! return value and partial derivatives at certain (x, y, z)
    storage::Pair< double, linal::Vector< double> > TricubicSpline::FdF( const double x, const double y, const double z) const
    {
      const int dimx( m_Values.NumberLayers());
      const int dimy( m_Values.GetNumberRows());
      const int dimz( m_Values.GetNumberCols());
      double fvalue( 0), dfdxvalue( 0), dfdyvalue( 0), dfdzvalue( 0);

      //check if argument is in range for non-periodic splines
      if
      (
        (
              ( m_Border[ coord::GetAxes().e_X] != e_Periodic)
          && ( ( x < m_Start[ coord::GetAxes().e_X]) || ( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X] < x))
        )
        ||
        (
              ( m_Border[ coord::GetAxes().e_Y] != e_Periodic)
          && ( ( y < m_Start[ coord::GetAxes().e_Y]) || ( m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y] < y))
        )
        ||
        (
              ( m_Border[ coord::GetAxes().e_Z] != e_Periodic)
          && ( ( z < m_Start[ coord::GetAxes().e_Z]) || ( m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z] < z))
        )
      )
      {
        BCL_MessageDbg( "argument out of range, using linear continuation");
        if( x < m_Start[ coord::GetAxes().e_X])
        {
          BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");
          dfdxvalue = dFdx( m_Start[ coord::GetAxes().e_X], y, z);
          dfdyvalue = dFdy( m_Start[ coord::GetAxes().e_X], y, z);
          dfdzvalue = dFdz( m_Start[ coord::GetAxes().e_X], y, z);
          fvalue    = F(    m_Start[ coord::GetAxes().e_X], y, z) + ( x - m_Start[ coord::GetAxes().e_X]) * dfdxvalue;
        }
        if( x > m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X])
        {
          BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");
          dfdxvalue = dFdx( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X] , y, z);
          dfdyvalue = dFdy( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X] , y, z);
          dfdzvalue = dFdz( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X] , y, z);
          fvalue    = F(    m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X] , y, z)
                      + ( x - m_Start[ coord::GetAxes().e_X] - ( dimx - 1) * m_Delta[ coord::GetAxes().e_X]) * dfdxvalue;
        }
        if( y < m_Start[ coord::GetAxes().e_Y])
        {
          BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");
          dfdxvalue = dFdx( x, m_Start[ coord::GetAxes().e_Y], z);
          dfdyvalue = dFdy( x, m_Start[ coord::GetAxes().e_Y], z);
          dfdzvalue = dFdz( x, m_Start[ coord::GetAxes().e_Y], z);
          fvalue    = F(    x, m_Start[ coord::GetAxes().e_Y], z) + ( y - m_Start[ coord::GetAxes().e_Y]) * dfdyvalue;
        }
        if( y > m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y])
        {
          BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");
          dfdxvalue = dFdx( x, m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y], z);
          dfdyvalue = dFdy( x, m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y], z);
          dfdzvalue = dFdz( x, m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y], z);
          fvalue    = F(    x, m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y], z)
                      + ( y - m_Start[ coord::GetAxes().e_Y] - ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y]) * dfdyvalue;
        }
        if( z < m_Start[ coord::GetAxes().e_Z])
        {
          BCL_Assert( m_LinCont[ coord::GetAxes().e_Z], "argument out of range for non-periodic spline!");
          dfdxvalue = dFdx( x, y, m_Start[ coord::GetAxes().e_Z]);
          dfdyvalue = dFdy( x, y, m_Start[ coord::GetAxes().e_Z]);
          dfdzvalue = dFdz( x, y, m_Start[ coord::GetAxes().e_Z]);
          fvalue    = F(    x, y, m_Start[ coord::GetAxes().e_Z]) + ( z - m_Start[ coord::GetAxes().e_Z]) * dfdzvalue;
        }
        if( z > m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z])
        {
          BCL_Assert( m_LinCont[ coord::GetAxes().e_Z], "argument out of range for non-periodic spline!");
          dfdxvalue = dFdx( x, y, m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z]);
          dfdyvalue = dFdy( x, y, m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z]);
          dfdzvalue = dFdz( x, y, m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z]);
          fvalue    = F(    x, y, m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z])
                      + ( z - m_Start[ coord::GetAxes().e_Z] - ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z]) * dfdzvalue;
        }
      }
      else
      {
        //determine i with m_Start[ coord::GetAxes().e_X]+(i-1)*m_Delta[ coord::GetAxes().e_X] < x < m_Start[ coord::GetAxes().e_X]+i*m_Delta[ coord::GetAxes().e_X]
        //for the correct supporting points
        int i( int( floor( ( x - m_Start[ coord::GetAxes().e_X]) / m_Delta[ coord::GetAxes().e_X])));
        while( m_Start[ coord::GetAxes().e_X] + i * m_Delta[ coord::GetAxes().e_X] < x) i++;

        //determine j with m_Start[ coord::GetAxes().e_Y]+(j-1)*m_Delta[ coord::GetAxes().e_Y] < y < m_Start[ coord::GetAxes().e_Y]+j*m_Delta[ coord::GetAxes().e_Y]
        //for the correct supporting points
        int j( int( floor( ( y - m_Start[ coord::GetAxes().e_Y]) / m_Delta[ coord::GetAxes().e_Y])));
        while( m_Start[ coord::GetAxes().e_Y] + j * m_Delta[ coord::GetAxes().e_Y] < y) j++;

        //determine k with m_Start[ coord::GetAxes().e_Z]+(k-1)*m_Delta[ coord::GetAxes().e_Z] < z < m_Start[ coord::GetAxes().e_Z]+k*m_Delta[ coord::GetAxes().e_Z]
        //for the correct supporting points
        int k( int( floor( ( z - m_Start[ coord::GetAxes().e_Z]) / m_Delta[ coord::GetAxes().e_Z])));
        while( m_Start[ coord::GetAxes().e_Z] + k * m_Delta[ coord::GetAxes().e_Z] < z) k++;

        //see method F(x,y,z) for detailed formula

        const double delta_aktx( x-m_Start[ coord::GetAxes().e_X]-(i-1)*m_Delta[ coord::GetAxes().e_X]);
        const double delta_akty( y-m_Start[ coord::GetAxes().e_Y]-(j-1)*m_Delta[ coord::GetAxes().e_Y]);
        const double delta_aktz( z-m_Start[ coord::GetAxes().e_Z]-(k-1)*m_Delta[ coord::GetAxes().e_Z]);

        const double dxp( delta_aktx / m_Delta[ coord::GetAxes().e_X]);
        const double dxm( 1 - dxp);
        const double dx3p( ( dxp * dxp * dxp - dxp) * Sqr( m_Delta[ coord::GetAxes().e_X]) / 6);
        const double dx3m( ( dxm * dxm * dxm - dxm) * Sqr( m_Delta[ coord::GetAxes().e_X]) / 6);

        const double dyp( delta_akty / m_Delta[ coord::GetAxes().e_Y]);
        const double dym( 1 - dyp);
        const double dy3p( ( dyp * dyp * dyp - dyp) * Sqr( m_Delta[ coord::GetAxes().e_Y]) / 6);
        const double dy3m( ( dym * dym * dym - dym) * Sqr( m_Delta[ coord::GetAxes().e_Y]) / 6);

        const double dzp( delta_aktz / m_Delta[ coord::GetAxes().e_Z]);
        const double dzm( 1 - dzp);
        const double dz3p( ( dzp * dzp * dzp - dzp) * Sqr( m_Delta[ coord::GetAxes().e_Z]) / 6);
        const double dz3m( ( dzm * dzm * dzm - dzm) * Sqr( m_Delta[ coord::GetAxes().e_Z]) / 6);

        //generate positive values to prevent some problems with the indizes
        while( i < 1) i += dimx;
        while( j < 1) j += dimy;
        while( k < 1) k += dimz;

        //compute values
        fvalue =
        dzm
        *(dxm*(dym*m_Values( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Values( (i-1)%dimx  , j%dimy, (k-1)%dimz))
          + dxp*(dym*m_Values(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Values(i%dimx      , j%dimy, (k-1)%dimz))
        +dx3m*(dym*m_Dsecox( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecox( (i-1)%dimx  , j%dimy, (k-1)%dimz))
          +dx3p*(dym*m_Dsecox(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecox(i%dimx      , j%dimy, (k-1)%dimz))
        + dxm*(dy3m*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoy( (i-1)%dimx , j%dimy, (k-1)%dimz))
          + dxp*(dy3m*m_Dsecoy(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoy(i%dimx     , j%dimy, (k-1)%dimz))
        +dx3m*(dy3m*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxy( (i-1)%dimx, j%dimy, (k-1)%dimz))
          +dx3p*(dy3m*m_Dsecoxy(i%dimx    , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxy(i%dimx    , j%dimy, (k-1)%dimz)))

        +dzp
        *(dxm*(dym*m_Values( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Values( (i-1)%dimx  , j%dimy, k%dimz))
          + dxp*(dym*m_Values(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Values(i%dimx      , j%dimy, k%dimz))
        +dx3m*(dym*m_Dsecox( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecox( (i-1)%dimx  , j%dimy, k%dimz))
          +dx3p*(dym*m_Dsecox(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecox(i%dimx      , j%dimy, k%dimz))
        + dxm*(dy3m*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoy( (i-1)%dimx , j%dimy, k%dimz))
          + dxp*(dy3m*m_Dsecoy(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoy(i%dimx     , j%dimy, k%dimz))
        +dx3m*(dy3m*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxy( (i-1)%dimx, j%dimy, k%dimz))
          +dx3p*(dy3m*m_Dsecoxy(i%dimx    , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxy(i%dimx    , j%dimy, k%dimz)))

        +dz3m
        *(dxm*(dym*m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoz( (i-1)%dimx  , j%dimy, (k-1)%dimz))
          + dxp*(dym*m_Dsecoz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoz(i%dimx      , j%dimy, (k-1)%dimz))
        +dx3m*(dym*m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoxz( (i-1)%dimx  , j%dimy, (k-1)%dimz))
          +dx3p*(dym*m_Dsecoxz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoxz(i%dimx      , j%dimy, (k-1)%dimz))
        + dxm*(dy3m*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoyz( (i-1)%dimx , j%dimy, (k-1)%dimz))
          + dxp*(dy3m*m_Dsecoyz(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoyz(i%dimx     , j%dimy, (k-1)%dimz))
        +dx3m*(dy3m*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxyz( (i-1)%dimx, j%dimy, (k-1)%dimz))
          +dx3p*(dy3m*m_Dsecoxyz(i%dimx    , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxyz(i%dimx    , j%dimy, (k-1)%dimz)))

        +dz3p
        *(dxm*(dym*m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecoz( (i-1)%dimx  , j%dimy, k%dimz))
          + dxp*(dym*m_Dsecoz(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecoz(i%dimx      , j%dimy, k%dimz))
        +dx3m*(dym*m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecoxz( (i-1)%dimx  , j%dimy, k%dimz))
          +dx3p*(dym*m_Dsecoxz(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecoxz(i%dimx      , j%dimy, k%dimz))
        + dxm*(dy3m*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoyz( (i-1)%dimx , j%dimy, k%dimz))
          + dxp*(dy3m*m_Dsecoyz(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoyz(i%dimx     , j%dimy, k%dimz))
        +dx3m*(dy3m*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxyz( (i-1)%dimx, j%dimy, k%dimz))
          +dx3p*(dy3m*m_Dsecoxyz(i%dimx    , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxyz(i%dimx    , j%dimy, k%dimz)))
        ;

        dfdxvalue =
        dzm
        * (
          -(dym*m_Values( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Values( (i-1)%dimx  , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dym*m_Values(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Values(i%dimx      , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecox( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecox( (i-1)%dimx  , j%dimy, (k-1)%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecox(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecox(i%dimx      , j%dimy, (k-1)%dimz))
          -(dy3m*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoy( (i-1)%dimx , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dy3m*m_Dsecoy(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoy(i%dimx     , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxy( (i-1)%dimx, j%dimy, (k-1)%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxy(i%dimx    , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxy(i%dimx    , j%dimy, (k-1)%dimz))
        )

        +dzp
        *(
          -(dym*m_Values( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Values( (i-1)%dimx  , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dym*m_Values(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Values(i%dimx      , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecox( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecox( (i-1)%dimx  , j%dimy, k%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecox(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecox(i%dimx      , j%dimy, k%dimz))
          -(dy3m*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoy( (i-1)%dimx , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dy3m*m_Dsecoy(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoy(i%dimx     , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxy( (i-1)%dimx, j%dimy, k%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxy(i%dimx    , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxy(i%dimx    , j%dimy, k%dimz))
        )

        +dz3m
        *(
          -(dym*m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoz( (i-1)%dimx  , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dym*m_Dsecoz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoz(i%dimx      , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoxz( (i-1)%dimx  , j%dimy, (k-1)%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecoxz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoxz(i%dimx      , j%dimy, (k-1)%dimz))
          -(dy3m*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoyz( (i-1)%dimx , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dy3m*m_Dsecoyz(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoyz(i%dimx     , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxyz( (i-1)%dimx, j%dimy, (k-1)%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxyz(i%dimx    , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxyz(i%dimx    , j%dimy, (k-1)%dimz))
        )

        +dz3p
        *(
          -(dym*m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecoz( (i-1)%dimx  , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dym*m_Dsecoz(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecoz(i%dimx      , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecoxz( (i-1)%dimx  , j%dimy, k%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecoxz(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecoxz(i%dimx      , j%dimy, k%dimz))
          -(dy3m*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoyz( (i-1)%dimx , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dy3m*m_Dsecoyz(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoyz(i%dimx     , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxyz( (i-1)%dimx, j%dimy, k%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxyz(i%dimx    , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxyz(i%dimx    , j%dimy, k%dimz))
        )
        ;

        dfdyvalue =
        dzm
        * (
            dxm*(-m_Values( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+m_Values( (i-1)%dimx  , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxp*(-m_Values(i%dimx      , (j-1)%dimy, (k-1)%dimz)+m_Values(i%dimx      , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3m*(-m_Dsecox( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+m_Dsecox( (i-1)%dimx  , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3p*(-m_Dsecox(i%dimx      , (j-1)%dimy, (k-1)%dimz)+m_Dsecox(i%dimx      , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxm*(-(3*dym*dym-1)*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoy( (i-1)%dimx , j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          + dxp*(-(3*dym*dym-1)*m_Dsecoy(i%dimx     , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoy(i%dimx     , j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3m*(-(3*dym*dym-1)*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoxy( (i-1)%dimx, j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3p*(-(3*dym*dym-1)*m_Dsecoxy(i%dimx    , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoxy(i%dimx    , j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
        )

        +dzp
        *(
            dxm*(-m_Values( (i-1)%dimx  , (j-1)%dimy, k%dimz)+m_Values( (i-1)%dimx  , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxp*(-m_Values(i%dimx      , (j-1)%dimy, k%dimz)+m_Values(i%dimx      , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3m*(-m_Dsecox( (i-1)%dimx  , (j-1)%dimy, k%dimz)+m_Dsecox( (i-1)%dimx  , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3p*(-m_Dsecox(i%dimx      , (j-1)%dimy, k%dimz)+m_Dsecox(i%dimx      , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxm*(-(3*dym*dym-1)*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoy( (i-1)%dimx , j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          + dxp*(-(3*dym*dym-1)*m_Dsecoy(i%dimx     , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoy(i%dimx     , j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3m*(-(3*dym*dym-1)*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoxy( (i-1)%dimx, j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3p*(-(3*dym*dym-1)*m_Dsecoxy(i%dimx    , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoxy(i%dimx    , j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
        )

        +dz3m
        *(
            dxm*(-m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+m_Dsecoz( (i-1)%dimx  , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxp*(-m_Dsecoz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+m_Dsecoz(i%dimx      , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3m*(-m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+m_Dsecoxz( (i-1)%dimx  , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3p*(-m_Dsecoxz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+m_Dsecoxz(i%dimx      , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxm*(-(3*dym*dym-1)*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoyz( (i-1)%dimx , j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          + dxp*(-(3*dym*dym-1)*m_Dsecoyz(i%dimx     , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoyz(i%dimx     , j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3m*(-(3*dym*dym-1)*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoxyz( (i-1)%dimx, j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3p*(-(3*dym*dym-1)*m_Dsecoxyz(i%dimx    , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoxyz(i%dimx    , j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
        )

        +dz3p
        *(
            dxm*(-m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+m_Dsecoz( (i-1)%dimx  , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxp*(-m_Dsecoz(i%dimx      , (j-1)%dimy, k%dimz)+m_Dsecoz(i%dimx      , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3m*(-m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+m_Dsecoxz( (i-1)%dimx  , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3p*(-m_Dsecoxz(i%dimx      , (j-1)%dimy, k%dimz)+m_Dsecoxz(i%dimx      , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxm*(-(3*dym*dym-1)*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoyz( (i-1)%dimx , j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          + dxp*(-(3*dym*dym-1)*m_Dsecoyz(i%dimx     , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoyz(i%dimx     , j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3m*(-(3*dym*dym-1)*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoxyz( (i-1)%dimx, j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3p*(-(3*dym*dym-1)*m_Dsecoxyz(i%dimx    , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoxyz(i%dimx    , j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
        )
        ;

        dfdzvalue =
        -(dxm*(dym*m_Values( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Values( (i-1)%dimx  , j%dimy, (k-1)%dimz))
          + dxp*(dym*m_Values(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Values(i%dimx      , j%dimy, (k-1)%dimz))
        +dx3m*(dym*m_Dsecox( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecox( (i-1)%dimx  , j%dimy, (k-1)%dimz))
          +dx3p*(dym*m_Dsecox(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecox(i%dimx      , j%dimy, (k-1)%dimz))
        + dxm*(dy3m*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoy( (i-1)%dimx , j%dimy, (k-1)%dimz))
          + dxp*(dy3m*m_Dsecoy(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoy(i%dimx     , j%dimy, (k-1)%dimz))
        +dx3m*(dy3m*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxy( (i-1)%dimx, j%dimy, (k-1)%dimz))
          +dx3p*(dy3m*m_Dsecoxy(i%dimx    , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxy(i%dimx    , j%dimy, (k-1)%dimz)))
        / m_Delta[ coord::GetAxes().e_Z]

        +(dxm*(dym*m_Values( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Values( (i-1)%dimx  , j%dimy, k%dimz))
          + dxp*(dym*m_Values(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Values(i%dimx      , j%dimy, k%dimz))
        +dx3m*(dym*m_Dsecox( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecox( (i-1)%dimx  , j%dimy, k%dimz))
          +dx3p*(dym*m_Dsecox(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecox(i%dimx      , j%dimy, k%dimz))
        + dxm*(dy3m*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoy( (i-1)%dimx , j%dimy, k%dimz))
          + dxp*(dy3m*m_Dsecoy(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoy(i%dimx     , j%dimy, k%dimz))
        +dx3m*(dy3m*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxy( (i-1)%dimx, j%dimy, k%dimz))
          +dx3p*(dy3m*m_Dsecoxy(i%dimx    , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxy(i%dimx    , j%dimy, k%dimz)))
        / m_Delta[ coord::GetAxes().e_Z]

        -(3*dzm*dzm-1)*m_Delta[ coord::GetAxes().e_Z]/ 6
        *(dxm*(dym*m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoz( (i-1)%dimx  , j%dimy, (k-1)%dimz))
          + dxp*(dym*m_Dsecoz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoz(i%dimx      , j%dimy, (k-1)%dimz))
        +dx3m*(dym*m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoxz( (i-1)%dimx  , j%dimy, (k-1)%dimz))
          +dx3p*(dym*m_Dsecoxz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoxz(i%dimx      , j%dimy, (k-1)%dimz))
        + dxm*(dy3m*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoyz( (i-1)%dimx , j%dimy, (k-1)%dimz))
          + dxp*(dy3m*m_Dsecoyz(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoyz(i%dimx     , j%dimy, (k-1)%dimz))
        +dx3m*(dy3m*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxyz( (i-1)%dimx, j%dimy, (k-1)%dimz))
          +dx3p*(dy3m*m_Dsecoxyz(i%dimx    , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxyz(i%dimx    , j%dimy, (k-1)%dimz)))

        +(3*dzp*dzp-1)*m_Delta[ coord::GetAxes().e_Z]/ 6
        *(dxm*(dym*m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecoz( (i-1)%dimx  , j%dimy, k%dimz))
          + dxp*(dym*m_Dsecoz(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecoz(i%dimx      , j%dimy, k%dimz))
        +dx3m*(dym*m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecoxz( (i-1)%dimx  , j%dimy, k%dimz))
          +dx3p*(dym*m_Dsecoxz(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecoxz(i%dimx      , j%dimy, k%dimz))
        + dxm*(dy3m*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoyz( (i-1)%dimx , j%dimy, k%dimz))
          + dxp*(dy3m*m_Dsecoyz(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoyz(i%dimx     , j%dimy, k%dimz))
        +dx3m*(dy3m*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxyz( (i-1)%dimx, j%dimy, k%dimz))
          +dx3p*(dy3m*m_Dsecoxyz(i%dimx    , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxyz(i%dimx    , j%dimy, k%dimz)))
        ;
      }

      return storage::Pair< double, linal::Vector< double> >( fvalue, linal::MakeVector< double>( dfdxvalue, dfdyvalue, dfdzvalue));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! read TricubicSpline from std::istream
    std::istream &TricubicSpline::Read( std::istream &ISTREAM)
    {
      // read parameters
      int border;
      ISTREAM >> border;
      m_Border[ coord::GetAxes().e_X] = SplineBorderType( border);
      ISTREAM >> border;
      m_Border[ coord::GetAxes().e_Y] = SplineBorderType( border);
      ISTREAM >> border;
      m_Border[ coord::GetAxes().e_Z] = SplineBorderType( border);
      ISTREAM >> m_Start[ coord::GetAxes().e_X];
      ISTREAM >> m_Start[ coord::GetAxes().e_Y];
      ISTREAM >> m_Start[ coord::GetAxes().e_Z];
      ISTREAM >> m_Delta[ coord::GetAxes().e_X];
      ISTREAM >> m_Delta[ coord::GetAxes().e_Y];
      ISTREAM >> m_Delta[ coord::GetAxes().e_Z];
      ISTREAM >> m_Values;
      ISTREAM >> m_Dsecox;
      ISTREAM >> m_Dsecoy;
      ISTREAM >> m_Dsecoxy;
      ISTREAM >> m_Dsecoz;
      ISTREAM >> m_Dsecoxz;
      ISTREAM >> m_Dsecoyz;
      ISTREAM >> m_Dsecoxyz;
      ISTREAM >> m_LinCont[ coord::GetAxes().e_X];
      ISTREAM >> m_LinCont[ coord::GetAxes().e_Y];
      ISTREAM >> m_LinCont[ coord::GetAxes().e_Z];
      ISTREAM >> m_FirstBe[ coord::GetAxes().e_X].First();
      ISTREAM >> m_FirstBe[ coord::GetAxes().e_X].Second();
      ISTREAM >> m_FirstBe[ coord::GetAxes().e_Y].First();
      ISTREAM >> m_FirstBe[ coord::GetAxes().e_Y].Second();
      ISTREAM >> m_FirstBe[ coord::GetAxes().e_Z].First();
      ISTREAM >> m_FirstBe[ coord::GetAxes().e_Z].Second();

      // end
      return ISTREAM;
    }

    //! write TricubicSpline into std::ostream
    std::ostream &TricubicSpline::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write
      OSTREAM << m_Border[ coord::GetAxes().e_X]           << '\n';
      OSTREAM << m_Border[ coord::GetAxes().e_Y]           << '\n';
      OSTREAM << m_Border[ coord::GetAxes().e_Z]           << '\n';
      OSTREAM << m_Start[ coord::GetAxes().e_X]            << '\n';
      OSTREAM << m_Start[ coord::GetAxes().e_Y]            << '\n';
      OSTREAM << m_Start[ coord::GetAxes().e_Z]            << '\n';
      OSTREAM << m_Delta[ coord::GetAxes().e_X]            << '\n';
      OSTREAM << m_Delta[ coord::GetAxes().e_Y]            << '\n';
      OSTREAM << m_Delta[ coord::GetAxes().e_Z]            << '\n';
      OSTREAM << m_Values                   << '\n';
      OSTREAM << m_Dsecox                   << '\n';
      OSTREAM << m_Dsecoy                   << '\n';
      OSTREAM << m_Dsecoxy                  << '\n';
      OSTREAM << m_Dsecoz                   << '\n';
      OSTREAM << m_Dsecoxz                  << '\n';
      OSTREAM << m_Dsecoyz                  << '\n';
      OSTREAM << m_Dsecoxyz                 << '\n';
      OSTREAM << m_LinCont[ coord::GetAxes().e_X]          << '\n';
      OSTREAM << m_LinCont[ coord::GetAxes().e_Y]          << '\n';
      OSTREAM << m_LinCont[ coord::GetAxes().e_Z]          << '\n';
      OSTREAM << m_FirstBe[ coord::GetAxes().e_X].First()  << '\n';
      OSTREAM << m_FirstBe[ coord::GetAxes().e_X].Second() << '\n';
      OSTREAM << m_FirstBe[ coord::GetAxes().e_Y].First()  << '\n';
      OSTREAM << m_FirstBe[ coord::GetAxes().e_Y].Second() << '\n';
      OSTREAM << m_FirstBe[ coord::GetAxes().e_Z].First()  << '\n';
      OSTREAM << m_FirstBe[ coord::GetAxes().e_Z].Second() << '\n';

      // end
      return OSTREAM;
    }

  } // namespace math
} // namespace bcl
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
#include "math/bcl_math_z_score.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ZScore::s_Instance
    (
      GetObjectInstances().AddInstance( new ZScore())
    );

  //////////
  // data //
  //////////

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &ZScore::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "z_score");

      return s_default_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ZScore::ZScore() :
      m_Mean( 0),
      m_Sigma( 0)
    {
    }

    //! @brief construct with Mean and Standard Deviation
    //! @param MEAN mean value of the considered distribution
    //! @param SIGMA standard deviation of the considered distribution
    //! @param SCHEME
    ZScore::ZScore
    (
      const double MEAN,
      const double SIGMA,
      const std::string &SCHEME
    ) :
      m_Mean( MEAN),
      m_Sigma( SIGMA),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new ZScore
    ZScore *ZScore::Clone() const
    {
      return new ZScore( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ZScore::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief access to Mean
    //! @return the mean used in the z-score
    double ZScore::GetMean() const
    {
      return m_Mean;
    };

    //! @brief access to Standard Deviation
    //! @return the sigma used in the z-score
    double ZScore::GetSigma() const
    {
      return m_Sigma;
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &ZScore::GetScheme() const
    {
      return m_Scheme;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculates the z-score of the value in the normal distribution with the given mean and sigma
    //! @param ARGUMENT value which z-score shall be calculated from
    //! @return z-score of the given value
    double ZScore::operator()( const double &ARGUMENT) const
    {
      return m_Sigma ? ( ARGUMENT - m_Mean) / m_Sigma : ARGUMENT - m_Mean;
    }

    //! @brief multiply z-score by given weight
    //! @param WEIGHT weight to multiply by
    //! @return reference to this
    ZScore &ZScore::operator *=( const double WEIGHT)
    {
      m_Sigma /= WEIGHT;
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ZScore::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Mean  , ISTREAM);
      io::Serialize::Read( m_Sigma , ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &ZScore::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Mean  , OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_Sigma , OSTREAM,      0) << '\t';
      io::Serialize::Write( m_Scheme, OSTREAM,      0);

      // end
      return OSTREAM;
    }

  } // namespace math
} // namespace bcl
