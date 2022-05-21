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
