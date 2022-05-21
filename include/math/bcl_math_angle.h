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

#ifndef BCL_MATH_ANGLE_H_
#define BCL_MATH_ANGLE_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Angle
    //! @brief a class bundling all functionalities related to angles
    //! @details characters and strings to represent the degree symbol; function to convert between degree and radian
    //!
    //! @see @link example_math_angle.cpp @endlink
    //! @author woetzen
    //! @date Aug 4, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Angle
    {

    public:

    //////////
    // data //
    //////////

      //! ASCII code for degree symbol
      static const unsigned char s_DegreeChar = 248;

      //! string used to represent the degree character in gnuplot scripts
      static const std::string s_DegreeSymbolGnuplot;

      enum Unit
      {
        e_Radian,
        e_Degree,
        e_Undefined,
        s_NumberUnitTypes
      };

      //! @brief conversion to a string from a Unit
      //! @param UNIT the Unit to get a string for
      //! @return a string representing that Unit
      static const std::string &GetUnitName( const Unit &UNIT);

      //! @brief enum class wrapper for Unit
      typedef util::WrapperEnum< Unit, &GetUnitName, s_NumberUnitTypes> UnitEnum;

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor not implemented, since there are just static functions in that class
      //! @brief the prevents others from creating an object of this class
      Angle();

    public:

    ////////////////
    // operations //
    ////////////////

      //! @brief convert an angle that is given in degree into radians
      //! @param DEGREE that value in degree
      //! @return the value in radian
      template< typename t_DataType>
      static t_DataType Radian( const t_DataType &DEGREE)
      {
        return t_DataType( DEGREE * g_Pi / t_DataType( 180));
      }

      //! @brief convert an angle that is given in radian into degree
      //! @param RADIAN that value in radian
      //! @return the value in degree
      template< typename t_DataType>
      static t_DataType Degree( const t_DataType &RADIAN)
      {
        return t_DataType( RADIAN / g_Pi * t_DataType( 180));
      }

    }; // class Angle

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_ANGLE_H_
