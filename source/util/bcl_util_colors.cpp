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
#include "util/bcl_util_colors.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_range.h"
#include "math/bcl_math_running_min_max.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct all AAClasses
    Colors::Colors() :
      e_White(   AddEnum( "White",   linal::Vector3D( 1.00, 1.00, 1.00))),
      e_Black(   AddEnum( "Black",   linal::Vector3D( 0.00, 0.00, 0.00))),
      e_Red(     AddEnum( "Red",     linal::Vector3D( 1.00, 0.00, 0.00))),
      e_Green(   AddEnum( "Green",   linal::Vector3D( 0.00, 1.00, 0.00))),
      e_Blue(    AddEnum( "Blue",    linal::Vector3D( 0.00, 0.00, 1.00))),
      e_Yellow(  AddEnum( "Yellow",  linal::Vector3D( 1.00, 1.00, 0.00))),
      e_Cyan(    AddEnum( "Cyan",    linal::Vector3D( 0.00, 1.00, 1.00))),
      e_Magenta( AddEnum( "Magenta", linal::Vector3D( 1.00, 0.00, 1.00)))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Colors::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief return vector of rainbow colors
    //! @return vector of rainbow
    const storage::Vector< Color> &Colors::GetRainbow() const
    {
      // initialize static vector
      static const storage::Vector< Color> s_rainbow
      (
        storage::Vector< Color>::Create( e_Blue, e_Green, e_Yellow, e_Red)
      );

      // end
      return s_rainbow;
    }

    //! @brief construct on access function for all Colors
    //! @return reference to only instances of Colors
    const Colors &GetColors()
    {
      return Colors::GetEnums();
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief convert an rgd color to a hsv color
    //! @param RGB_COLOR the color as rgb values 0-1
    //! @return hsv color
    linal::Vector3D Colors::ConvertRGBToHSV( const linal::Vector3D &RGB_COLOR)
    {
      math::RunningMinMax< double> min_max;
      min_max += RGB_COLOR.X();
      min_max += RGB_COLOR.Y();
      min_max += RGB_COLOR.Z();

      const math::Range< double> range( min_max.GetMin(), min_max.GetMax());

      double s(  0);
      double h( -1);

      // s
      if( range.GetMax() != 0)
      {
        s = range.GetWidth() / range.GetMax();
      }

      if( RGB_COLOR.X() == range.GetMax())
      {
        h = ( RGB_COLOR.Y() - RGB_COLOR.Z() ) / range.GetWidth();   // between yellow & magenta
      }
      else if( RGB_COLOR.Y() == range.GetMax())
      {
        h = 2 + ( RGB_COLOR.Z() - RGB_COLOR.X()) / range.GetWidth(); // between cyan & yellow
      }
      else
      {
        h = 4 + ( RGB_COLOR.X() - RGB_COLOR.Y()) / range.GetWidth(); // between magenta & cyan
      }

      h /= 6.0; // normalize to [0,1]
      if( h < 0)
      {
        h += 1.0;
      }

      // end
      return linal::Vector3D( h, s, range.GetMax());
    }

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< linal::Vector3D, Colors>;

  } // namespace util
} // namespace bcl
