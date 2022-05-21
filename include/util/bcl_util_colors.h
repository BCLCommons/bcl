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

#ifndef BCL_UTIL_COLORS_H_
#define BCL_UTIL_COLORS_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_util_enumerate.h"
#include "linal/bcl_linal_vector_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Colors
    //! @brief enumeration of classes as defined by their RGB values
    //! @details Each color is represented with an RGB value, which is a vector of 3 double each with range [0..1]
    //! where the first number corresponds to red, second number corresponds to blue and the third number corresponds
    //! to blue contribution.
    //!
    //! @see @link example_util_colors.cpp @endlink
    //! @author karakam
    //! @date Nov 6, 2010
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Colors :
      public Enumerate< linal::Vector3D, Colors>
    {

      friend class Enumerate< linal::Vector3D, Colors>;

    public:

    //////////
    // data //
    //////////

      // declare commonly used colors
      const Color e_White;   //!< white
      const Color e_Black;   //!< black
      const Color e_Red;     //!< red
      const Color e_Green;   //!< green
      const Color e_Blue;    //!< blue
      const Color e_Yellow;  //!< yellow
      const Color e_Cyan;    //!< cyan
      const Color e_Magenta; //!< magenta

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Colors();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief return vector of rainbow colors
      //! @return vector of rainbow
      const storage::Vector< Color> &GetRainbow() const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief convert an rgd color to a hsv color
      //! @param RGB_COLOR the color as rgb values 0-1
      //! @return hsv color
      static linal::Vector3D ConvertRGBToHSV( const linal::Vector3D &RGB_COLOR);

    }; // class Colors

    //! @brief construct on access function for all Colors
    //! @return reference to only instances of Colors
    BCL_API
    const Colors &GetColors();

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< linal::Vector3D, Colors>;

  } // namespace util
} // namespace bcl

#endif // BCL_UTIL_COLORS_H_ 
