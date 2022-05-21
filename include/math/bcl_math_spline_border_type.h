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

#ifndef BCL_MATH_SPLINE_BORDER_TYPE_H_
#define BCL_MATH_SPLINE_BORDER_TYPE_H_

// include the namespace header
#include "bcl_math.h"

namespace bcl
{
  namespace math
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @file bcl_math_spline_border_type.h
    //! @brief enum for different border types (boundary conditions) for cubic splines
    //!
    //! @remarks example unnecessary
    //! @author mueller
    //! @date Jul 22, 2010
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Boundary condition for a spline
    enum SplineBorderType
    {
      e_Natural,   //!< 2nd derivative = 0 at endpoints
      e_Periodic,  //!< 1st and 2nd derivatives equal at endpoints
      e_FirstDer,  //!< First derivatives constant at endpoints (equal to given slope)
      e_NotAKnot   //!< Requires that y''' be continuous across the second and penultimate points
    };

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_SPLINE_BORDER_TYPE_H_

