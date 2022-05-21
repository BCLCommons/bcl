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

#ifndef BCL_MATH_SMOOTH_DATA_H_
#define BCL_MATH_SMOOTH_DATA_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SmoothData
    //! @brief for a collection of data smoothing functions
    //! @details this class provides a set of functions for smoothing data presented in: Vector, Matrix, Tensor
    //!
    //! @see @link example_math_smooth_data.cpp @endlink
    //! @author mueller
    //! @date Nov 10, 2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SmoothData
    {
    public:

    //////////
    // data //
    //////////

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    private:

      //! @brief default constructor
      SmoothData()
      {
      }

    public:

      //! @brief a linear smoothing
      //! @param DATA a Vector of T containg the data to be smoothed
      //! @param FRACTION a value between 0 and 1 - how much of the orginal data should be used vs the fraction of the neigbors
      //! @param SMOOTH_BORDERS bool on whether to smooth borders or not
      template< typename t_DataType>
      static linal::Vector< t_DataType> SmoothVector( const linal::VectorConstInterface< t_DataType> &DATA, const double FRACTION, const bool SMOOTH_BORDERS = false)
      {
        BCL_Assert( 0.0 <= FRACTION && FRACTION <= 1.0, "FRACTION must be between 0. and 1.");

        //for datasets smaller than two, no smoothing
        if( DATA.GetSize() < 2)
        {
          return DATA;
        }

        const double neighborpairweight( ( 1 - FRACTION) / 2.0);
        const double neighborweight( 1 - FRACTION);

        //make a copy of the DATA
        linal::Vector< t_DataType> smoothed_data( DATA);

        //smooth the borders
        if( SMOOTH_BORDERS)
        {
          smoothed_data( 0) = FRACTION * DATA.First() + neighborweight * DATA( 1);
          smoothed_data( DATA.GetSize() - 1) = FRACTION * DATA.Last() + neighborweight * DATA( DATA.GetSize() - 2);
        }

        const t_DataType *data_ptr( DATA.Begin() + 1), *data_ptr_end( DATA.End() - 1);

        //iterate over reaming values - take an average between the left and the right one and combine wight by FRACTION
        for
        (
          t_DataType *smooth_ptr( smoothed_data.Begin() + 1);
          data_ptr != data_ptr_end;
          ++data_ptr, ++smooth_ptr
        )
        {
          *smooth_ptr = FRACTION * ( *data_ptr) + neighborpairweight * ( ( *( data_ptr - 1)) + ( *( data_ptr + 1)));
        }

        //end
        return smoothed_data;
      }

      template< typename t_DataType>
      static linal::Matrix< t_DataType> SmoothMatrix
      (
        const linal::MatrixConstInterface< t_DataType> &DATA,
        const double FRACTION,
        const bool SMOOTH_BORDERS = false
      )
      {
        BCL_Assert( 0.0 <= FRACTION && FRACTION <= 1.0, "FRACTION must be between 0. and 1.");

        //for datasets smaller than two, no smoothing
        if( DATA.GetNumberRows() < 2 && DATA.GetNumberCols() < 2)
        {
          return DATA;
        }

        //preprocessing of data
        linal::Matrix< t_DataType> smoothed_data( DATA);

        // smooth borders
        if( SMOOTH_BORDERS)
        {
          //corners
          smoothed_data( 0, 0) = FRACTION * DATA( 0, 0) + ( 1 - FRACTION) * ( DATA( 0, 1) + DATA( 1, 0) + DATA( 1, 1)) / 3;
          smoothed_data( DATA.GetNumberRows() - 1, 0) = FRACTION * DATA( DATA.GetNumberRows() - 1, 0) + ( 1 - FRACTION) * ( DATA( DATA.GetNumberRows() - 2, 0) + DATA( DATA.GetNumberRows() - 2, 1) + DATA( DATA.GetNumberRows() - 1, 1)) / 3;
          smoothed_data( 0, DATA.GetNumberCols() - 1) = FRACTION * DATA( 0, DATA.GetNumberCols() - 1) + ( 1 - FRACTION) * ( DATA( 0, DATA.GetNumberCols() - 2) + DATA( 1, DATA.GetNumberCols() - 2) + DATA( 1, DATA.GetNumberCols() - 1)) / 3;
          smoothed_data( DATA.GetNumberRows() - 1, DATA.GetNumberCols() - 1) = FRACTION * DATA( DATA.GetNumberRows() - 1, DATA.GetNumberCols() - 1) + ( 1 - FRACTION) * ( DATA( DATA.GetNumberRows() - 2, DATA.GetNumberCols() - 2) + DATA( DATA.GetNumberRows() - 1, DATA.GetNumberCols() - 2) + DATA( DATA.GetNumberRows() - 2, DATA.GetNumberCols() - 1)) / 3;

          //sides
          for( size_t r( 1); r < DATA.GetNumberRows() - 1; ++r)
          {
            smoothed_data( r, 0) = FRACTION * DATA( r, 0) + ( 1 - FRACTION) * ( DATA( r - 1, 0) + DATA( r + 1, 0) + DATA( r - 1, 1) + DATA( r, 1) + DATA( r + 1, 1)) / 5;
            smoothed_data( r, DATA.GetNumberCols() - 1) = FRACTION * DATA( r, DATA.GetNumberCols() - 1) + ( 1 - FRACTION) * ( DATA( r - 1, DATA.GetNumberCols() - 1) + DATA( r + 1, DATA.GetNumberCols() - 1) + DATA( r - 1, DATA.GetNumberCols() - 2) + DATA( r, DATA.GetNumberCols() - 2) + DATA( r + 1, DATA.GetNumberCols() - 2)) / 5;
          }
          for( size_t c( 1); c < DATA.GetNumberCols() - 1; ++c)
          {
            smoothed_data(0,c) = FRACTION*DATA(0,c)+(1-FRACTION)*(DATA(0,c-1)+DATA(0,c+1)+DATA(1,c-1)+DATA(1,c)+DATA(1,c+1))/5;
            smoothed_data(DATA.GetNumberRows()-1,c) = FRACTION*DATA(DATA.GetNumberRows()-1,c)+(1-FRACTION)*(DATA(DATA.GetNumberRows()-1,c-1)+DATA(DATA.GetNumberRows()-1,c+1)+DATA(DATA.GetNumberRows()-2, c-1)+DATA(DATA.GetNumberRows()-2, c)+DATA(DATA.GetNumberRows()-2, c+1))/5;
          }
        }

        //the middle elements
        for( size_t r( 1); r < DATA.GetNumberRows() - 1; ++r)
        {
          for( size_t c( 1); c < DATA.GetNumberCols() - 1; ++c)
          {
            smoothed_data( r, c) = FRACTION * DATA( r, c) + ( 1 - FRACTION) * ( DATA( r - 1, c - 1) + DATA( r - 1, c) + DATA( r - 1, c + 1) + DATA( r, c - 1) + DATA( r, c + 1) + DATA( r + 1, c - 1) + DATA( r + 1, c) + DATA( r + 1, c + 1)) / 8;
          }
        }

        //end
        return smoothed_data;
      }

    }; // class SmoothData

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_SMOOTH_DATA_H_
