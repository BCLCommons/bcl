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
#include "linal/bcl_linal_operations_interface.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix.h"
#include "util/bcl_util_message.h"
#include "util/bcl_util_si_ptr_vector.h"

using bcl::command::Parameter;
using bcl::command::FlagStatic;

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief creates a matrix from a SiPtrVector of VectorInterfaces
    //! @param LIST_VECTORS the SiPtrVector of VectorInterfaces to be transferred into a matrix
    //! @return the matrix created from the SiPtrList of VectorInterfaces
    template< typename t_DataType>
    Matrix< t_DataType> OperationsInterface< t_DataType>::CopyListToMatrix
    (
      const util::SiPtrVector< const VectorConstInterface< t_DataType> > &LIST_VECTORS
    )
    {
      // check if list is empty
      if( LIST_VECTORS.IsEmpty())
      {
        return Matrix< t_DataType>();
      }

      // getting matrix dimensions
      const size_t num_rows( LIST_VECTORS.GetSize());
      const size_t num_cols( LIST_VECTORS.FirstElement()->GetSize());

      // allocate matrix
      Matrix< t_DataType> matrix( num_rows, num_cols, t_DataType( 0));
      t_DataType *row_itr( matrix.Begin());

      // iterate over list and rows
      for
      (
          typename util::SiPtrVector< const VectorConstInterface< t_DataType> >::const_iterator
          list_itr( LIST_VECTORS.Begin()),
          list_itr_end( LIST_VECTORS.End());
          list_itr != list_itr_end;
          ++list_itr, row_itr += num_cols
      )
      {
        // get current vector and verify size consistency
        const VectorConstInterface< t_DataType> &current_vector( **list_itr);
        if( current_vector.GetSize() != num_cols)
        {
          BCL_MessageCrt( "Supplied list of vectors have different lengths!");
          return Matrix< t_DataType>();
        }

        // copy vector into row
        std::copy( current_vector.Begin(), current_vector.End(), row_itr);
      }

      // end
      return matrix;
    }

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API OperationsInterface< float>;
    template class BCL_API OperationsInterface< double>;

  } // namespace linal
} // namespace bcl
