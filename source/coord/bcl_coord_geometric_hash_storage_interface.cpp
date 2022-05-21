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
#include "coord/bcl_coord_geometric_hash_storage_interface.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_transformation_matrix_3d.h"
#include "storage/bcl_storage_triplet.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

    //! This function converts triplet of ints to a hashkey by shifting each integer by certain amount, so that in the binary representation these three integers do not overlap, given that the ints are smaller than pow(2,s_SingleKeyShift)-1
    size_t GeometricHashStorageInterface::ConvertTripletToKey( const storage::Triplet< int, int, int> &HASHKEYS)
    {
      return ( size_t( HASHKEYS.First()) << ( 2 * s_SingleKeyShift)) ^
             ( size_t( HASHKEYS.Second()) <<       s_SingleKeyShift)  ^
               size_t( HASHKEYS.Third());
    }

    //! checks if the TRANSFORMATIONMATRIX3D is similar to any stored transformation matrix according to DIFFERENCE
    bool GeometricHashStorageInterface::IsSimilarTransformation
    (
      const util::ShPtrVector< math::TransformationMatrix3D> &TRANSFORAMTIONS,
      const math::TransformationMatrix3D &TRANSFORMATIONMATRIX3D,
      const storage::VectorND< 2, double> &DIFFERENCE_ROT_TRANS
    )
    {
      // iterate over all stored transformation matrices
      for
      (
        util::ShPtrVector< math::TransformationMatrix3D>::const_iterator
          matrix_itr( TRANSFORAMTIONS.Begin()), matrix_itr_end( TRANSFORAMTIONS.End());
        matrix_itr != matrix_itr_end;
        ++matrix_itr
      )
      {
        // check if current matrix is similar within given tolerance to given TRANSFORMATIONMATRIX3D
        if
        (
          math::SimilarWithinTolerance
          (
            **matrix_itr,
            TRANSFORMATIONMATRIX3D,
            DIFFERENCE_ROT_TRANS.Second(),
            DIFFERENCE_ROT_TRANS.First()
          )
        )
        {
          return true;
        }
      }

      // no similar transformation matrix found
      return false;
    }

  } // namespace coord
} // namespace bcl
