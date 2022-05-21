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
#include "coord/bcl_coord_moment_of_inertia.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "math/bcl_math_statistics.h"
#include "math/bcl_math_transformation_matrix_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MomentOfInertia::s_Instance
    (
      GetObjectInstances().AddInstance( new MomentOfInertia())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MomentOfInertia::MomentOfInertia()
    {
    }

    //! @brief Clone function
    //! @return pointer to new MomentOfInertia
    MomentOfInertia *MomentOfInertia::Clone() const
    {
      return new MomentOfInertia( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MomentOfInertia::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate transformation that translates into the center of weights and that sorts principal axes of inertia according to principal moments of inertia x - smallest, z - largest
    //! @param COORDINATES_WEIGHT_MATRIX matrix with coordinate and weight (4 cols) in number points rows
    //! @return transformation matrix and principal moments of inertias
    storage::Pair< math::TransformationMatrix3D, linal::Vector3D>
    MomentOfInertia::TransformationAndMoments
    (
      const linal::MatrixConstInterface< double> &COORDINATES_WEIGHT_MATRIX
    ) const
    {
      BCL_Assert( COORDINATES_WEIGHT_MATRIX.GetNumberCols() == 4, "need 4 columns");

      const storage::Pair< math::TransformationMatrix3D, linal::Vector3D> error_return
      (
        math::TransformationMatrix3D( util::UndefinedObject()),
        linal::Vector3D( util::GetUndefined< double>())
      );

      // center of mass
      const linal::Vector3D center_of_mass( CenterOfMass( COORDINATES_WEIGHT_MATRIX));

      // error
      if( !center_of_mass.IsDefined())
      {
        return error_return;
      }

      double ixx( 0);
      double iyy( 0);
      double izz( 0);
      double ixy( 0);
      double ixz( 0);
      double iyz( 0);

      // iterate over all amino acids
      for
      (
        const double *row( COORDINATES_WEIGHT_MATRIX.Begin()), *row_end( COORDINATES_WEIGHT_MATRIX.End());
        row != row_end;
        row += 4
      )
      {
        linal::Vector3D current_coord( row);
        current_coord -= center_of_mass;

        const double current_weight( row[ 3]);

        // diagonal
        ixx += current_weight * ( math::Sqr( current_coord.Y()) + math::Sqr( current_coord.Z()));
        iyy += current_weight * ( math::Sqr( current_coord.X()) + math::Sqr( current_coord.Z()));
        izz += current_weight * ( math::Sqr( current_coord.X()) + math::Sqr( current_coord.Y()));

        // off diagonal; products of inertia
        ixy -= current_weight * current_coord.X() * current_coord.Y();
        ixz -= current_weight * current_coord.X() * current_coord.Z();
        iyz -= current_weight * current_coord.Y() * current_coord.Z();
      }

      // return the tensor
      const double elements[ 3 * 3] = { ixx, ixy, ixz, ixy, iyy, iyz, ixz, iyz, izz};
      linal::Matrix3x3< double> moment( elements);

      linal::Matrix3x3< double> eigen_vectors( 0.0);
      linal::Vector< double> eigenvalues( 3);
      moment.EigenVectorsSymmetric( eigen_vectors, eigenvalues);
      // transpose so that the eigenvectors are in the rows
      eigen_vectors.Transpose();
      // sort eigenvectors by eigenvalues
      eigen_vectors.SortRowsAndVector( eigenvalues);

      eigen_vectors.SwapRows( 0, 2);
      std::swap( eigenvalues( 0), eigenvalues( 2));

      // orthogonalize
      eigen_vectors.ReplaceRow
      (
        2,
        linal::CrossProduct( linal::Vector3D( eigen_vectors[ 0]), linal::Vector3D( eigen_vectors[ 1]))
      );

      // shift and rotate molecule
      math::TransformationMatrix3D transformation;
      transformation( math::RotationMatrix3D( eigen_vectors));
      transformation( center_of_mass);

      // end
      return
        storage::Pair< math::TransformationMatrix3D, linal::Vector3D>
        (
          math::Inverse( transformation),
          linal::Vector3D( eigenvalues)
        );
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that returns the ordered principal moments of inertia
    //! @param COORDINATE_WEIGHT_MATRIX 4*n matrix with three coordinates and one weight as ciolumns, and n rows
    //! @return ordered principal moments of inertia
    linal::Vector3D MomentOfInertia::operator()( const linal::MatrixConstInterface< double> &COORDINATE_WEIGHT_MATRIX) const
    {
      return TransformationAndMoments( COORDINATE_WEIGHT_MATRIX).Second();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MomentOfInertia::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MomentOfInertia::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculate the center of mass
    //! @param COORDINATES_WEIGHT_MATRIX matrix with coordinate and weight (4 cols) in number points rows
    //! @return the center of mass
    linal::Vector3D
    MomentOfInertia::CenterOfMass( const linal::MatrixConstInterface< double> &COORDINATES_WEIGHT_MATRIX)
    {
      BCL_Assert( COORDINATES_WEIGHT_MATRIX.GetNumberCols() == 4, "need 4 columns: 3 coordinates and 1 weight");

      // sum of weight
      double weight_sum( 0);

      // sum of weighted positions
      linal::Vector3D weighted_positions( 0.0);

      // iterate over all rows
      for
      (
        const double *row( COORDINATES_WEIGHT_MATRIX.Begin()), *row_end( COORDINATES_WEIGHT_MATRIX.End());
        row != row_end;
        row += 4
      )
      {
        // current coordinate
        linal::Vector3D current_coord( row);

        // current weight
        const double current_weight( row[ 3]);

        // weigh position
        current_coord *= current_weight;

        // add position
        weighted_positions += current_coord;

        // sum weights
        weight_sum += current_weight;
      }

      // average
      return weighted_positions / weight_sum;
    }

  } // namespace coord
} // namespace bcl
