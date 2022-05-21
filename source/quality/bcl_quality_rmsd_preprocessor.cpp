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
#include "quality/bcl_quality_rmsd_preprocessor.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_operations.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "linal/bcl_linal_vector_operations.h"
#include "linal/bcl_linal_vector_reference.h"
#include "math/bcl_math_comparisons.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "math/bcl_math_running_average_sd.h"
#include "math/bcl_math_transformation_matrix_3d.h"
#include "quality/bcl_quality_rmsd.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> RMSDPreprocessor::s_Instance
    (
      GetObjectInstances().AddInstance( new RMSDPreprocessor())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a boolean to whether superimpose coordinates (set to true by default)
    RMSDPreprocessor::RMSDPreprocessor( const util::SiPtrVector< const linal::Vector3D> &COORDINATES, const bool &RECENTER) :
      m_SquareCenteredNorm( 0.0),
      m_CoordinatesTransposed( size_t( 3), COORDINATES.GetSize()),
      m_Center( size_t( 3), float( 0.0)),
      m_RecenterForSuperimposition( RECENTER)
    {
      size_t coord_id( 0);
      for
      (
        auto itr_coords( COORDINATES.Begin()), itr_coords_end( COORDINATES.End());
        itr_coords != itr_coords_end;
        ++itr_coords, ++coord_id
      )
      {
        linal::Vector< float> coords( ( *itr_coords)->Begin(), ( *itr_coords)->End());
        m_Center += coords;
        m_CoordinatesTransposed( 0, coord_id) = coords( 0);
        m_CoordinatesTransposed( 1, coord_id) = coords( 1);
        m_CoordinatesTransposed( 2, coord_id) = coords( 2);
      }
      m_Center /= float( COORDINATES.GetSize());
      if( m_RecenterForSuperimposition)
      {
        auto x_trans( m_CoordinatesTransposed.GetRow( 0));
        auto y_trans( m_CoordinatesTransposed.GetRow( 1));
        auto z_trans( m_CoordinatesTransposed.GetRow( 2));
        x_trans -= m_Center( 0);
        y_trans -= m_Center( 1);
        z_trans -= m_Center( 2);
        m_SquareCenteredNorm = m_CoordinatesTransposed.AsVector().SquareNorm();
      }
    }

    //! @brief virtual copy constructor
    //! @return pointer to new RMSD
    RMSDPreprocessor *RMSDPreprocessor::Clone() const
    {
      return new RMSDPreprocessor( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &RMSDPreprocessor::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RMSDPreprocessor::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Center, ISTREAM);
      io::Serialize::Read( m_CoordinatesTransposed, ISTREAM);
      m_SquareCenteredNorm = m_CoordinatesTransposed.AsVector().SquareNorm();

      // end
      return ISTREAM;
    }

    //! @brief read from std::ostream
    //! @param OSTREAM input stream
    //! @param INDENT indentation
    //! @return ostream which was read from
    std::ostream &RMSDPreprocessor::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Center, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_CoordinatesTransposed, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief determine the transformation matrix to optimally (lowest RMSD_ superimpose two set of coordinates)
    //! @param COORDINATES_A set of coordinates A
    //! @param COORDINATES_B set of coordinates B
    //! @return Transformation matrix that superimposes B onto A
    math::TransformationMatrix3D RMSDPreprocessor::SuperimposeCoordinates
    (
      const RMSDPreprocessor &COORDINATES_B
    ) const
    {
      BCL_Assert
      (
        m_RecenterForSuperimposition && COORDINATES_B.m_RecenterForSuperimposition,
        "Tried to superimpose coordinates without recentering; indicates constructor was called improperly"
      );
      // check size of vectors
      if( m_CoordinatesTransposed.GetNumberCols() != COORDINATES_B.m_CoordinatesTransposed.GetNumberCols())
      {
        BCL_MessageCrt
        (
          "number of points differs: " + util::Format()( m_CoordinatesTransposed.GetNumberCols()) +
          " != " + util::Format()( COORDINATES_B.m_CoordinatesTransposed.GetNumberCols())
        );

        return math::TransformationMatrix3D( util::UndefinedObject());
      }

      // check for minimal size
      if( m_CoordinatesTransposed.GetNumberCols() < 3)
      {
        return math::TransformationMatrix3D( util::UndefinedObject());
      }

      // Calculate the covariance matrix
      linal::Matrix3x3< float> moment( BuildCovarianceMatrix( COORDINATES_B));
      if( !moment.IsDefined())
      {
        BCL_MessageVrb( "covariance matrix is undefined");
        return math::TransformationMatrix3D( util::UndefinedObject());
      }

      linal::Matrix3x3< double> moment_dbl;
      std::copy( moment.Begin(), moment.End(), moment_dbl.Begin());
      linal::Vector3D center_dbl( m_Center.Begin(), m_Center.End());
      linal::Vector3D b_center_dbl( COORDINATES_B.m_Center.Begin(), COORDINATES_B.m_Center.End());
      return RMSD::CovarianceToTransformationMatrix( moment_dbl, center_dbl, b_center_dbl);
    }

    //! @brief calculate the rmsd of two sets of coordinates if they are optimally
    //! they are never superimposed - the rmsd is derived from the eigenvalues of the SVD of the covariance matrix
    //! @param COORDINATES_A set of coordinates A
    //! @param COORDINATES_B set of coordinates B
    //! @return the rmsd of the coordinates
    double RMSDPreprocessor::SuperimposedRMSD
    (
      const RMSDPreprocessor &COORDINATES_B
    ) const
    {
      BCL_Assert
      (
        m_RecenterForSuperimposition && COORDINATES_B.m_RecenterForSuperimposition,
        "Tried to superimpose coordinates without recentering; indicates constructor was called improperly"
      );
      if( !GetSize() || !COORDINATES_B.GetSize())
      {
        return util::GetUndefined< double>();
      }

      // Calculate the covariance matrix
      linal::Matrix3x3< float> moment( BuildCovarianceMatrix( COORDINATES_B));

      if( !moment.IsDefined())
      {
        BCL_MessageCrt( "covariance matrix is undefined");

        return util::GetUndefinedDouble();
      }

      // determine sign of last element
      static const float s_chi_threshold( 1e-10);
      const int chi( moment.Determinant() < s_chi_threshold ? -1 : 1);

      moment = linal::MatrixTimesItselfTransposed( moment);
      // sort diagonal
      linal::Vector< float> eigenvalues( moment.EigenValues());
      // handle numerical issues that could cause one of the eigenvalues to become slightly less than 0; also take
      // square root
      for
      (
        linal::Vector< float>::iterator itr( eigenvalues.Begin()), itr_end( eigenvalues.End());
        itr != itr_end; ++itr
      )
      {
        *itr = math::Sqrt( std::max( *itr, float( 0.0)));
      }
      std::sort( eigenvalues.Begin(), eigenvalues.End());
      eigenvalues( 0) *= chi;

      // calculate the square deviation
      float square_deviation( m_SquareCenteredNorm + COORDINATES_B.m_SquareCenteredNorm - 2.0 * eigenvalues.Sum());

      // root mean and return
      return math::Sqrt( std::max( square_deviation, float( 0)) / float( m_CoordinatesTransposed.GetNumberCols()));
    }

    //! @brief calculate the rmsd of two sets of coordinates if they are optimally
    //! they are never superimposed - the rmsd is derived from the eigenvalues of the SVD of the covariance matrix
    //! @param COORDINATES_A set of coordinates A
    //! @param COORDINATES_B set of coordinates B
    //! @return the rmsd of the coordinates
    double RMSDPreprocessor::RMSD( const RMSDPreprocessor &COORDINATES_B) const
    {
      BCL_Assert
      (
        !m_RecenterForSuperimposition && !COORDINATES_B.m_RecenterForSuperimposition,
        "Tried to superimpose coordinates but did recentering; indicates constructor was called improperly"
      );
      // check size of vectors
      if( GetSize() != COORDINATES_B.GetSize())
      {
        return util::GetUndefined< double>();
      }

      // compute RMSD
      double sum_square_rmsd( 0);
      for
      (
        auto itr_a( m_CoordinatesTransposed.Begin()), itr_a_end( m_CoordinatesTransposed.End()),
             itr_b( COORDINATES_B.m_CoordinatesTransposed.Begin());
        itr_a != itr_a_end;
        ++itr_a, ++itr_b
      )
      {
        sum_square_rmsd += math::Sqr( *itr_a - *itr_b);
      }

      // end
      return math::Sqrt( sum_square_rmsd / float( 3 * GetSize()));
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief compute covariance matrix of two sets of coordinates COORDINATES_A on COORDINATES_B
    //! both coordinate sets are translated to the center of mass
    linal::Matrix3x3< float> RMSDPreprocessor::BuildCovarianceMatrix( const RMSDPreprocessor &COORDINATES_B) const
    {
      // check the centers are defined
      if( !m_Center.IsDefined() || !COORDINATES_B.m_Center.IsDefined())
      {
        BCL_MessageVrb
        (
          "given centers are undefined:\n" + util::Format()( m_Center) + "\n" + util::Format()( COORDINATES_B.m_Center)
        );

        // return empty matrix
        return linal::Matrix3x3< float>( util::GetUndefined< float>());
      }
      BCL_Assert
      (
        m_CoordinatesTransposed.GetNumberCols() == COORDINATES_B.m_CoordinatesTransposed.GetNumberCols(),
        "Non-equal sized sets of coordinates cannot be used to build a covariance matrix!"
      );

      // make cross moments matrix
      const linal::Matrix< float> covariance_matrix
      (
        linal::MatrixTimesMatrixTranspose( m_CoordinatesTransposed, COORDINATES_B.m_CoordinatesTransposed)
      );

      // end
      return covariance_matrix;
    }

  } // namespace quality
} // namespace bcl
