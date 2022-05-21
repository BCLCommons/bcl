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
#include "quality/bcl_quality_rmsd.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_operations.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "linal/bcl_linal_vector_operations.h"
#include "linal/bcl_linal_vector_reference.h"
#include "math/bcl_math_comparisons.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "math/bcl_math_running_average_sd.h"
#include "math/bcl_math_transformation_matrix_3d.h"
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
    const util::SiPtr< const util::ObjectInterface> RMSD::s_Instance
    (
      GetObjectInstances().AddInstance( new RMSD())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a boolean to whether superimpose coordinates (set to true by default)
    //! @param SUPERIMPOSE_COORDINATES boolean to whether superimpose coordinates before calculating RMSD
    //! @param IGNORE_Z_COORDINATES boolean whether to superimpose using the Z-coordinates
    RMSD::RMSD( const bool SUPERIMPOSE_COORDINATES, const bool IGNORE_Z_COORDINATES) :
      m_SuperimposeCoordinates( SUPERIMPOSE_COORDINATES),
      m_IgnoreZCoordinates( IGNORE_Z_COORDINATES)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to new RMSD
    RMSD *RMSD::Clone() const
    {
      return new RMSD( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &RMSD::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the comparison function for better quality
    //! @return binary function to compare two quality measure values
    const util::BinaryFunctionInterface< double, double, bool> &RMSD::GetComparisonFunction() const
    {
      return **math::Comparisons< double>::GetEnums().e_Less;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief creates and returns a new coordinate set with all Z coordinates set to 0 for the given coordinate set
    //! @param COORDINATES vector of coordinates of interest
    //! @return new vector of coordinates with all Z coordinates set to 0
    storage::Vector< linal::Vector3D> RMSD::RemoveZCoordinates
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES
    ) const
    {
      // initialize vector
      storage::Vector< linal::Vector3D> xy_coordinates;
      xy_coordinates.AllocateMemory( COORDINATES.GetSize());

      // iterate over given coordinates
      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator
          itr( COORDINATES.Begin()), itr_end( COORDINATES.End());
        itr != itr_end; ++itr
      )
      {
        // create new coordinate set without the Z coordinate and push it back
        xy_coordinates.PushBack( linal::Vector3D( ( *itr)->X(), ( *itr)->Y(), 0.0));
      }

      // end
      return xy_coordinates;
    }

    //! @brief calculates root mean square deviation between given coordinates
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return root mean square deviation between given coordinates
    double RMSD::CalculateMeasure
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // if superimpose coordinates is set
      if( m_SuperimposeCoordinates)
      {
        // if Z coordinates are to be used for the superimposition
        if( !m_IgnoreZCoordinates)
        {
          return SuperimposedRMSD( COORDINATES, REFERENCE_COORDINATES);
        }
        // Z coordinates are not to be used for the superimposition
        else
        {
          // create vectors for the new coordinates
          storage::Vector< linal::Vector3D> xy_coordinates( RemoveZCoordinates( COORDINATES));
          storage::Vector< linal::Vector3D> reference_xy_coordinates( RemoveZCoordinates( REFERENCE_COORDINATES));

          // calculate the superimposition
          const math::TransformationMatrix3D transform
          (
            SuperimposeCoordinates
            (
              util::ConvertToConstSiPtrVector( reference_xy_coordinates),
              util::ConvertToConstSiPtrVector( xy_coordinates)
            )
          );

          // create a vector to hold the transformed coordinates
          storage::Vector< linal::Vector3D> new_coordinates;

          // apply the transformation to the b coordinates
          for
          (
            util::SiPtrVector< const linal::Vector3D>::const_iterator coord_itr( COORDINATES.Begin()),
              coord_itr_end( COORDINATES.End());
            coord_itr != coord_itr_end; ++coord_itr
          )
          {
            linal::Vector3D this_coordinate( **coord_itr);
            this_coordinate.Transform( transform);
            new_coordinates.PushBack( this_coordinate);
          }

          // calculate the RMSD and return
          return RMSD::RealSpaceRMSD( util::ConvertToConstSiPtrVector( new_coordinates), REFERENCE_COORDINATES);
        }
      }
      // don't superimpose coordinates
      else
      {
        return RMSD::RealSpaceRMSD( COORDINATES, REFERENCE_COORDINATES);
      }
    }

    //! @brief calculates the transformation matrix that superimposes COORDINATES onto REFERENCE_COORDINATES
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return the transformation matrix that superimposes COORDINATES onto REFERENCE_COORDINATES
    math::TransformationMatrix3D RMSD::CalculateSuperimposition
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // if Z coordinates are to be ignore for the superimposition
      if( m_IgnoreZCoordinates)
      {
        // make copies of the coordinates without Z coordinates
         storage::Vector< linal::Vector3D> xy_coordinates( RemoveZCoordinates( COORDINATES));
         storage::Vector< linal::Vector3D> reference_xy_coordinates( RemoveZCoordinates( REFERENCE_COORDINATES));

         // calculate the transformation matrix and return it
         return
           SuperimposeCoordinates
           (
             util::ConvertToConstSiPtrVector( reference_xy_coordinates),
             util::ConvertToConstSiPtrVector( xy_coordinates)
           );
      }
      // otherwise do regular superimposition
      else
      {
        return SuperimposeCoordinates( REFERENCE_COORDINATES, COORDINATES);
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RMSD::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SuperimposeCoordinates, ISTREAM);
      io::Serialize::Read( m_IgnoreZCoordinates, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief read from std::ostream
    //! @param OSTREAM input stream
    //! @param INDENT indentation
    //! @return ostream which was read from
    std::ostream &RMSD::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SuperimposeCoordinates, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_IgnoreZCoordinates, OSTREAM, INDENT);

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
    math::TransformationMatrix3D RMSD::SuperimposeCoordinates
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES_A,
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES_B
    )
    {
      // check size of vectors
      if( COORDINATES_A.GetSize() != COORDINATES_B.GetSize())
      {
        BCL_MessageCrt
        (
          "number of points differs: " + util::Format()( COORDINATES_A.GetSize()) +
          " != " + util::Format()( COORDINATES_B.GetSize())
        );

        return math::TransformationMatrix3D( util::UndefinedObject());
      }

      const size_t number( COORDINATES_A.GetSize());

      // check for minimal size
      if( number < 3)
      {
//        BCL_Message
//        (
//          util::Message::e_Critical,
//          "number points should be at least 3, otherwise ambiguous solution for super imposition"
//        );

        return math::TransformationMatrix3D( util::UndefinedObject());
      }

      // calculate the center of mass
      const linal::Vector3D shift_a( coord::CenterOfMass( COORDINATES_A, false));
      const linal::Vector3D shift_b( coord::CenterOfMass( COORDINATES_B, false));

      // Calculate the covariance matrix
      linal::Matrix3x3< double> moment( BuildCovarianceMatrix( COORDINATES_A, COORDINATES_B, shift_a, shift_b));
      if( !moment.IsDefined())
      {
        BCL_MessageVrb( "covariance matrix is undefined");
        return math::TransformationMatrix3D( util::UndefinedObject());
      }

      return CovarianceToTransformationMatrix( moment, shift_a, shift_b);
    }

    //! @brief calculate the real space rmsd of two sets of coordinates
    //! uses the coordinates as they are passed
    //! @param COORDINATES_A set of coordinates A
    //! @param COORDINATES_B set of coordinates B
    //! @return the rmsd between the passed coordinates
    double RMSD::RealSpaceRMSD
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES_A,
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES_B
    )
    {
      // check size of vectors
      if( COORDINATES_A.GetSize() != COORDINATES_B.GetSize())
      {
        return util::GetUndefined< double>();
      }

      // compute RMSD
      double sum_square_rmsd( 0);
      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator
          itr_a( COORDINATES_A.Begin()), itr_a_end( COORDINATES_A.End()),
          itr_b( COORDINATES_B.Begin()), itr_b_end( COORDINATES_B.End());
        itr_a != itr_a_end && itr_b != itr_b_end;
        ++itr_a, ++itr_b
      )
      {
        sum_square_rmsd += ( **itr_a - **itr_b).SquareNorm();
      }

      // end
      return math::Sqrt( sum_square_rmsd / COORDINATES_A.GetSize());
    }

    //! @brief calculate the real space rmsd for a set of coordinates - compares each coordinate with every other
    //! uses the coordinates as they are passed
    //! @param COORDINATES set of coordinates that will be compared with themselves
    //! @return first the rmsd between the passed coordinates and the standard deviation of the rmsd
    //!         the RunningAverageSD< double> has the mean distance and the standard deviation in the distances
    storage::Pair< storage::VectorND< 2, double>, math::RunningAverageSD< double> > RMSD::RealSpaceRMSDPairwise
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES
    )
    {
      // true if only one coordinate or less
      if( COORDINATES.GetSize() < 2)
      {
        return storage::Pair< storage::VectorND< 2, double>, math::RunningAverageSD< double> >
        (
          storage::VectorND< 2, double>( util::GetUndefinedDouble(), util::GetUndefinedDouble()),
          math::RunningAverageSD< double>()
        );
      }

      math::RunningAverageSD< double> mean_sd;

      // compute RMSD
      double sum_square_rmsd( 0);

      // iterate over coordinates
      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator
          itr_a( COORDINATES.Begin()), itr_a_end( COORDINATES.End());
        itr_a != itr_a_end;
        ++itr_a
      )
      {
        // iterate over coordinates
        for
        (
          util::SiPtrVector< const linal::Vector3D>::const_iterator itr_b( itr_a + 1);
          itr_b != itr_a_end;
          ++itr_b
        )
        {
          const double distance( ( **itr_a - **itr_b).Norm());
          BCL_MessageDbg( " distance " + util::Format()( distance));
          sum_square_rmsd += math::Sqr( distance);
          mean_sd += distance;
        }
      }

      // number pairwise comparisons
      const size_t num_comparisons( COORDINATES.GetSize() * ( COORDINATES.GetSize() - 1) / 2);
      BCL_MessageDbg( " sum_square_rmsd " + util::Format()( sum_square_rmsd));
      BCL_MessageDbg( " num_comparisons " + util::Format()( num_comparisons));
      // calculate standard deviation of rmsd
      const double rmsd_std_dev
      (
        math::Sqrt( std::abs( sum_square_rmsd / num_comparisons - math::Sqr( mean_sd.GetAverage())))
      );

      // calculate the rmsd
      const double rmsd( math::Sqrt( sum_square_rmsd / num_comparisons));

      // end
      return storage::Pair< storage::VectorND< 2, double>, math::RunningAverageSD< double> >
      (
        storage::VectorND< 2, double>( rmsd, rmsd_std_dev), mean_sd
      );
    }

    //! @brief calculate the rmsd of two sets of coordinates if they are optimally
    //! they are never superimposed - the rmsd is derived from the eigenvalues of the SVD of the covariance matrix
    //! @param COORDINATES_A set of coordinates A
    //! @param COORDINATES_B set of coordinates B
    //! @return the rmsd of the coordinates
    double RMSD::SuperimposedRMSD
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES_A,
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES_B
    )
    {
      if( COORDINATES_A.IsEmpty() || COORDINATES_B.IsEmpty())
      {
        return util::GetUndefined< double>();
      }

      double square_norm_centered_a( 0.0);
      double square_norm_centered_b( 0.0);

      // Calculate the covariance matrix
      linal::Matrix3x3< double> moment
      (
        BuildCovarianceMatrix
        (
          COORDINATES_A,
          COORDINATES_B,
          coord::CenterOfMass( COORDINATES_A, false),
          coord::CenterOfMass( COORDINATES_B, false),
          &square_norm_centered_a,
          &square_norm_centered_b
        )
      );

      if( !moment.IsDefined())
      {
        BCL_MessageCrt( "covariance matrix is undefined");

        return util::GetUndefinedDouble();
      }

      // determine sign of last element
      static const double s_chi_threshold( 1e-10);
      const int chi( moment.Determinant() < s_chi_threshold ? -1 : 1);

      moment = linal::MatrixTimesItselfTransposed( moment);
      // sort diagonal
      linal::Vector< double> eigenvalues( moment.EigenValues());
      // handle numerical issues that could cause one of the eigenvalues to become slightly less than 0; also take
      // square root
      for
      (
        linal::Vector< double>::iterator itr( eigenvalues.Begin()), itr_end( eigenvalues.End());
        itr != itr_end; ++itr
      )
      {
        *itr = math::Sqrt( std::max( *itr, double( 0.0)));
      }
      std::sort( eigenvalues.Begin(), eigenvalues.End());
      eigenvalues( 0) *= chi;

      // calculate the square deviation
      double square_deviation( 0.0);
      square_deviation += square_norm_centered_a;
      square_deviation += square_norm_centered_b;
      square_deviation -= 2 * eigenvalues.Sum();

      // root mean and return
      return math::Sqrt( std::max( square_deviation, double( 0)) / double( COORDINATES_A.GetSize()));
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief compute covariance matrix of two sets of coordinates COORDINATES_A on COORDINATES_B
    //! both coordinate sets are translated to the center of mass
    //! @param COORDINATES_A set of coordinates
    //! @param COORDINATES_B set of coordinates
    //! @param CENTER_A the center of COORDINATES_A
    //! @param CENTER_B the center of COORDINATES_B
    //! @param SQUARE_NORM_CENTERED_COORDINATES_A optional pointer to which the square norm of the centered coordinates a will be deposited
    //! @param SQUARE_NORM_CENTERED_COORDINATES_B optional pointer to which the square norm of the centered coordinates b will be deposited
    //! @return COORDINATES_A * COORDINATES_B
    linal::Matrix3x3< double> RMSD::BuildCovarianceMatrix
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES_A,
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES_B,
      const linal::Vector3D &CENTER_A,
      const linal::Vector3D &CENTER_B,
      double *SQUARE_NORM_CENTERED_COORDINATES_A,
      double *SQUARE_NORM_CENTERED_COORDINATES_B
    )
    {
      // check the centers are defined
      if( !CENTER_A.IsDefined() || !CENTER_B.IsDefined())
      {
        BCL_MessageVrb
        (
          "given centers are undefined:\n" + util::Format()( CENTER_A) + "\n" + util::Format()( CENTER_B)
        );

        // return empty matrix
        return linal::Matrix3x3< double>( util::GetUndefined< double>());
      }
      BCL_Assert
      (
        COORDINATES_A.GetSize() == COORDINATES_B.GetSize(),
        "Non-equal sized sets of coordinates cannot be used to build a covariance matrix!"
      );

      const size_t number( COORDINATES_A.GetSize());

      // initialize 2 matrices with coordinates and 2 to keep the start coordinates
      linal::Matrix< double> coord_a( number, 3);
      linal::Matrix< double> coord_b( number, 3);

      // copy and shift coordinates

      for( size_t row( 0); row < number; ++row)
      {
        linal::VectorReference< double> row_a( coord_a.GetRow( row)), row_b( coord_b.GetRow( row));
        row_a.CopyValues( *COORDINATES_A( row));
        row_a -= CENTER_A;
        row_b.CopyValues( *COORDINATES_B( row));
        row_b -= CENTER_B;
      }

      // make cross moments matrix
      const linal::Matrix< double> covariance_matrix( linal::MatrixTransposeTimesMatrix( coord_a, coord_b));

      // set the argument centered coordinates if desired
      if( SQUARE_NORM_CENTERED_COORDINATES_A != NULL)
      {
        *SQUARE_NORM_CENTERED_COORDINATES_A = coord_a.AsVector().SquareNorm();
      }
      if( SQUARE_NORM_CENTERED_COORDINATES_B != NULL)
      {
        *SQUARE_NORM_CENTERED_COORDINATES_B = coord_b.AsVector().SquareNorm();
      }

      // end
      return covariance_matrix;
    }

    //! @brief Transformation matrix from Covariance matrix
    //! @param MOMENT covariance matrix
    //! @param CENTER_COORDINATES of coordinates
    //! @param CENTER_REFERENCE_COORDINATES center of reference coordinates
    //! @return transformation matrix
    math::TransformationMatrix3D RMSD::CovarianceToTransformationMatrix
    (
      const linal::Matrix3x3< double> &MOMENT,
      const linal::Vector3D &CENTER_COORDINATES,
      const linal::Vector3D &CENTER_REFERENCE_COORDINATES
    )
    {
      // diagonalization
      linal::Matrix3x3< double> rotate( MOMENT);
      linal::Matrix3x3< double> moment_t( MOMENT);
      moment_t.Transpose();
      rotate *= moment_t;

      // solve Eigensystem
      linal::Vector< double> eigenvalues( 3, 0.0);
      linal::Matrix3x3< double> eigenvectors;
      if( !rotate.EigenVectorsSymmetric( eigenvectors, eigenvalues))
      {
        BCL_MessageCrt( "Non-symmetric eigenmatrix! Check: " + util::Format()( eigenvectors) + "; should be A^T*A of " + util::Format()( MOMENT));
        return math::TransformationMatrix3D( util::UndefinedObject());
      }

      if( !util::IsDefined( eigenvalues( 0)))
      {
        BCL_MessageCrt( "Undefined principle eigenvalue (degenerate solutions)!");
        return math::TransformationMatrix3D( util::UndefinedObject());
      } // error

      // sort eigenvectors by eigenvalues
      eigenvectors.Transpose();
      eigenvectors.SortRowsAndVector( eigenvalues);

      // check second eigenvalue
      if( eigenvalues( 1) <= 0.0 || eigenvalues( 0) <= 0.0)
      {
        BCL_MessageCrt( "Undefined secondary eigenvalue (degenerate solutions)!");
        return math::TransformationMatrix3D();
      } // error

      // build rotation matrix
      eigenvectors.ReplaceRow( 2, linal::CrossProduct( linal::Vector3D( eigenvectors[ 0]), linal::Vector3D( eigenvectors[ 1])));
      rotate = eigenvectors * MOMENT;

      std::transform( eigenvalues.Begin(), eigenvalues.End(), eigenvalues.Begin(), &math::Sqrt< double>);
      std::transform( eigenvalues.Begin(), eigenvalues.End(), eigenvalues.Begin(), std::bind1st( std::divides< double>(), 1.0));
      for( double *ptr( rotate.Begin()), *ptr_e( eigenvalues.Begin()), *ptr_end( rotate.End()); ptr < ptr_end; ptr += rotate.GetNumberCols(), ++ptr_e)
      {
        std::transform( ptr, ptr + rotate.GetNumberCols(), ptr, std::bind2nd( std::multiplies< double>(), *ptr_e));
      }

      rotate.ReplaceRow( 2, linal::CrossProduct( linal::Vector3D( rotate[ 0]), linal::Vector3D( rotate[ 1])));
      rotate.Transpose();
      rotate *= eigenvectors;

      // shift and rotate molecule
      math::TransformationMatrix3D transform;
      transform( -CENTER_REFERENCE_COORDINATES);
      transform( math::RotationMatrix3D( rotate));
      transform( CENTER_COORDINATES);

      // return transformation matrix calculated
      return transform;
    }

  } // namespace quality
} // namespace bcl
