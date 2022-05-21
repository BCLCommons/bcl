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

#ifndef BCL_QUALITY_RMSD_H_
#define BCL_QUALITY_RMSD_H_

// include the namespace header
#include "bcl_quality.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_quality_superimpose_interface.h"
#include "linal/bcl_linal_matrix3x3.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RMSD
    //! @brief RMSD is a class which calculates the root mean square deviation between two sets of coordinates
    //! @details RMSD class calculates the root mean square deviation between two sets of coordinates. If superimposition
    //! flag is given, then the returned value is the RMSD between superimposed coordinates.
    //! An additional flags provides, superimposition without Z coordinates which can be used in comparing structures
    //! for membrane proteins, where the superimposition should only be in XY coordinates.
    //!
    //! @see @link example_quality_rmsd.cpp @endlink
    //! @author alexanns, staritrd, woetzen, karakam
    //! @date 04/07/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RMSD :
      public SuperimposeInterface
    {

    private:

    //////////
    // data //
    //////////

      //! boolean to whether superimpose coordinates before calculating RMSD
      bool m_SuperimposeCoordinates;

      //! boolean whether to superimpose not using the Z-coordinates, only checked when m_SuperImposeCoordinates is true
      bool m_IgnoreZCoordinates;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a boolean to whether superimpose coordinates (set to true by default)
      //! @param SUPERIMPOSE_COORDINATES boolean to whether superimpose coordinates before calculating RMSD
      //! @param IGNORE_Z_COORDINATES boolean whether to superimpose using the Z-coordinates
      RMSD
      (
        const bool SUPERIMPOSE_COORDINATES = true,
        const bool IGNORE_Z_COORDINATES = false
      );

      //! @brief virtual copy constructor
      //! @return pointer to new RMSD
      RMSD *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get name of the current class
      //! @return name of the class
      const std::string &GetClassIdentifier() const;

      //! @brief return the optimal value for that quality measurement
      //! @return the best value by which two sets of coordinates can agree
      double OptimalValue() const
      {
        return 0.0;
      }

      //! @brief return the comparison function for better quality
      //! @return binary function to compare two quality measure values
      const util::BinaryFunctionInterface< double, double, bool> &GetComparisonFunction() const;

      //! @brief whether superimpose coordinates boolean is set
      //! @return whether superimpose coordinates boolean is set
      bool GetSuperimposeCoordinates() const
      {
        return m_SuperimposeCoordinates;
      }

      //! @brief whether to ignore Z coordinates when superimposing
      //! @return whether to ignore Z coordinates when superimposing
      bool GetIgnoreZCoordinates() const
      {
        return m_IgnoreZCoordinates;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief creates and returns a new coordinate set with all Z coordinates set to 0 for the given coordinate set
      //! @param COORDINATES vector of coordinates of interest
      //! @return new vector of coordinates with all Z coordinates set to 0
      storage::Vector< linal::Vector3D> RemoveZCoordinates
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES
      ) const;

      //! @brief calculates root mean square deviation between given coordinates
      //! @param COORDINATES coordinates of interest
      //! @param REFERENCE_COORDINATES reference coordinates
      //! @return root mean square deviation between given coordinates
      double CalculateMeasure
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
      ) const;

      //! @brief calculates the transformation matrix that superimposes COORDINATES onto REFERENCE_COORDINATES
      //! @param COORDINATES coordinates of interest
      //! @param REFERENCE_COORDINATES reference coordinates
      //! @return the transformation matrix that superimposes COORDINATES onto REFERENCE_COORDINATES
      math::TransformationMatrix3D CalculateSuperimposition
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief read from std::ostream
      //! @param OSTREAM input stream
      //! @param INDENT number of indentations
      //! @return ostream which was read from
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

    ////////////////
    // operations //
    ////////////////

      //! @brief determine the transformation matrix to optimally (lowest RMSD_ superimpose two set of coordinates)
      //! @param COORDINATES_A set of coordinates A
      //! @param COORDINATES_B set of coordinates B
      //! @return Transformation matrix that superimposes B onto A
      static math::TransformationMatrix3D SuperimposeCoordinates
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES_A,
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES_B
      );

      //! @brief calculate the real space rmsd of two sets of coordinates
      //! uses the coordinates as they are passed
      //! @param COORDINATES_A set of coordinates A
      //! @param COORDINATES_B set of coordinates B
      //! @return the rmsd between the passed coordinates
      static double RealSpaceRMSD
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES_A,
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES_B
      );

      //! @brief calculate the real space rmsd for a set of coordinates - compares each coordinate with every other
      //! uses the coordinates as they are passed
      //! @param COORDINATES set of coordinates that will be compared with themselves
      //! @return first the rmsd between the passed coordinates and the standard deviation of the rmsd
      //!         the RunningAverageSD< double> has the mean distance and the standard deviation in the distances
      static storage::Pair< storage::VectorND< 2, double>, math::RunningAverageSD< double> > RealSpaceRMSDPairwise
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES
      );

      //! @brief calculate the rmsd of two sets of coordinates if they are optimally
      //! they are never superimposed - the rmsd is derived from the eigenvalues of the SVD of the covariance matrix
      //! @param COORDINATES_A set of coordinates A
      //! @param COORDINATES_B set of coordinates B
      //! @return the rmsd of the coordinates
      static double SuperimposedRMSD
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES_A,
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES_B
      );

    private:

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
      static linal::Matrix3x3< double> BuildCovarianceMatrix
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES_A,
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES_B,
        const linal::Vector3D &CENTER_A,
        const linal::Vector3D &CENTER_B,
        double *SQUARE_NORM_CENTERED_COORDINATES_A = NULL,
        double *SQUARE_NORM_CENTERED_COORDINATES_B = NULL
      );

    public:

      //! @brief Transformation matrix from Covariance matrix
      //! @param MOMENT covariance matrix
      //! @param CENTER_COORDINATES of coordinates
      //! @param CENTER_REFERENCE_COORDINATES center of reference coordinates
      //! @return transformation matrix
      static math::TransformationMatrix3D CovarianceToTransformationMatrix
      (
        const linal::Matrix3x3< double> &MOMENT,
        const linal::Vector3D &CENTER_COORDINATES,
        const linal::Vector3D &CENTER_REFERENCE_COORDINATES
      );

    }; // class RMSD

  } // namespace quality
} // namespace bcl

#endif //BCL_QUALITY_RMSD_H_
