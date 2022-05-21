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

#ifndef BCL_OPENCL_QUALITY_RMSD_H_
#define BCL_OPENCL_QUALITY_RMSD_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_opencl_command_queue.h"
#include "quality/bcl_quality_superimpose_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class QualityRMSD
    //! @brief QualityRMSD is a class which calculates the root mean square deviation between two sets of coordinates
    //! @details QualityRMSD class calculates the root mean square deviation between two sets of coordinates. If superimposition
    //! flag is given, then the returned value is the RMSD between superimposed coordinates.
    //! An additional flags provides, superimposition without Z coordinates which can be used in comparing structures
    //! for membrane proteins, where the superimposition should only be in XY coordinates.
    //!
    //! @see @link example_opencl_quality_rmsd.cpp @endlink
    //! @author woetzen
    //! @date 01/15/2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API QualityRMSD :
      public quality::SuperimposeInterface
    {

    private:

    //////////
    // data //
    //////////

      //! opencl queue
      CommandQueue m_Queue;

      //! opencl program
      cl::Program m_Program;

      //! boolean to whether superimpose coordinates before calculating RMSD
      bool         m_SuperimposeCoordinates;

    public:

      static const size_t s_BlockSize;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a boolean to whether superimpose coordinates (set to true by default)
      //! @param SUPERIMPOSE_COORDINATES boolean to whether superimpose coordinates before calculating RMSD
      QualityRMSD
      (
        const bool SUPERIMPOSE_COORDINATES = true
      );

      //! @brief constructor from a boolean to whether superimpose coordinates (set to true by default)
      //! @param QUEUE command queue
      //! @param SUPERIMPOSE_COORDINATES boolean to whether superimpose coordinates before calculating RMSD
      QualityRMSD
      (
        const CommandQueue &QUEUE,
        const bool SUPERIMPOSE_COORDINATES = true
      );

      //! @brief virtual copy constructor
      //! @return pointer to new QualityRMSD
      QualityRMSD *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get name of the current class
      //! @return name of the class
      const std::string &GetClassIdentifier() const;

      //! @brief access to the program
      const cl::Program &GetProgram() const
      {
        return m_Program;
      }

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

    ////////////////
    // operations //
    ////////////////

      //! @brief is this class compatible with given command queue
      //! @param COMMAND_QUEUE the command queue this object would operate on
      //! @return bool true, if this command queue can be used, e.g. the devices support the necessary extensions
      bool IsCompatible( const CommandQueue &COMMAND_QUEUE) const;

      //! @brief initialize this class
      //! @brief COMMAND_QUEUE queue to use
      //! @return bool true is initialization was successful, false otherwise (when program could not be compiled ...)
      bool Initialize( const CommandQueue &COMMAND_QUEUE);

      //! @brief calculates root mean square deviation between given coordinates
      //! @param COORDINATES coordinates of interest
      //! @param REFERENCE_COORDINATES reference coordinates
      //! @return root mean square deviation between given coordinates
      double CalculateMeasure
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
      ) const;

      //! @brief calculates root mean square deviation between given coordinates
      //! @param COORDINATES matrix of coordinates of interest
      //! @param REFERENCE_COORDINATES matrix of reference coordinates
      //! @return root mean square deviation between given coordinates
      double CalculateMeasure
      (
        const Matrix< double> &COORDINATES,
        const Matrix< double> &REFERENCE_COORDINATES
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

      //! @brief calculates the transformation matrix that superimposes COORDINATES onto REFERENCE_COORDINATES
      //! @param COORDINATES matrix of coordinates of interest
      //! @param REFERENCE_COORDINATES matrix of reference coordinates
      //! @return the transformation matrix that superimposes COORDINATES onto REFERENCE_COORDINATES
      math::TransformationMatrix3D CalculateSuperimposition
      (
        const Matrix< double> &COORDINATES,
        const Matrix< double> &REFERENCE_COORDINATES
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

      //! @brief determine the transformation matrix to optimally (lowest RMSD) superimpose two sets of coordinates
      //! @param COORDINATES matrix of coordinates of interest
      //! @param REFERENCE_COORDINATES matrix of reference coordinates
      //! @return Transformation matrix that superimposes B onto A
      Matrix< double> SuperimposeCoordinates
      (
        const Matrix< double> &REFERENCE_COORDINATES,
        const Matrix< double> &COORDINATES
      ) const;

      //! @brief calculate the real space rmsd of two sets of coordinates
      //! uses the coordinates as they are passed
      //! @param COORDINATES matrix of coordinates of interest
      //! @param REFERENCE_COORDINATES matrix of reference coordinates
      //! @return the rmsd between the passed coordinates
      double RealSpaceRMSD
      (
        const Matrix< double> &COORDINATES,
        const Matrix< double> &REFERENCE_COORDINATES
      ) const;

      //! @brief calculate the rmsd of two sets of coordinates if they are optimally
      //! they are never superimposed - the rmsd is derived from the eigenvalues of the SVD of the covariance matrix
      //! @param COORDINATES matrix of coordinates of interest
      //! @param REFERENCE_COORDINATES matrix of reference coordinates
      //! @return the rmsd of the coordinates
      double SuperimposedRMSD
      (
        const Matrix< double> &COORDINATES,
        const Matrix< double> &REFERENCE_COORDINATES
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief create a padded matrix for a vector of coordinates
      //! @param COORDINATES the coordinates
      //! @param BLOCK_SIZE block size for kernels
      //! @return Matrix that are padded to desired blocksize and contain the coordinates as rows
      Matrix< double> MatrixFromCoordinates( const util::SiPtrVector< const linal::Vector3D> &COORDINATES, const size_t BLOCK_SIZE) const;

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief compute covariance matrix of two sets of coordinates COORDINATES_A on COORDINATES_B
      //! both coordinate sets are translated to the center of mass
      //! @param COORDINATES matrix of coordinates of interest
      //! @param REFERENCE_COORDINATES matrix of reference coordinates
      //! @param CENTER_A the center of COORDINATES_A
      //! @param CENTER_B the center of COORDINATES_B
      //! @param SQUARE_NORM_CENTERED_COORDINATES_A optional pointer to which the square norm of the centered coordinates a will be depsosited
      //! @param SQUARE_NORM_CENTERED_COORDINATES_B optional pointer to which the square norm of the centered coordinates b will be depsosited
      //! @return COORDINATES_A * COORDINATES_B
      Matrix3x3< double> BuildCovarianceMatrix
      (
        const Matrix< double> &COORDINATES,
        const Matrix< double> &REFERENCE_COORDINATES,
        const Vector< double> &CENTER_A,
        const Vector< double> &CENTER_B,
        Vector< double> &SQUARE_NORM_CENTERED_COORDINATES_A,
        Vector< double> &SQUARE_NORM_CENTERED_COORDINATES_B
      ) const;

      //! @brief calculate the center of a given matrix of coordinates
      //! @param COORDINATES matrix of coordinates in rows
      //! @return the center as vector
      Vector< double> Center
      (
        const Matrix< double> &COORDINATES
      ) const;

      //! @brief Transformation matrix from Covariance matrix
      //! @param MOMENT covariance matrix
      //! @param CENTER_COORDINATES of coordinates
      //! @param CENTER_REFERENCE_COORDINATES center of reference coordinates
      //! @return transformation matrix
      Matrix< double> CovarianceToTransformationMatrix
      (
        const Matrix3x3< double> &MOMENT,
        const Vector< double> &CENTER_COORDINATES,
        const Vector< double> &CENTER_REFERENCE_COORDINATES
      ) const;

    }; // class QualityRMSD

  } // namespace opencl
} // namespace bcl

#endif //BCL_OPENCL_QUALITY_RMSD_H_
