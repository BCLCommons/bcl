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

#ifndef BCL_OPENCL_QUALITY_GDT_H_
#define BCL_OPENCL_QUALITY_GDT_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opencl_quality_lcs.h"
#include "quality/bcl_quality_gdt.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class QualityGDT
    //! @brief opencl implementation of GDT quality measure
    //! @details matrix of coordinates is transferred to the device and almost all operations are performed there
    //!
    //! @see @link example_opencl_quality_gdt.cpp @endlink
    //! @author woetzen
    //! @date Jan 13, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API QualityGDT :
      public quality::SuperimposeInterface
    {

    private:

    //////////
    // data //
    //////////

      //! opencl queue
      CommandQueue m_Queue;

      //! quality lcs object
      QualityLCS m_QualityLCS;

    public:

      static const size_t s_BlockSize = 16;
      static const size_t s_NumberRowsTransformationMatrix;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from a single distance cutoff and seed length
      //! @param QUEUE command queue
      //! @param DISTANCE_CUTOFF distance cutoff to be used
      //! @param SEED_LENGTH length of seed
      QualityGDT
      (
        const double &DISTANCE_CUTOFF,
        const size_t SEED_LENGTH = quality::GDT::GetDefaultSeedLength()
      );

      //! @brief construct from a single distance cutoff and seed length
      //! @param QUEUE command queue
      //! @param DISTANCE_CUTOFF distance cutoff to be used
      //! @param SEED_LENGTH length of seed
      QualityGDT
      (
        const CommandQueue &QUEUE,
        const double &DISTANCE_CUTOFF,
        const size_t SEED_LENGTH = quality::GDT::GetDefaultSeedLength()
      );

      //! @brief Clone function
      //! @return pointer to new QualityGDT
      QualityGDT *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the optimal value for that quality measurement
      //! @return the best value by which two sets of coordinates can agree
      double OptimalValue() const
      {
        return 100.0;
      }

      //! @brief return the comparison function for better quality
      //! @return binary function to compare two quality measure values
      const util::BinaryFunctionInterface< double, double, bool> &GetComparisonFunction() const;

      //! @brief get seed length
      //! @return seed length
      size_t GetSeedLength() const;

      //! @brief get the distance cutoff
      //! @return distance cutoff
      double GetDistanceCutoff() const;

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

      //! @brief calculates GDT between COORDINATES and REFERENCE_COORDINATES
      //! @param COORDINATES coordinates of interest
      //! @param REFERENCE_COORDINATES reference coordinates
      //! @return GDT between COORDINATES and REFERENCE_COORDINATES
      double CalculateMeasure
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
      ) const;

      //! @brief calculates the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
      //! @param COORDINATES coordinates of interest
      //! @param REFERENCE_COORDINATES reference coordinates
      //! @return the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
      math::TransformationMatrix3D CalculateSuperimposition
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
      ) const;

      //! @brief calculates GDT and the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
      //! @param COORDINATES coordinates of interest
      //! @param REFERENCE_COORDINATES reference coordinates
      //! @return pair of GDT and the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
      storage::Pair< double, math::TransformationMatrix3D> CalculateGDTAndSuperimposition
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
      ) const;

      //! @brief calculates GDT and the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
      //! @param COORDINATES matrix of coordinates of interest
      //! @param REFERENCE_COORDINATES matrix of reference coordinates
      //! @return pair of GDT and the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
      storage::Pair< double, Matrix< double> > CalculateGDTAndSuperimposition
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

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief returns all possible seed ranges as selections in cols (not necessarily valid)
      //! @param NUMBER_OF_COORDINATES number of coordinates and rows in the selection matrix
      //! @param NUMBER_OF_FRAGMENTS_RND rounded number of fragments - number of cols in selections matrix
      //! @return Matrix of selections for coordinate vectors
      Matrix< int> GetSeedSelections
      (
        const size_t NUMBER_OF_COORDINATES,
        const size_t NUMBER_OF_FRAGMENTS_RND
      ) const;

      //! @brief calculate the transformations for superimposition
      //! @param COORDINATES matrix of coordinates of interest
      //! @param REFERENCE_COORDINATES matrix of reference coordinates
      //! @param SELECTIONS matrix of elements indicating which rows to use for each fragment
      //! @param TRANSFORMATIONS transformation matrices for all selections
      //! @param NUMBER_OF_FRAGMENTS number of fragments
      void CalculateTransformations
      (
        const Matrix< double> &COORDINATES,
        const Matrix< double> &REFERENCE_COORDINATES,
        const Matrix< int>    &SELECTIONS,
              Matrix< double> &TRANSFORMATIONS,
        const size_t          NUMBER_OF_FRAGMENTS
      ) const;

      //! @brief update the selection with the given transformation
      //! @param COORDINATES matrix of coordinates of interest
      //! @param REFERENCE_COORDINATES matrix of reference coordinates
      //! @param SELECTIONS matrix of elements indicating which rows to use for each fragment
      //! @param TRANSFORMATIONS transformation matrices for all selections
      //! @param NUMBER_OF_FRAGMENTS number of fragments
      void UpdateSelections
      (
        const Matrix< double> &COORDINATES,
        const Matrix< double> &REFERENCE_COORDINATES,
              Matrix< int>    &SELECTIONS,
        const Matrix< double> &TRANSFORMATIONS,
        const size_t NUMBER_OF_FRAGMENTS
      ) const;

    }; // class QualityGDT

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_QUALITY_GDT_H_ 
