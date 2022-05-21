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

#ifndef BCL_QUALITY_MAXSUB_H_
#define BCL_QUALITY_MAXSUB_H_

// include the namespace header
#include "bcl_quality.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_quality_superimpose_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MaxSub
    //! @brief evaluates the percentage of coordinates that can be superimposed below a threshold
    //! @details Maxsub provides an alternative way to measure the similarity of two coordinate sets of same size.
    //! Unlike RMSD, it tries to find the percentage of coordinates that can be superimposed below a given distance threshold.
    //!
    //! @see Siew et al, MaxSub: an automated measure for the assessment of protein structure prediction quality. Bioinformatics. 2000;16:776–785.
    //!
    //! @see @link example_quality_maxsub.cpp @endlink
    //! @author karakam
    //! @date Oct 8, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MaxSub :
      public SuperimposeInterface
    {

    private:

    //////////
    // data //
    //////////

      //! RMSD cutoff for the longest subset of superimposed coordinates
      double m_RMSDCutoff;

      //! length of the seed subset of coordinates
      size_t m_SeedLength;

      //! number of iterations to extend an individual seed subset
      size_t m_NumberIterations;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief returns default RMSD cutoff
      //! @return default RMSD cutoff
      static double GetDefaultRMSDCutoff();

      //! @brief returns default length of the seed subset of coordinates
      //! @return default seed length of the seed subset of coordinates
      static size_t GetDefaultSeedLength();

      //! @brief returns default number of iterations to extend an individual seed subset
      //! @return default number of iterations to extend an individual seed subset
      static size_t GetDefaultNumberIterations();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a RMSD cutoff, seed length and number of iterations
      //! @param RMSD_CUTOFF RMSD cutoff for the longest subset of superimposed coordinates
      //! @param SEED_LENGTH length of the seed subset of coordinates
      //! @param NUMBER_ITERATIONS number of iterations to extend an individual seed subset
      MaxSub
      (
        const double RMSD_CUTOFF = GetDefaultRMSDCutoff(),
        const size_t SEED_LENGTH = GetDefaultSeedLength(),
        const size_t NUMBER_ITERATIONS = GetDefaultNumberIterations()
      );

      //! @brief Clone function
      //! @return pointer to new MaxSub
      MaxSub *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
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

      //! @brief get RMSD cutoff
      //! @return RMSD cutoff
      double GetRMSDCutoff() const
      {
        return m_RMSDCutoff;
      }

      //! @brief get seed length
      //! @return seed length
      size_t GetSeedLength() const
      {
        return m_SeedLength;
      }

      //! @brief get number of iterations
      //! @return number of iterations
      size_t GetNumberIterations() const
      {
        return m_NumberIterations;
      }

    ////////////////
    // operations //
    ////////////////

    private:

      //! @brief find a larger subset by extending the given one that has a RMSD below the cutoff
      //! @param SUBSET subset of coordinates that are used to be as seed and to be extended
      //! @param COORDINATES coordinates of interest
      //! @param REFERENCE_COORDINATES reference coordinates
      //! @return transformation for the largest subset
      math::TransformationMatrix3D ExtendSubset
      (
        storage::List< size_t> &SUBSET,
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
      ) const;

      //! @brief collects the subset of coordinates specified by the given list of indices
      //! @param SUBSET indices of coordinates that define the subset of coordinates to be collected
      //! @param COORDINATES util::SiPtrVector< Vector3D> which will be measured for agreement with COORDINATES_B
      //! @return vector of coordinates that correspond to the requested subset
      util::SiPtrVector< const linal::Vector3D> CollectCoordinatesSubset
      (
        const storage::List< size_t> &SUBSET,
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES
      ) const;

    ////////////////
    // operations //
    ////////////////

    public:

      //! @brief calculates MaxSub between given coordinates
      //! @param COORDINATES coordinates of interest
      //! @param REFERENCE_COORDINATES reference coordinates
      //! @return MaxSub between given coordinates
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

      //! @brief calculates the MaxSub and the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
      //! @param COORDINATES coordinates of interest
      //! @param REFERENCE_COORDINATES reference coordinates
      //! @return pair of the MaxSub and the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
      storage::Pair< double, math::TransformationMatrix3D> CalculateMaxSubAndSuperimposition
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

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class MaxSub

  } // namespace quality
} // namespace bcl

#endif // BCL_QUALITY_MAXSUB_H_ 
