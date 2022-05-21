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

#ifndef BCL_QUALITY_LCS_H_
#define BCL_QUALITY_LCS_H_

// include the namespace header
#include "bcl_quality.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_quality_superimpose_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LCS
    //! @brief LCS (longest continuous segment) returns the length of the longest continuous segment that can be
    //!        superimposed below given cutoff RMSD
    //! @details LCS algorithm starts with set of overlapping segments of given seed length and tries to elongate them
    //!          while making sure when superimposed the distance for the segment is still below the given RMSD cutoff.
    //!
    //! @see Zemla,A. (2003) LGA—a method for finding 3D similarities in protein structures. Nucleic Acids Res., 31, 3370–3374
    //! @see @link example_quality_lcs.cpp @endlink
    //! @author karakam
    //! @date Oct 8, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LCS :
      public SuperimposeInterface
    {

    private:

    //////////
    // data //
    //////////

      //! RMSD cutoff
      double m_RmsdCutoff;

      //! seed length
      size_t m_SeedLength;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief returns the default Rmsd cutoff
      //! @return the default Rmsd cutoff
      static double GetDefaultRmsdCutoff();

      //! @brief returns the default seed length
      //! @return the default seed length
      static size_t GetDefaultSeedLength();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from a RMSD cutoff and a seed length
      //! @param RMSD_CUTOFF distance cutoff
      //! @param SEED_LENGTH length of seeds
      LCS
      (
        const double RMSD_CUTOFF = GetDefaultRmsdCutoff(),
        const size_t SEED_LENGTH = GetDefaultSeedLength()
      );

      //! @brief Clone function
      //! @return pointer to new LCS
      LCS *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return the optimal value for that quality measurement
      //! @return the best value by which two sets of coordinates can agree
      double OptimalValue() const;

      //! @brief return the comparison function for better quality
      //! @return binary function to compare two quality measure values
      const util::BinaryFunctionInterface< double, double, bool> &GetComparisonFunction() const;

      //! @brief return rmsd cutoff
      //! @return rmsd cutoff
      double GetCutoff() const;

      //! @brief get seed length
      //! @return seed length
      size_t GetSeedLength() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief find a larger range by extending the given one that has a RMSD below the cutoff
      //! @param RANGE range of coordinates that are used to be as seed and to be extended
      //! @param COORDINATES coordinates of interest
      //! @param REFERENCE_COORDINATES reference coordinates
      //! @return extended range
      math::Range< size_t> ExtendRange
      (
        const math::Range< size_t> &RANGE,
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
      ) const;

      //! @brief returns the ranges of longest continuous segments that can be superimposed below cutoff for given coordinates
      //! @param COORDINATES coordinates of interest
      //! @param REFERENCE_COORDINATES reference coordinates
      //! @return the ranges of longest continuous segments that can be superimposed below cutoff for given coordinates
      storage::List< math::Range< size_t> > CalculateRanges
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
      ) const;

      //! @brief returns the indices to the coordinates of longest continuous segments that can be superimposed below cutoff for given coordinates
      //! @param COORDINATES coordinates of interest
      //! @param REFERENCE_COORDINATES reference coordinates
      //! @return the indices of longest continuous segments that can be superimposed below cutoff for given coordinates
      storage::List< storage::List< size_t> > CalculateIndices
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
      ) const;

      //! @brief calculates LCS between given coordinates
      //! @param COORDINATES coordinates of interest
      //! @param REFERENCE_COORDINATES reference coordinates
      //! @return LCS between COORDINATES and REFERENCE_COORDINATES
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

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief check if a given range superimposes below the cutoff
      //! @param COORDINATES coordinates of interest
      //! @param REFERENCE_COORDINATES reference coordinates
      //! @return true, if coordinates within the given range are superimposable below the cutoff
      bool IsGoodRange
      (
        const math::Range< size_t> &RANGE,
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
      ) const;

      //! @brief converts a range to list of indices
      //! @param RANGE range of indices
      //! @return list that contains the indices in the range
      static storage::List< size_t> ConvertRangeToIndices
      (
        const math::Range< size_t> &RANGE
      );

      //! @brief converts a vector of ranges to list of list of indices
      //! @param RANGES vector of range of indices
      //! @return vector of lists that contains the indices in the range vector
      static storage::List< storage::List< size_t> > ConvertRangesToLists
      (
        const storage::List< math::Range< size_t> > &RANGES
      );

    }; // class LCS

  } // namespace quality
} // namespace bcl

#endif // BCL_QUALITY_LCS_H_ 
