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

#ifndef BCL_QUALITY_GDT_H_
#define BCL_QUALITY_GDT_H_

// include the namespace header
#include "bcl_quality.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"
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
    //! @class GDT
    //! @brief GDT (Global Distance Test) is a measure to compare two sets of coordinates
    //! @details GDT provides an alternative coordinate comparison measure to RMSD. Unlike RMSD, it focuses on superimposing the
    //! best alignable parts of the coordinate sets, instead of trying to align all coordinates. This can be useful
    //! in protein structure-structure comparison to identify the best superimposable regions.
    //! The algorithm tries to find the largest subset of coordinate pairs that can be superimposed below a given RMSD
    //! cutoff. The alignment is not required to be continuous.
    //! This implementation also allows calculation of GDT_TS and GDT_HA, which are the averages of 4 different cutoffs
    //! 1.0, 2.0, 4.0, 8.0 Angstrom RMSD and 0.5, 1.0, 2.0, 4.0 Angstrom RMSD respectively.
    //! The iterative process starts with a set of seed fragments, short continuous stretches and with each iteration
    //! uses the previous superimposition to extend the subset of coordinate pairs that can be superimposed below the
    //! given cutoff.
    //! In addition to consecutive seeds, GDT also uses LCS to find longest continuous segments and also use these as
    //! additional seeds.
    //!
    //! @see Zemla,A. (2003) LGA—a method for finding 3D similarities in protein structures. Nucleic Acids Res., 31, 3370–3374
    //! @see http://predictioncenter.org/casp/casp7/public/doc/LCS_GDT.README
    //! @see @link example_quality_gdt.cpp @endlink
    //! @author karakam, woetzen
    //! @date Oct 8, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API GDT :
      public SuperimposeInterface
    {

    private:

    //////////
    // data //
    //////////

      //! distance cutoff
      double m_DistanceCutoff;

      //! length of seed
      size_t m_SeedLength;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief returns the default seed length
      //! @return the default seed length
      static const size_t GetDefaultSeedLength();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from a single distance cutoff and seed length
      //! @param DISTANCE_CUTOFF distance cutoff to be used
      //! @param SEED_LENGTH length of seed
      GDT
      (
        const double &DISTANCE_CUTOFF,
        const size_t SEED_LENGTH = GetDefaultSeedLength()
      );

      //! @brief Clone function
      //! @return pointer to new GDT
      GDT *Clone() const;

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

      //! @brief get seed length
      //! @return seed length
      size_t GetSeedLength() const;

      //! @brief get distance cutoff
      //! @return distance cutoff
      double GetDistanceCutoff() const;

    ////////////////
    // operations //
    ////////////////

    private:

      //! @brief collects the subset of coordinates specified by the given list of indices
      //! @param SUBSET indices of coordinates that define the subset of coordinates to be collected
      //! @param COORDINATES coordinates of interest
      //! @return vector of coordinates that correspond to the requested subset
      static util::SiPtrVector< const linal::Vector3D> CollectCoordinatesSubset
      (
        const storage::List< size_t> &SUBSET,
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES
      );

    ////////////////
    // operations //
    ////////////////

    public:

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

      //! @brief create an average GDT object with given set of distance cutoffs
      //! @param DISTANCE_CUTOFFS set of distance cutoffs
      //! @param SEED_LENGTH length of seed
      //! @return Average object
      static Average CreateAverageGDT
      (
        const storage::Set< double> &DISTANCE_CUTOFFS,
        const size_t SEED_LENGTH = GetDefaultSeedLength()
      );

    private:

      //! @brief the iterative superimposition procedure (ISP)
      //! @param COORDINATES coordinates of interest
      //! @param REFERENCE_COORDINATES reference coordinates
      //! @param START_INDICES the selected coordinates as seed, do not have to fullfill the CUTOFF
      //! @return the newly selected coordinates and its transformation
      storage::Pair< storage::List< size_t>, math::TransformationMatrix3D> IterativeSuperimpositionProcedure
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES,
        const storage::List< size_t> &START_INDICES
      ) const;

      //! @brief returns all possible seed ranges (not necessarily valid)
      //! @param NUMBER_COORDIANTES number of coordinates
      //! @return continuous seed fragments
      storage::List< math::Range< size_t> > GetSeedFragments
      (
        const size_t NUMBER_COORDINATES
      ) const;

      //! @brief check if two coordinates and the transformation will bring them close below cutoff
      //! @param COORDINATE coord to transform
      //! @param REFERENCE_COORDINATE coord to compare against
      //! @param TRANSFORMATION transformation to apply
      //! @param DISTANCE_CUTOFF_SQUARED the distance below which they are considered good
      //! @return true if the distance between the transformed coord and the reference is below cutoff
      static bool IsBelowCutoff
      (
        const linal::Vector3D &COORDINATE,
        const linal::Vector3D &REFERENCE_COORDINATE,
        const math::TransformationMatrix3D &TRANSFORMATION,
        const double DISTANCE_CUTOFF_SQUARED
      );

      //! @brief collect all indices below the cutoff
      //! @param COORDINATES coordinates of interest
      //! @param REFERENCE_COORDINATES reference coordinates
      //! @param TRANSFORMATION transformation to apply before comparison
      //! @param DISTANCE_CUTOFF distance cutoff
      //! @param list of coordinate indices that are below the cutoff for that transformation
      storage::List< size_t> CollectIndicesBelowCutoff
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES,
        const math::TransformationMatrix3D &TRANSFORMATION
      ) const;

      //! @brief transformation for a new set of indices
      //! @param COORDINATES coordinates of interest
      //! @param REFERENCE_COORDINATES reference coordinates
      //! @param INDICES the selected coordinates
      //! @return the transformation matrix to superimpose the coordinates selected by the indices
      static math::TransformationMatrix3D CalculateTransformation
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES,
        const storage::List< size_t> &INDICES
      );

    }; // class GDT

  } // namespace quality
} // namespace bcl

#endif // BCL_QUALITY_GDT_H_ 
