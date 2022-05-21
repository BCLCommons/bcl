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

#ifndef BCL_ASSEMBLE_AA_SASA_OLS_H_
#define BCL_ASSEMBLE_AA_SASA_OLS_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_assemble_aa_exposure_interface.h"
#include "math/bcl_math_range.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AASasaOLS
    //! @brief This is a AAExposureInterface derived class for determining the overlapping spheres (ols) solven accessible surface area for an AA
    //! @details Using the AANeighborList uses that overlapping spheres sasa algorithm .
    //!
    //! @see @link example_assemble_aa_sasa_ols.cpp @endlink
    //! @author woetzen
    //! @date 12.10.2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AASasaOLS :
      public AAExposureInterface
    {

    private:

    //////////
    // data //
    //////////

      //! scheme to be used in outputting
      std::string m_Scheme;

      //! the unified radius for spheres around the cb atom
      math::Range< double> m_SphereRadiusThreshold;

      //! @brief minimal sequence separation
      size_t m_MininalSequenceSeparation;

    public:

    //////////
    // data //
    //////////

      //! @brief returns the default sphere radius
      //! @return the default sphere radius
      static double GetDefaultSphereRadius();

      //! @brief default minimal sequence separation
      static const size_t s_DefaultMinimalSequenceSeparation = 0;

      //! @brief returns default file where the statistics and in consequence the energy potentials are read from
      //! @return default file where the statistics and in consequence the energy potentials are read from
      static const std::string &GetDefaultHistogramFilename();

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! @param SCHEME scheme to be used
      AASasaOLS( const std::string &SCHEME = GetDefaultScheme());

      //! @brief virtual copy constructor
      //! @return pointer to a new AASasaOLS object copied from this one
      AASasaOLS *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

      //! @brief min and max exposure measure
      //! @return range in which exposure can be
      const math::Range< double> &GetRange() const;

      //! @brief get threshold range
      //! @details threshold are used, to have a continuous function for the exposure measure, instead of a stepwise
      //! @return the default thresholds used
      const math::Range< double> &GetThresholdRange() const
      {
        return m_SphereRadiusThreshold;
      }

      //! @brief set the threshold range
      //! @param RANGE the range in which other neighbors that are considered for the exposture do not count full
      void SetThresholdRange( const math::Range< double> &RANGE)
      {
        m_SphereRadiusThreshold = RANGE;
      }

      //! @brief access to the minimal sequence separation
      //! @return minimal sequence separation in sequence distance
      size_t GetMinimalSequenceSeparation() const
      {
        return m_MininalSequenceSeparation;
      }

      //! @brief access to the minimal sequence separation
      //! @param MINIMAL_SEQUENCE_SEPARATION in sequence distance
      void SetMinimalSequenceSeparation( const size_t MINIMAL_SEQUENCE_SEPARATION)
      {
        m_MininalSequenceSeparation = MINIMAL_SEQUENCE_SEPARATION;
      }

      //! @brief return the histogram filename where statistics are stored
      //! @return the filename containing the thresholds, sequence separation and histograms for each environment and
      //!         aa type from a databank of proteins
      const std::string &GetHistogramFileName() const
      {
        return GetDefaultHistogramFilename();
      }

      //! @brief is direct measure - exposed surface correlates with the measure or actually measures buried surface
      //! @return true, if exposure measure correlates with exposed surface, false if it correlates to buried surface
      bool IsDirect() const
      {
        return true;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief calculate SASA for a given amino acid and its AANeighborList
      //! @return SASA for a given amino acid and its AANeighborList
      double operator()( const AANeighborList &AA_NEIGHBOR_LIST) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from istream
      //! @param ISTREAM is the input stream
      //! @return returns the input stream
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to ostream
      //! @param OSTREAM is the output stream
      //! @param INDENT number of indentations
      //! @return returns the output stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; //class AASasaOLS

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_AA_SASA_OLS_H_
