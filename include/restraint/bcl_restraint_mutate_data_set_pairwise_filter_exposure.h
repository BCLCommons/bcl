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

#ifndef BCL_RESTRAINT_MUTATE_DATA_SET_PAIRWISE_FILTER_EXPOSURE_H_
#define BCL_RESTRAINT_MUTATE_DATA_SET_PAIRWISE_FILTER_EXPOSURE_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_atom_coordinates_interface.h"
#include "math/bcl_math_mutate_interface.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateDataSetPairwiseFilterExposure
    //! @brief Mutates a DataSetPairwise by removing data pairs that contain points not meeting an exposure criteria
    //! @details If either of the two data points do not have an exposure greater than the provided cutoff then they
    //!          are removed from the data set. The average exposure of the data point is used when an ensemble of
    //!          structures is provided. Currently neighbor count is used as the exposure measure, so data points
    //!          must have an exposure less than the provided threshold to be kept.
    //!
    //! @see @link example_restraint_mutate_data_set_pairwise_filter_exposure.cpp @endlink
    //! @author alexanns
    //! @date May 10, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateDataSetPairwiseFilterExposure :
      public math::MutateInterface< DataSetPairwise>
    {

    private:

      // typedef for exposure map which maps data point to the exposures calculated for it from an ensemble
      typedef storage::Map
      <
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>,
        storage::Vector< double>, //< exposure values, one per structure in ensemble
        assemble::LocatorAtomCoordinatesInterface::PtrLessThan
      > ExposureMap;

    //////////
    // data //
    //////////

      //! the amount of exposure a data point must have to not be removed
      double m_ExposureCutoff;

      //! maps a data point to the exposure measures calculated for it over the ensemble
      ExposureMap m_ExposureMap;

      //! the scheme of this mutate
      std::string m_Scheme;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor taking optional scheme
      //! @param SCHEME the scheme for this scoring function
      explicit MutateDataSetPairwiseFilterExposure( const std::string &SCHEME = GetDefaultScheme());

      //! @brief constructor taking members
      //! @param EXPOSURE_CUTOFF data involving residues outside of this exposure amount will be penalized
      //! @param ENSEMBLE ensemble for which exposures will be calculated
      //! @param SCHEME the scheme for this scoring function
      MutateDataSetPairwiseFilterExposure
      (
        const double &EXPOSURE_CUTOFF,
        const assemble::ProteinEnsemble &ENSEMBLE,
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief Clone function
      //! @return pointer to new MutateDataSetPairwiseFilterExposure
      MutateDataSetPairwiseFilterExposure *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking an DataSetPairwise and returning a mutated object
      //! @param DATA_SET DataSetPairwise of interest that will be mutated
      //! @return MutateResult that results from mutating to the argument DATA_SET
      math::MutateResult< DataSetPairwise> operator()( const DataSetPairwise &DATA_SET) const;

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

    }; // class MutateDataSetPairwiseFilterExposure

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_MUTATE_DATA_SET_PAIRWISE_FILTER_EXPOSURE_H_ 
