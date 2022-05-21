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

#ifndef BCL_SCORE_DATA_SET_PAIRWISE_STRUCTURAL_EXPOSURE_H_
#define BCL_SCORE_DATA_SET_PAIRWISE_STRUCTURAL_EXPOSURE_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "restraint/bcl_restraint.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_atom_coordinates_interface.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DataSetPairwiseStructuralExposure
    //! @brief Favors DataSetPairwise with data points whose exposure meets a given desired minimum exposure
    //! @details If both of the data points in the data pair have average (if ensemble) exposure greater than the
    //!          provided threshold, then the data pair is considered favorable.
    //!
    //! @see @link example_score_data_set_pairwise_structural_exposure.cpp @endlink
    //! @author alexanns
    //! @date May 8, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DataSetPairwiseStructuralExposure :
      public math::FunctionInterfaceSerializable< restraint::DataSetPairwise, double>
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
      explicit DataSetPairwiseStructuralExposure( const std::string &SCHEME = GetDefaultScheme());

      //! @brief constructor taking members
      //! @param EXPOSURE_CUTOFF data involving residues outside of this exposure amount will be penalized
      //! @param ENSEMBLE ensemble for which exposures will be calculated
      //! @param SCHEME the scheme for this scoring function
      DataSetPairwiseStructuralExposure
      (
        const double &EXPOSURE_CUTOFF,
        const assemble::ProteinEnsemble &ENSEMBLE,
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief Clone function
      //! @return pointer to new DataSetPairwiseStructuralExposure
      DataSetPairwiseStructuralExposure *Clone() const;

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

      //! @brief calculate the score of a data set
      //! @param DATA data set to be scored
      //! @return the score of the current data set
      double operator()( const restraint::DataSetPairwise &DATA) const;

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

      //! @brief fills the given exposure map based on the exposures calculated from a given ensemble
      //! @param ENSEMBLE ensemble of models exposures will be calculated for
      //! @param EXPOSURE_MAP the map of exposures that will be filled
      static void FillExposureMap( const assemble::ProteinEnsemble &ENSEMBLE, ExposureMap &EXPOSURE_MAP);

      //! @brief calculates the score that a given data point has given its exposures and the exposure cutoff
      //! @param LOCATOR the data point the score will be calculated for
      //! @param EXPOSURE_MAP the map of exposures
      //! @param EXPOSURE_CUTOFF data involving residues outside of this exposure amount will be penalized
      //! @return double which is the score of LOCATOR
      static double CalculateExposureScore
      (
        const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &LOCATOR,
        const ExposureMap &EXPOSURE_MAP,
        const double EXPOSURE_CUTOFF
      );

    }; // class DataSetPairwiseStructuralExposure

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_DATA_SET_PAIRWISE_STRUCTURAL_EXPOSURE_H_ 
