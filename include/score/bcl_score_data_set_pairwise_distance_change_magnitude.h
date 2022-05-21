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

#ifndef BCL_SCORE_DATA_SET_PAIRWISE_DISTANCE_CHANGE_MAGNITUDE_H_
#define BCL_SCORE_DATA_SET_PAIRWISE_DISTANCE_CHANGE_MAGNITUDE_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "restraint/bcl_restraint.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_function_interface_serializable.h"
#include "math/bcl_math_running_average_sd.h"
#include "restraint/bcl_restraint_data_pairwise.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DataSetPairwiseDistanceChangeMagnitude
    //! @brief Scores a DataSetPairwise to select for datapairs that give large distance changes
    //! @details Given two ensembles of two states that differ (i.e. different conformations), the score will favor
    //!          data pairs that give larger distance changes according to how the distance are changing in the
    //!          ensembles.
    //!
    //! @see @link example_score_data_set_pairwise_distance_change_magnitude.cpp @endlink
    //! @author alexanns
    //! @date May 10, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DataSetPairwiseDistanceChangeMagnitude :
      public math::FunctionInterfaceSerializable< restraint::DataSetPairwise, double>
    {

    private:

    //////////
    // data //
    //////////

      //! holds the stat for the distance change associated each datapair
      storage::Map< restraint::DataPairwise, math::RunningAverageSD< double> > m_DistanceChanges;

      //! the largest magnitude change seen between the two states
      double m_MaxMagnitudeChange;

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
      explicit DataSetPairwiseDistanceChangeMagnitude( const std::string &SCHEME = GetDefaultScheme());

      //! @brief constructor taking members
      //! @param ENSEMBLE_A ensemble for which distances will be calculated
      //! @param ENSEMBLE_B ensemble for which distances will be calculated
      //! @param FULL_DATA_SET data set distance changes should be calculated over
      //! @param SCHEME the scheme for this scoring function
      DataSetPairwiseDistanceChangeMagnitude
      (
        const assemble::ProteinEnsemble &ENSEMBLE_A,
        const assemble::ProteinEnsemble &ENSEMBLE_B,
        const restraint::DataSetPairwise &FULL_DATA_SET,
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief Clone function
      //! @return pointer to new DataSetPairwiseDistanceChangeMagnitude
      DataSetPairwiseDistanceChangeMagnitude *Clone() const;

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

    private:

      //! @brief for two ensembles calculates the mean and stddev distance changes observed for each data pair
      //! @param ENSEMBLE_A ensemble of first state of the protein
      //! @param ENSEMBLE_B ensemble of the second state of the protein
      //! @param DATA_SET the data set for which distance changes will be calculated
      //! @return pair of map with each datapair in DATA_SET and the associated mean and sd distance change observed in
      //!         the two ensembles
      //!         and a RunningMinMax< double> indicating the min and max mean distance change observed
      static storage::Pair
      <
        storage::Map< restraint::DataPairwise, math::RunningAverageSD< double> >, math::RunningMinMax< double>
      > GetDistanceChangeStatistics
      (
        const assemble::ProteinEnsemble &ENSEMBLE_A, const assemble::ProteinEnsemble &ENSEMBLE_B,
        const restraint::DataSetPairwise &DATA_SET
      );

    }; // class DataSetPairwiseDistanceChangeMagnitude

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_DATA_SET_PAIRWISE_DISTANCE_CHANGE_MAGNITUDE_H_ 
