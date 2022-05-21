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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "score/bcl_score_data_set_pairwise_distance_change_magnitude.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "restraint/bcl_restraint_data_set_pairwise.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> DataSetPairwiseDistanceChangeMagnitude::s_Instance
    (
      GetObjectInstances().AddInstance( new DataSetPairwiseDistanceChangeMagnitude())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &DataSetPairwiseDistanceChangeMagnitude::GetDefaultScheme()
    {
      static const std::string s_scheme( "dist_change_magnitude");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor taking optional scheme
    //! @param SCHEME the scheme for this scoring function
    DataSetPairwiseDistanceChangeMagnitude::DataSetPairwiseDistanceChangeMagnitude( const std::string &SCHEME) :
      m_DistanceChanges(),
      m_MaxMagnitudeChange(),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking members
    //! @param ENSEMBLE_A ensemble for which distances will be calculated
    //! @param ENSEMBLE_B ensemble for which distances will be calculated
    //! @param FULL_DATA_SET data set distance changes should be calculated over
    //! @param SCHEME the scheme for this scoring function
    DataSetPairwiseDistanceChangeMagnitude::DataSetPairwiseDistanceChangeMagnitude
    (
      const assemble::ProteinEnsemble &ENSEMBLE_A,
      const assemble::ProteinEnsemble &ENSEMBLE_B,
      const restraint::DataSetPairwise &FULL_DATA_SET,
      const std::string &SCHEME
    ) :
      m_DistanceChanges(),
      m_MaxMagnitudeChange(),
      m_Scheme( SCHEME)
    {
      // get distance change statistics
      const storage::Pair
      <
        storage::Map< restraint::DataPairwise, math::RunningAverageSD< double> >, math::RunningMinMax< double>
      > stats
      (
        GetDistanceChangeStatistics( ENSEMBLE_A, ENSEMBLE_B, FULL_DATA_SET)
      );

      m_DistanceChanges = stats.First();

      // the max magnitude is the absolute value of either the minimum or maximum distance change
      m_MaxMagnitudeChange =
          std::max( math::Absolute( stats.Second().GetMax()), math::Absolute( stats.Second().GetMin()));
    }

    //! @brief Clone function
    //! @return pointer to new DataSetPairwiseDistanceChangeMagnitude
    DataSetPairwiseDistanceChangeMagnitude *DataSetPairwiseDistanceChangeMagnitude::Clone() const
    {
      return new DataSetPairwiseDistanceChangeMagnitude( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DataSetPairwiseDistanceChangeMagnitude::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &DataSetPairwiseDistanceChangeMagnitude::GetScheme() const
    {
      return m_Scheme;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate the score of a data set
    //! @param DATA data set to be scored
    //! @return the score of the current data set
    double DataSetPairwiseDistanceChangeMagnitude::operator()( const restraint::DataSetPairwise &DATA) const
    {
      double score( 0);

      // iterate through the data set
      for
      (
        restraint::DataSetPairwise::const_iterator data_itr( DATA.Begin()), data_itr_end( DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // find the pairwise data in m_DistanceChanges
        storage::Map< restraint::DataPairwise, math::RunningAverageSD< double> >::const_iterator
          found_data_itr( m_DistanceChanges.Find( *data_itr));

        // true if data was not found
        if( found_data_itr == m_DistanceChanges.End())
        {
          continue;
        }

        // get the current average distance change
        double mean_distance_change_magnitude( math::Absolute( found_data_itr->second.GetAverage()));

        // divide mean_distance_change by the maximum distance change
        mean_distance_change_magnitude /= m_MaxMagnitudeChange;

        // give bonus for changes that are larger fraction of m_MaxMagnitudeChange
        score -= mean_distance_change_magnitude;
      }

      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DataSetPairwiseDistanceChangeMagnitude::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_DistanceChanges, ISTREAM);
      io::Serialize::Read( m_MaxMagnitudeChange, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DataSetPairwiseDistanceChangeMagnitude::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_DistanceChanges, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MaxMagnitudeChange, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief for two ensembles calculates the mean and stddev distance changes observed for each data pair
    //! @param ENSEMBLE_A ensemble of first state of the protein
    //! @param ENSEMBLE_B ensemble of the second state of the protein
    //! @param DATA_SET the data set for which distance changes will be calculated
    //! @return pair of map with each datapair in DATA_SET and the associated mean and sd distance change observed in
    //!         the two ensembles
    //!         and a RunningMinMax< double> indicating the min and max mean distance change observed
    storage::Pair
    <
      storage::Map< restraint::DataPairwise, math::RunningAverageSD< double> >, math::RunningMinMax< double>
    > DataSetPairwiseDistanceChangeMagnitude::GetDistanceChangeStatistics
    (
      const assemble::ProteinEnsemble &ENSEMBLE_A, const assemble::ProteinEnsemble &ENSEMBLE_B,
      const restraint::DataSetPairwise &DATA_SET
    )
    {
      // hold mean and std dev distance change for each data pair as calculated from ENSEMBLE_A and ENSEMBLE_B
      storage::Map< restraint::DataPairwise, math::RunningAverageSD< double> > data_mean;

      // will hold the minimum and maximum distance change observed
      math::RunningMinMax< double> min_max;

      // iterate through the data set
      for
      (
        restraint::DataSetPairwise::const_iterator data_itr( DATA_SET.Begin()), data_itr_end( DATA_SET.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // get the distance change mean and standard deviation for the current data pair
        const math::RunningAverageSD< double> distance_changes_mean_sd
        (
          ENSEMBLE_A.GetDistanceChangesMeanSD( *data_itr, ENSEMBLE_B)
        );

        // make a pair out of the current data pair and its statistics object
        const std::pair< restraint::DataPairwise, math::RunningAverageSD< double> > pair
        (
          *data_itr, distance_changes_mean_sd
        );

        // insert the current data pair and its statistics object and make sure insertion was successful
        const bool successful_insertion( data_mean.Insert( pair).second);
        BCL_Assert( successful_insertion, "could not insert");

        // add the mean distance change to min max statistics object
        min_max += distance_changes_mean_sd.GetAverage();
      }

      // return the data of mean, stddev for each data pair and min, max ever observed
      const storage::Pair
       <
         storage::Map< restraint::DataPairwise, math::RunningAverageSD< double> >, math::RunningMinMax< double>
       > pair( data_mean, min_max);

      return pair;
    }

  } // namespace score

} // namespace bcl
