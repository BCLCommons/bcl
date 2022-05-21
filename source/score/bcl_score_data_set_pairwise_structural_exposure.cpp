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
#include "score/bcl_score_data_set_pairwise_structural_exposure.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_count.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "math/bcl_math_statistics.h"
#include "math/bcl_math_trigonometric_transition.h"
#include "restraint/bcl_restraint_data_set_pairwise.h"
#include "restraint/bcl_restraint_locator_coordinates_first_side_chain_atom.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> DataSetPairwiseStructuralExposure::s_Instance
    (
      GetObjectInstances().AddInstance( new DataSetPairwiseStructuralExposure())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &DataSetPairwiseStructuralExposure::GetDefaultScheme()
    {
      static const std::string s_scheme( "exposure");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor taking optional scheme
    //! @param SCHEME the scheme for this scoring function
    DataSetPairwiseStructuralExposure::DataSetPairwiseStructuralExposure( const std::string &SCHEME) :
      m_ExposureCutoff(),
      m_ExposureMap(),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking members
    //! @param EXPOSURE_CUTOFF data involving residues outside of this exposure amount will be penalized
    //! @param ENSEMBLE ensemble for which exposures will be calculated
    //! @param SCHEME the scheme for this scoring function
    DataSetPairwiseStructuralExposure::DataSetPairwiseStructuralExposure
    (
      const double &EXPOSURE_CUTOFF,
      const assemble::ProteinEnsemble &ENSEMBLE,
      const std::string &SCHEME
    ) :
      m_ExposureCutoff( EXPOSURE_CUTOFF),
      m_ExposureMap(),
      m_Scheme( SCHEME)
    {
      FillExposureMap( ENSEMBLE, m_ExposureMap);
    }

    //! @brief Clone function
    //! @return pointer to new DataSetPairwiseStructuralExposure
    DataSetPairwiseStructuralExposure *DataSetPairwiseStructuralExposure::Clone() const
    {
      return new DataSetPairwiseStructuralExposure( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DataSetPairwiseStructuralExposure::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &DataSetPairwiseStructuralExposure::GetScheme() const
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
    double DataSetPairwiseStructuralExposure::operator()( const restraint::DataSetPairwise &DATA) const
    {
      if( DATA.IsEmpty())
      {
        return 0;
      }

      double score( 0);

      // iterate through the data set
      for
      (
        restraint::DataSetPairwise::const_iterator data_itr( DATA.Begin()), data_itr_end( DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        const double score_first( CalculateExposureScore( data_itr->First(), m_ExposureMap, m_ExposureCutoff));
        const double score_second( CalculateExposureScore( data_itr->Second(), m_ExposureMap, m_ExposureCutoff));

        // true if the first data point is defined
        if( util::IsDefined( score_first))
        {
          score += score_first;
        }

        // true if the second data point is defined
        if( util::IsDefined( score_second))
        {
          score += score_second;
        }
      }

      score /= double( DATA.GetSize());

      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DataSetPairwiseStructuralExposure::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ExposureCutoff, ISTREAM);
      io::Serialize::Read( m_ExposureMap, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DataSetPairwiseStructuralExposure::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ExposureCutoff, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ExposureMap, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief fills the given exposure map based on the exposures calculated from a given ensemble
    //! @param ENSEMBLE ensemble of models exposures will be calculated for
    //! @param EXPOSURE_MAP the map of exposures that will be filled
    void DataSetPairwiseStructuralExposure::FillExposureMap
    (
      const assemble::ProteinEnsemble &ENSEMBLE, ExposureMap &EXPOSURE_MAP
    )
    {
      assemble::AANeighborCount neighbor_counter;

      // iterate through the ensemble
      for
      (
        assemble::ProteinEnsemble::const_iterator ensemble_itr( ENSEMBLE.Begin()), ensemble_itr_end( ENSEMBLE.End());
        ensemble_itr != ensemble_itr_end;
        ++ensemble_itr
      )
      {
        // get the neighbor list container for the current protein
        const assemble::AANeighborListContainer nl_container
        (
          ( *ensemble_itr)->GetAminoAcids(),
          neighbor_counter.GetDistanceCutoff(),
          neighbor_counter.GetMinimalSequenceSeparation(),
          true
        );

        // iterate through the neighbor list container
        for
        (
          assemble::AANeighborListContainer::const_iterator
            nl_itr( nl_container.Begin()), nl_itr_end( nl_container.End());
          nl_itr != nl_itr_end;
          ++nl_itr
        )
        {
          // get the center residue
          const util::SiPtr< const biol::AABase> &center_aa( nl_itr->second.GetCenterAminoAcid());

          BCL_Assert
          (
            center_aa != util::SiPtr< const biol::AABase>( assemble::AANeighborList::GetDefaultCenterAA()),
            "center residue undefined"
          );

          // make locator to current center residue
          const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> residue_locator
          (
            new restraint::LocatorCoordinatesFirstSideChainAtom
            (
              center_aa->GetChainID(),                      center_aa->GetSeqID(),
              center_aa->GetFirstSidechainAtom().GetType(), center_aa->GetType()
            )
          );

          BCL_MessageDbg
          (
            "doing nc for residue " + nl_itr->second.GetCenterAminoAcid()->GetIdentification() +
            " center_aa identification is " + center_aa->GetIdentification()
          );

          // get the neighbor count for the current residue
          const double nc( neighbor_counter( nl_itr->second));

          BCL_MessageDbg
          (
            "current neighbor count is " + util::Format()( nc) + " for " +
            residue_locator->GetIdentification()
          );

          // try to find the residue locator
          ExposureMap::iterator resi_locator_itr( EXPOSURE_MAP.Find( residue_locator));

          // true if locator was not found
          if( resi_locator_itr == EXPOSURE_MAP.End())
          {
            // insert the locator with a vector containing the current nc
            std::pair< util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, storage::Vector< double> > pair
            (
              residue_locator, storage::Vector< double>( 1, nc)
            );
            BCL_MessageDbg
            (
              "inserted residue key for " + residue_locator->GetIdentification() +
              " with nc " + util::Format()( nc)
            );

            BCL_Assert
            (
              EXPOSURE_MAP.Insert( pair).second, "could not insert " + residue_locator->GetIdentification() +
              " with nc " + util::Format()( nc)
            );
          }
          else //< residue locator already exists in map so just add nc to the vector
          {
            // add the current neighbor count to the list of neighbor counts for the current residue
            resi_locator_itr->second.PushBack( nc);
            BCL_MessageDbg
            (
              "inserted nc value for " + residue_locator->GetIdentification() +
              " with nc " + util::Format()( nc)
            );
          }
        } //< iterate through current neighbor list container
      } //< iterate through ensemble
    } //< function FillExposureMap

    //! @brief calculates the score that a given data point has given its exposures and the exposure cutoff
    //! @param LOCATOR the data point the score will be calculated for
    //! @param EXPOSURE_MAP the map of exposures
    //! @param EXPOSURE_CUTOFF data involving residues outside of this exposure amount will be penalized
    //! @return double which is the score of LOCATOR
    double DataSetPairwiseStructuralExposure::CalculateExposureScore
    (
      const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &LOCATOR,
      const ExposureMap &EXPOSURE_MAP,
      const double EXPOSURE_CUTOFF
    )
    {
      // try to find the datapoint in the exposure map
      ExposureMap::const_iterator itr_first( EXPOSURE_MAP.Find( LOCATOR));

      // true if the data point is not found in the exposure map
      if( itr_first == EXPOSURE_MAP.End())
      {
        BCL_MessageDbg( "could not find " + LOCATOR->GetIdentification());
        return util::GetUndefinedDouble();
      }

      // get the average exposure for the current locator
      const double mean_exposure( math::Statistics::Mean( itr_first->second.Begin(), itr_first->second.End()));

      BCL_MessageDbg
      (
        "mean exposure for " + LOCATOR->GetIdentification() + " is " + util::Format()( mean_exposure)
      );

      // percent the mean exposure is of the desired exposure cutoff
      const double percent_exposure( ( mean_exposure - EXPOSURE_CUTOFF) / EXPOSURE_CUTOFF);

      // threshold for how much larger the mean exposure can be than the desired exposure cutoff
      // ( as a % of desired cutoff) before maximum penalty is given
      static const double transition_amount( 0.5);

      // the max score i.e. worst score
      static const double max_penalty( 1.0);

      // the most favorable score
      static const double best_score( 0);

      // true if mean exposure is more than allowed by transition
      if( percent_exposure >= transition_amount)
      {
        // return max penalty score of max_penalty
        return max_penalty;
      }

      // true if mean exposure is within the transition region
      if( percent_exposure < transition_amount && mean_exposure > EXPOSURE_CUTOFF)
      {
        const double score
        (
          math::TrigonometricTransition
          (
            EXPOSURE_CUTOFF, //< x value of start of transition
            EXPOSURE_CUTOFF + transition_amount * EXPOSURE_CUTOFF, //< x value at end of transition
            best_score, //< score value at start of transition region
            max_penalty //< score value at end of transition region
          )
          ( mean_exposure) //< current value
        );

        return score;
      }

      // at this point best score can be returned
      return best_score;
    }

  } // namespace score

} // namespace bcl
