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
#include "restraint/bcl_restraint_mutate_data_set_pairwise_filter_exposure.h"

// includes from bcl - sorted alphabetically
#include "restraint/bcl_restraint_data_set_pairwise.h"
#include "score/bcl_score_data_set_pairwise_structural_exposure.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateDataSetPairwiseFilterExposure::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateDataSetPairwiseFilterExposure())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &MutateDataSetPairwiseFilterExposure::GetDefaultScheme()
    {
      static const std::string s_scheme( "filter_exposure");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor taking optional scheme
    //! @param SCHEME the scheme for this scoring function
    MutateDataSetPairwiseFilterExposure::MutateDataSetPairwiseFilterExposure( const std::string &SCHEME) :
      m_ExposureCutoff(),
      m_ExposureMap(),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking members
    //! @param EXPOSURE_CUTOFF data involving residues outside of this exposure amount will be penalized
    //! @param ENSEMBLE ensemble for which exposures will be calculated
    //! @param SCHEME the scheme for this scoring function
    MutateDataSetPairwiseFilterExposure::MutateDataSetPairwiseFilterExposure
    (
      const double &EXPOSURE_CUTOFF,
      const assemble::ProteinEnsemble &ENSEMBLE,
      const std::string &SCHEME
    ) :
      m_ExposureCutoff( EXPOSURE_CUTOFF),
      m_ExposureMap(),
      m_Scheme( SCHEME)
    {
      score::DataSetPairwiseStructuralExposure::FillExposureMap( ENSEMBLE, m_ExposureMap);
    }

    //! @brief Clone function
    //! @return pointer to new MutateDataSetPairwiseFilterExposure
    MutateDataSetPairwiseFilterExposure *MutateDataSetPairwiseFilterExposure::Clone() const
    {
      return new MutateDataSetPairwiseFilterExposure( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateDataSetPairwiseFilterExposure::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &MutateDataSetPairwiseFilterExposure::GetScheme() const
    {
      return m_Scheme;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an DataSetPairwise and returning a mutated object
    //! @param DATA_SET DataSetPairwise of interest that will be mutated
    //! @return MutateResult that results from mutating to the argument DATA_SET
    math::MutateResult< DataSetPairwise>
    MutateDataSetPairwiseFilterExposure::operator()( const DataSetPairwise &DATA_SET) const
    {
      static util::ShPtr< DataSetPairwise> s_empty;

      if( DATA_SET.IsEmpty())
      {
        BCL_MessageStd( "MutateDataSetPairwiseFilterExposure data set is empty");
        return math::MutateResult< DataSetPairwise>( s_empty, *this);
      }

      // the dataset filtered
      util::ShPtr< DataSetPairwise> filtered_data( new DataSetPairwise());

      // iterate through the data set to filter out data pairs that don't meet exposure criteria
      for
      (
        DataSetPairwise::const_iterator data_itr( DATA_SET.Begin()), data_itr_end( DATA_SET.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // score of first data point
        const double score_first
        (
          score::DataSetPairwiseStructuralExposure::CalculateExposureScore
          (
            data_itr->First(), m_ExposureMap, m_ExposureCutoff
          )
        );

        // score of second data point
        const double score_second
        (
          score::DataSetPairwiseStructuralExposure::CalculateExposureScore
          (
            data_itr->Second(), m_ExposureMap, m_ExposureCutoff
          )
        );

        BCL_MessageDbg
        (
          "first score " + util::Format()( score_first) +
          " second score " + util::Format()( score_second) + " for " + data_itr->GetIdentification()
        );

        // true if the first and second data points are defined and score is zero meaning they meet the cutoff
        if( util::IsDefined( score_first) && !score_first && util::IsDefined( score_second) && !score_second)
        {
          // add data pair to filtered_data
          filtered_data->Insert( *data_itr);
        }
      }

      BCL_MessageDbg
      (
        "MutateDataSetPairwiseFilterExposure num unique locators " +
        util::Format()( filtered_data->GetUniqueDataPoints().GetSize())
      );

      return math::MutateResult< DataSetPairwise>( filtered_data, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateDataSetPairwiseFilterExposure::Read( std::istream &ISTREAM)
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
    std::ostream &MutateDataSetPairwiseFilterExposure::Write( std::ostream &OSTREAM, const size_t INDENT) const
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

  } // namespace restraint
  
} // namespace bcl
