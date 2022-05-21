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
#include "score/bcl_score_data_set_pairwise_euclidian_distance.h"

// includes from bcl - sorted alphabetically
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
    const util::SiPtr< const util::ObjectInterface> DataSetPairwiseEuclidianDistance::s_Instance
    (
      GetObjectInstances().AddInstance( new DataSetPairwiseEuclidianDistance())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &DataSetPairwiseEuclidianDistance::GetDefaultScheme()
    {
      static const std::string s_scheme( "distance_range");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor taking optional scheme
    //! @param SCHEME the scheme for this scoring function
    DataSetPairwiseEuclidianDistance::DataSetPairwiseEuclidianDistance( const std::string &SCHEME) :
      m_DistanceRange(),
      m_Ensemble(),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking optional scheme
    //! @param RANGE desired range of euclidian separation between data points
    //! @param ENSEMBLE ensemble distances will be calculated from
    //! @param SCHEME the scheme for this scoring function
    DataSetPairwiseEuclidianDistance::DataSetPairwiseEuclidianDistance
    (
      const math::Range< double> &RANGE,
      const util::ShPtr< assemble::ProteinEnsemble> &ENSEMBLE,
      const std::string &SCHEME
    ) :
      m_DistanceRange( RANGE),
      m_Ensemble( ENSEMBLE),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking optional scheme
    //! @param RANGE desired range of euclidian separation between data points
    //! @param MODEL protein model distances will be calculated from
    //! @param SCHEME the scheme for this scoring function
    DataSetPairwiseEuclidianDistance::DataSetPairwiseEuclidianDistance
    (
      const math::Range< double> &RANGE,
      const assemble::ProteinModel &MODEL,
      const std::string &SCHEME
    ) :
      m_DistanceRange( RANGE),
      m_Ensemble(),
      m_Scheme( SCHEME)
    {
      // make ensemble from single model and set m_Ensemble to it
      util::ShPtr< assemble::ProteinEnsemble> ensemble( new assemble::ProteinEnsemble());
      ensemble->InsertElement( util::ShPtr< assemble::ProteinModel>( MODEL.Clone()));
      m_Ensemble = ensemble;
    }

    //! @brief Clone function
    //! @return pointer to new DataSetPairwiseEuclidianDistance
    DataSetPairwiseEuclidianDistance *DataSetPairwiseEuclidianDistance::Clone() const
    {
      return new DataSetPairwiseEuclidianDistance( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DataSetPairwiseEuclidianDistance::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &DataSetPairwiseEuclidianDistance::GetScheme() const
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
    double DataSetPairwiseEuclidianDistance::operator()( const restraint::DataSetPairwise &DATA) const
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
        // get mean distance - it is undefined if no counts were available
        const math::RunningAverageSD< double> mean_sd( data_itr->EuclidianDistance( *m_Ensemble).First());
        const double mean_distance( mean_sd.GetWeight() ? mean_sd.GetAverage() : util::GetUndefinedDouble());

        // true if size is within the desired range and defined
        if( m_DistanceRange.IsWithin( mean_distance) || !util::IsDefined( mean_distance))
        {
          continue;
        }

        // score is the distance outside of the range
        const double current_score
        (
          math::Absolute( mean_distance - m_DistanceRange.GetMiddle()) - ( m_DistanceRange.GetWidth() * 0.5)
        );

        score += current_score;
      }

      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DataSetPairwiseEuclidianDistance::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_DistanceRange, ISTREAM);
      io::Serialize::Read( m_Ensemble, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DataSetPairwiseEuclidianDistance::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_DistanceRange, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Ensemble, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace score
  
} // namespace bcl
