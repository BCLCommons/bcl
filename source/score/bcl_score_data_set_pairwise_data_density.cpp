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
#include "score/bcl_score_data_set_pairwise_data_density.h"

// includes from bcl - sorted alphabetically
#include "restraint/bcl_restraint_data_set_pairwise.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> DataSetPairwiseDataDensity::s_Instance
    (
      GetObjectInstances().AddInstance( new DataSetPairwiseDataDensity())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &DataSetPairwiseDataDensity::GetDefaultScheme()
    {
      static const std::string s_scheme( "data_density");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor taking optional scheme
    //! @param SCHEME the scheme for this scoring function
    DataSetPairwiseDataDensity::DataSetPairwiseDataDensity( const std::string &SCHEME) :
      m_SequenceLength(),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking scheme
    //! @param SEQUENCE_SIZE the sequence size the score will be normalized against
    //! @param SCHEME the scheme for this scoring function
    DataSetPairwiseDataDensity::DataSetPairwiseDataDensity
    (
      const size_t SEQUENCE_SIZE, const std::string &SCHEME
    ) :
      m_SequenceLength( SEQUENCE_SIZE),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new DataSetPairwiseDataDensity
    DataSetPairwiseDataDensity *DataSetPairwiseDataDensity::Clone() const
    {
      return new DataSetPairwiseDataDensity( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DataSetPairwiseDataDensity::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &DataSetPairwiseDataDensity::GetScheme() const
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
    double DataSetPairwiseDataDensity::operator()( const restraint::DataSetPairwise &DATA) const
    {
      // get the set of all unique individual data points
      util::ShPtrList< assemble::LocatorAtomCoordinatesInterface> data_points( DATA.GetDataPoints());

      // the number of unique data point positions (number of spin labels)
      const size_t l( data_points.GetSize());

      // the number of intervals
      const size_t n( l + 1);

      // optimal interval between spin labels
      const size_t I( m_SequenceLength / n);

      const double label_density_score( -CalculateDensityScore( data_points, I, n));

      return label_density_score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DataSetPairwiseDataDensity::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SequenceLength, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DataSetPairwiseDataDensity::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SequenceLength, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculates the density score
    //! @param DATA_POINTS the set of all unique individual data points
    //! @param INTERVAL optimal interval between individual data points, "I" in reference
    //! @param NUM_INTERVALS number of intervals between unique individual data points, "n" in reference
    //! @return double which is the density score
    double DataSetPairwiseDataDensity::CalculateDensityScore
    (
      const util::ShPtrList< assemble::LocatorAtomCoordinatesInterface> &DATA_POINTS,
      const size_t INTERVAL,
      const size_t NUM_INTERVALS
    ) const
    {
      double score( 0);

      if( DATA_POINTS.IsEmpty())
      {
        return 0;
      }

      BCL_Assert
      (
        DATA_POINTS.GetSize() > 1,
        "datapoints size should never be less than 2 at this point. The size is " + util::Format()( DATA_POINTS.GetSize())
      );

      // for nterminal interval
      {
        const size_t sequence_separation( ( *DATA_POINTS.Begin())->GetSeqID() - 1);
        const double current_score
        (
          math::Pow( double( math::Pow( sequence_separation - INTERVAL, size_t( 2)) + 1), double( -1))
        );
        score += current_score;
      }

      // for cterminal interval
      {
        const size_t sequence_separation( m_SequenceLength - ( *--DATA_POINTS.End())->GetSeqID());
        const double current_score
        (
          math::Pow( double( math::Pow( sequence_separation - INTERVAL, size_t( 2)) + 1), double( -1))
        );
        score += current_score;
      }

      // iterate through the data points two by two to sum up score
      for
      (
          util::ShPtrList< assemble::LocatorAtomCoordinatesInterface>::const_iterator
          data_itr( DATA_POINTS.Begin()),
          data_itr_b( ++DATA_POINTS.Begin()),
          data_itr_end( DATA_POINTS.End());
        data_itr != data_itr_end && data_itr_b != data_itr_end;
        ++data_itr, ++data_itr_b
      )
      {
        // get the current separation
        const size_t sequence_separation
        (
          restraint::DataPairwise::CalculateSequenceSeparation
          (
            ( *data_itr)->GetChainID(), ( *data_itr)->GetSeqID(),
            ( *data_itr_b)->GetChainID(), ( *data_itr_b)->GetSeqID()
          )
        );

        // true if the sequence separation is undefined
        if( !util::IsDefined( sequence_separation))
        {
          // continue to next pair
          continue;
        }

        const double current_score
        (
          math::Pow( double( math::Pow( sequence_separation - INTERVAL, size_t( 2)) + 1), double( -1))
        );

        score += current_score;
      }

      // divide by the number of intervals
      score /= double( NUM_INTERVALS);

      return score;
    }

  } // namespace score
  
} // namespace bcl
