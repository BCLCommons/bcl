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
#include "score/bcl_score_data_set_pairwise_sequence_separation.h"

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
    const util::SiPtr< const util::ObjectInterface> DataSetPairwiseSequenceSeparation::s_Instance
    (
      GetObjectInstances().AddInstance( new DataSetPairwiseSequenceSeparation())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &DataSetPairwiseSequenceSeparation::GetDefaultScheme()
    {
      static const std::string s_scheme( "seq_sep");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor taking optional scheme
    //! @param SCHEME the scheme for this scoring function
    DataSetPairwiseSequenceSeparation::DataSetPairwiseSequenceSeparation( const std::string &SCHEME) :
      m_SequenceLength(),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking scheme
    //! @param SEQUENCE_SIZE the sequence size the score will be normalized against
    //! @param SCHEME the scheme for this scoring function
    DataSetPairwiseSequenceSeparation::DataSetPairwiseSequenceSeparation
    (
      const size_t SEQUENCE_SIZE, const std::string &SCHEME
    ) :
      m_SequenceLength( SEQUENCE_SIZE),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new DataSetPairwiseSequenceSeparation
    DataSetPairwiseSequenceSeparation *DataSetPairwiseSequenceSeparation::Clone() const
    {
      return new DataSetPairwiseSequenceSeparation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DataSetPairwiseSequenceSeparation::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &DataSetPairwiseSequenceSeparation::GetScheme() const
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
    double DataSetPairwiseSequenceSeparation::operator()( const restraint::DataSetPairwise &DATA) const
    {
      double seq_sep_sum( 0);

      if( DATA.IsEmpty())
      {
        return seq_sep_sum;
      }

      // will keep track of the number of data pairs whose sequence separations can be calculated
      // sequence separation is undefined for pairs on different chains
      size_t count( 0);

      // iterate through the data set
      for
      (
        restraint::DataSetPairwise::const_iterator data_itr( DATA.Begin()), data_itr_end( DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // get sequence separation
        const size_t sequence_separation( data_itr->SequenceSeparation());

        // true if the sequence separation is defined
        if( util::IsDefined( sequence_separation))
        {
          // take the log of the sequence separation and add it to the running sum
          seq_sep_sum += std::log( double( sequence_separation));

          // increment the count of pairs being considered
          ++count;
        }
      }

      // calculate the score
      const double score( -seq_sep_sum / ( double( count) * std::log( double( m_SequenceLength))));

      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DataSetPairwiseSequenceSeparation::Read( std::istream &ISTREAM)
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
    std::ostream &DataSetPairwiseSequenceSeparation::Write( std::ostream &OSTREAM, const size_t INDENT) const
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

    //! @brief determines the sequence separation between two chain and sequence ids
    //! @param CHAIN_ID_A first chain id
    //! @param SEQ_ID_A first seq id
    //! @param CHAIN_ID_B second chain id
    //! @param SEQ_ID_B second seq id
    //! @return size_t which is the sequence separation between the first and second
    size_t DataSetPairwiseSequenceSeparation::CalculateSequenceSeparation
    (
      const char CHAIN_ID_A, const int SEQ_ID_A, const char CHAIN_ID_B, const int SEQ_ID_B
    )
    {
      // true if the chain ids don't match
      if( CHAIN_ID_A != CHAIN_ID_B)
      {
        return util::GetUndefinedSize_t();
      }

      // return difference between SEQ_ID_A and SEQ_ID_B
      return math::Absolute( SEQ_ID_A - SEQ_ID_B);
    }

  } // namespace score
  
} // namespace bcl
