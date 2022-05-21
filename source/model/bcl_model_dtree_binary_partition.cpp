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
#include "model/bcl_model_dtree_binary_partition.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> DtreeBinaryPartition::s_Instance
    (
      GetObjectInstances().AddInstance( new DtreeBinaryPartition())
    );

    //! @brief Data as string
    //! @param DATA the data of interest
    //! @return the string for the data
    const std::string &DtreeBinaryPartition::GetDataName( const DtreeBinaryPartition::Data &DATA)
    {
      static const std::string s_names[ s_NumberData + 1] =
      {
        "SplitRating",
        "InitialNumIncorrect",
        "RatingTimesInitialNumIncorrect",
        "InitialIncorrectPlusFinalCorrect",
        GetStaticClassName< DtreeBinaryPartition::Data>()
      };
      return s_names[ size_t( DATA)];
    }

    //! @brief default constructor
    DtreeBinaryPartition::DtreeBinaryPartition() :
      m_FeatureIndex( util::GetUndefined< size_t>()),
      m_SplitValue( util::GetUndefined< float>()),
      m_SplitRating( -std::numeric_limits< float>::max()),
      m_Size( 0),
      m_InitialIncorrect( 0),
      m_FinalIncorrect( 0)
    {
    }

    //! @brief constructor from members
    DtreeBinaryPartition::DtreeBinaryPartition
    (
      const size_t &FEATURE_INDEX,
      const float  &SPLIT_VALUE,
      const float  &SPLIT_RATING,
      const size_t &SIZE,
      const size_t &NUM_INCORRECT
    ) :
      m_FeatureIndex( FEATURE_INDEX),
      m_SplitValue( SPLIT_VALUE),
      m_SplitRating( SPLIT_RATING),
      m_Size( SIZE),
      m_InitialIncorrect( NUM_INCORRECT),
      m_FinalIncorrect( 0)
    {
    }

    //! @brief Clone function
    //! @return pointer to new DtreeBinaryPartition
    DtreeBinaryPartition *DtreeBinaryPartition::Clone() const
    {
      return new DtreeBinaryPartition( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DtreeBinaryPartition::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the value specified by the enum
    //! @param DATA the data to retrieve from the partition
    //! @return the value specified by the enum
    float DtreeBinaryPartition::GetData( const DtreeBinaryPartition::Data &DATA) const
    {
      switch( DATA)
      {
        case e_SplitRating:
          return m_SplitRating;
        case e_InitialNumIncorrect:
          return m_InitialIncorrect;
        case e_RatingTimesInitialNumIncorrect:
          return m_InitialIncorrect;
        case e_InitialIncorrectPlusFinalCorrect:
          return 2 * m_InitialIncorrect - m_FinalIncorrect;
        default:
          BCL_Exit( "No value for data " + GetDataName( DATA), -1);
      }
      return 0.0;
    }

    //! @brief add incorrect estimate info to a binary partition
    //! @param STATE_COUNTS total state counts (from determine total state counts on the post-partition)
    void DtreeBinaryPartition::DetermineAccuracy( const storage::Vector< FeatureResultAndState> &DATA)
    {
      m_Size = DATA.GetSize();
      if( !util::IsDefined( m_FeatureIndex))
      {
        m_FinalIncorrect = m_InitialIncorrect = 0;
        return;
      }

      // initialize count of elements on the left and right
      linal::Vector< size_t> lhs_state_counts( DATA.FirstElement().GetResultState().GetSize(), size_t( 0));
      linal::Vector< size_t> rhs_state_counts( lhs_state_counts);

      // track the left and right hand sides
      size_t lhs_size( 0);

      //iterate over candidate set, which includes boundaries between possible split values,
      for
      (
        storage::Vector< FeatureResultAndState>::const_iterator itr( DATA.Begin()), itr_end( DATA.End());
        itr != itr_end;
        ++itr
      )
      {
        // determine whether this feature is on the lhs or rhs
        if( itr->GetFeature()( m_FeatureIndex) <= m_SplitValue)
        {
          lhs_state_counts += itr->GetResultState();
          ++lhs_size;
        }
        else
        {
          rhs_state_counts += itr->GetResultState();
        }
      }
      const size_t rhs_size( m_Size - lhs_size);

      m_FinalIncorrect
        = EstimateNumberIncorrect( lhs_size, lhs_state_counts) + EstimateNumberIncorrect( rhs_size, rhs_state_counts);
      m_InitialIncorrect = EstimateNumberIncorrect( m_Size, lhs_state_counts + rhs_state_counts);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DtreeBinaryPartition::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_FeatureIndex, ISTREAM);
      io::Serialize::Read( m_SplitValue, ISTREAM);
      io::Serialize::Read( m_SplitRating, ISTREAM);
      io::Serialize::Read( m_InitialIncorrect, ISTREAM);
      io::Serialize::Read( m_FinalIncorrect, ISTREAM);
      io::Serialize::Read( m_Size, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DtreeBinaryPartition::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_FeatureIndex, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SplitValue, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SplitRating, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_InitialIncorrect, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_FinalIncorrect, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Size, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief estimate # incorrect
    //! @param SIZE the total # of elements represented by the state counts vector
    //! @param STATE_COUNTS the count of features in each state
    size_t DtreeBinaryPartition::EstimateNumberIncorrect
    (
      const size_t &SIZE,
      const linal::Vector< size_t> &STATE_COUNTS
    )
    {
      size_t estimate( 0);
      // handle the common case where there is only 1 result.  In this case, the incorrect counts are exactly right
      if( STATE_COUNTS.GetSize() == size_t( 1))
      {
        estimate = std::min( STATE_COUNTS( 0), SIZE - STATE_COUNTS( 0));
      }
      else if( STATE_COUNTS.GetSize() > size_t( 0))
      {
        // estimate the # of incorrectly classified elements in each state
        // the upper bound is given by the sum of state counts in impure states
        // the lower bound is given by the maximum state count in an impure state
        size_t incorrect_count_upper_bound( 0);
        size_t incorrect_count_lower_bound( 0);
        for
        (
          const size_t *itr( STATE_COUNTS.Begin()), *itr_end( STATE_COUNTS.End());
          itr != itr_end;
          ++itr
        )
        {
          const size_t incorrect_this_state( std::min( *itr, SIZE - *itr));
          incorrect_count_upper_bound += incorrect_this_state;
          incorrect_count_lower_bound = std::max( incorrect_count_lower_bound, incorrect_this_state);
        }

        incorrect_count_upper_bound = std::min( incorrect_count_upper_bound, SIZE);
        estimate = ( incorrect_count_lower_bound + incorrect_count_upper_bound) / 2;
      }
      return estimate;
    }

  } // namespace model
} // namespace bcl
