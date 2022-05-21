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
#include "score/bcl_score_data_set_pairwise_sse_size.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_running_average_sd.h"
#include "math/bcl_math_running_min_max.h"
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
    const util::SiPtr< const util::ObjectInterface> DataSetPairwiseSSESize::s_Instance
    (
      GetObjectInstances().AddInstance( new DataSetPairwiseSSESize())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &DataSetPairwiseSSESize::GetDefaultScheme()
    {
      static const std::string s_scheme( "sse_size");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor taking optional scheme
    //! @param SCHEME the scheme for this scoring function
    DataSetPairwiseSSESize::DataSetPairwiseSSESize( const std::string &SCHEME) :
      m_SSEPool(),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking optional scheme
    //! @param SSE_POOL pool to use as sse definitions
    //! @param SCHEME the scheme for this scoring function
    DataSetPairwiseSSESize::DataSetPairwiseSSESize
    (
      const util::ShPtr< assemble::SSEPool> &SSE_POOL, const std::string &SCHEME
    ) :
      m_SSEPool( SSE_POOL),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new DataSetPairwiseSSESize
    DataSetPairwiseSSESize *DataSetPairwiseSSESize::Clone() const
    {
      return new DataSetPairwiseSSESize( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DataSetPairwiseSSESize::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &DataSetPairwiseSSESize::GetScheme() const
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
    double DataSetPairwiseSSESize::operator()( const restraint::DataSetPairwise &DATA) const
    {
      // get random non overlapping set of sses from the pool
      const storage::Set< util::SiPtr< const assemble::SSE>, assemble::SSELessThanNoOverlap> sse_set
      (
        m_SSEPool->GetRandomNonOverlappingSet()
      );

      // get siptrvector from the set of sses
      util::SiPtrVector< const assemble::SSE> sses( sse_set.Begin(), sse_set.End());

      // keep track of the average size of sse with data point in it
      math::RunningAverageSD< double> average_size;

      // keep track of the max size of any sse
      math::RunningMinMax< double> max_size;

      // iterate through the sses to get the max size
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator sse_itr( sses.Begin()), sse_itr_end( sses.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // consider the size of this sse to see if it is largest
        max_size += double( ( *sse_itr)->GetSize());

        // iterate through the data set to see if any data points are in the current sse
        for
        (
          restraint::DataSetPairwise::const_iterator data_itr( DATA.Begin()), data_itr_end( DATA.End());
          data_itr != data_itr_end;
          ++data_itr
        )
        {
          // true if the first data point is in the sse
          if( data_itr->First()->IsWithin( **sse_itr))
          {
            average_size += ( *sse_itr)->GetSize();
          }
          // true if the second data point is in the sse
          if( data_itr->Second()->IsWithin( **sse_itr))
          {
            average_size += ( *sse_itr)->GetSize();
          }
        }
      }

      // score is fraction of average compared to the maximum size
      const double score( -average_size.GetAverage() / max_size.GetMax());

      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DataSetPairwiseSSESize::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SSEPool, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DataSetPairwiseSSESize::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SSEPool, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace score

} // namespace bcl
