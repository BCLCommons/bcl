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
#include "restraint/bcl_restraint_mutate_data_set_pairwise_filter_sse_size.h"

// includes from bcl - sorted alphabetically
#include "restraint/bcl_restraint_data_set_pairwise.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateDataSetPairwiseFilterSSESize::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateDataSetPairwiseFilterSSESize())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &MutateDataSetPairwiseFilterSSESize::GetDefaultScheme()
    {
      static const std::string s_scheme( "filter_sse_size");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor taking optional scheme
    //! @param SCHEME the scheme for this scoring function
    MutateDataSetPairwiseFilterSSESize::MutateDataSetPairwiseFilterSSESize( const std::string &SCHEME) :
      m_SSEPool(),
      m_MinSSELengths(),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking optional scheme
    //! @param SSE_POOL pool to use as sse definitions
    //! @param MIN_SSE_LENGTHS minimum lengths each sse size should have for a data point to be kept if within an sse
    //! @param SCHEME the scheme for this scoring function
    MutateDataSetPairwiseFilterSSESize::MutateDataSetPairwiseFilterSSESize
    (
      const util::ShPtr< assemble::SSEPool> &SSE_POOL,
      const storage::Map< biol::SSType, size_t> &MIN_SSE_LENGTHS,
      const std::string &SCHEME
    ) :
      m_SSEPool( SSE_POOL),
      m_MinSSELengths( MIN_SSE_LENGTHS),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateDataSetPairwiseFilterSSESize
    MutateDataSetPairwiseFilterSSESize *MutateDataSetPairwiseFilterSSESize::Clone() const
    {
      return new MutateDataSetPairwiseFilterSSESize( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateDataSetPairwiseFilterSSESize::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &MutateDataSetPairwiseFilterSSESize::GetScheme() const
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
    MutateDataSetPairwiseFilterSSESize::operator()( const DataSetPairwise &DATA_SET) const
    {
      static util::ShPtr< DataSetPairwise> s_empty;

      if( DATA_SET.IsEmpty())
      {
        BCL_MessageStd( "MutateDataSetPairwiseFilterSSESize  data set is empty");
        return math::MutateResult< DataSetPairwise>( s_empty, *this);
      }

      // get random non overlapping set of sses from the pool
      const storage::Set< util::SiPtr< const assemble::SSE>, assemble::SSELessThanNoOverlap> sse_set
      (
        m_SSEPool->GetRandomNonOverlappingSet()
      );

      // the dataset filtered
      util::ShPtr< DataSetPairwise> filtered_data( new DataSetPairwise());

      // iterate through the data set to filter out data pairs that don't meet sse length requirement
      for
      (
        DataSetPairwise::const_iterator data_itr( DATA_SET.Begin()), data_itr_end( DATA_SET.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // locate the sses that each data point falls in
        const util::SiPtr< const assemble::SSE> sse_a( data_itr->First()->LocateSSE( sse_set));
        const util::SiPtr< const assemble::SSE> sse_b( data_itr->Second()->LocateSSE( sse_set));

        // true if either of the sses are not located
        if( !sse_a.IsDefined() || !sse_b.IsDefined())
        {
          // go to next data pair
          continue;
        }

        // find the sse types in m_MinSSELengths
        storage::Map< biol::SSType, size_t>::const_iterator
          ss_type_itr_a( m_MinSSELengths.Find( sse_a->GetType())), ss_type_itr_b( m_MinSSELengths.Find( sse_b->GetType()));

        // true if either of the sses types are not in the min sse size map
        if( ss_type_itr_a == m_MinSSELengths.End() || ss_type_itr_b == m_MinSSELengths.End())
        {
          // go to next data pair
          continue;
        }

        // check whether or not the sses meet the size requirement
        const bool sse_a_meet_size_requirement( sse_a->GetSize() >= ss_type_itr_a->second ? true : false);
        const bool sse_b_meet_size_requirement( sse_b->GetSize() >= ss_type_itr_b->second ? true : false);

        // true if both the sses meet the size requirement
        if( sse_a_meet_size_requirement && sse_b_meet_size_requirement)
        {
          // add the current data pair to the filtered dataset
          filtered_data->Insert( *data_itr);
        }
      }

      BCL_MessageStd
      (
        "MutateDataSetPairwiseFilterSSESize num unique locators " +
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
    std::istream &MutateDataSetPairwiseFilterSSESize::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SSEPool, ISTREAM);
      io::Serialize::Read( m_MinSSELengths, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateDataSetPairwiseFilterSSESize::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SSEPool, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MinSSELengths, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint
  
} // namespace bcl
