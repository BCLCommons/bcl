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
#include "assemble/bcl_assemble_sse_pool_insert_coil_into_sse.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_pool.h"
#include "math/bcl_math_mutate_result.h"
#include "math/bcl_math_running_average_sd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSEPoolInsertCoilIntoSSE::s_Instance
    (
      GetObjectInstances().AddInstance
      (
        new SSEPoolInsertCoilIntoSSE
        (
          sspred::GetMethods().e_Undefined,
          util::ShPtr< find::LocatorInterface< util::SiPtr< const SSE>, DomainInterface> >(),
          1
        )
      )
    );

    //! @brief the default scheme for this class
    const std::string &SSEPoolInsertCoilIntoSSE::GetDefaultScheme()
    {
      static const std::string s_default_scheme( "pool_insert_coil_into_sse");
      return s_default_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from locator and mutate
    //! @param SS_METHOD the method to use to locate the largest drop in the ss prediction for the located sse
    //! @param SP_LOCATE_SSE picker of sse from domain
    //! @param COIL_LENGTH length of the coil to insert
    //! @param SCHEME the scheme
    SSEPoolInsertCoilIntoSSE::SSEPoolInsertCoilIntoSSE
    (
      const sspred::Method &SS_METHOD,
      const util::ShPtr< find::LocatorInterface< util::SiPtr< const SSE>, DomainInterface> > &SP_LOCATE_SSE,
      const size_t COIL_LENGTH,
      const std::string &SCHEME
    ) :
      m_Method( SS_METHOD),
      m_SSELocator( SP_LOCATE_SSE),
      m_CoilLength( COIL_LENGTH),
      m_Scheme( SCHEME)
    {
      BCL_Assert( m_CoilLength > 0, "coil length needs to be at least 1, but is: " + util::Format()( m_CoilLength));
    }

    //! @brief Clone function
    //! @return pointer to new SSEPoolInsertCoilIntoSSE
    SSEPoolInsertCoilIntoSSE *SSEPoolInsertCoilIntoSSE::Clone() const
    {
      return new SSEPoolInsertCoilIntoSSE( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEPoolInsertCoilIntoSSE::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this mutate
    //! @return the scheme for this mutate
    const std::string &SSEPoolInsertCoilIntoSSE::GetScheme() const
    {
      return m_Scheme;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that mutates a pool by splitting a single sse
    //! @param SSE_POOL
    //! @return MutateResult that results from mutating to the SSE_POOL
    math::MutateResult< SSEPool> SSEPoolInsertCoilIntoSSE::operator()( const SSEPool &SSE_POOL) const
    {
      // pick a sse from the pool
      util::SiPtrList< const SSE> pool_sses( SSE_POOL.Begin(), SSE_POOL.End());
      const util::SiPtr< const SSE> picked_sse( m_SSELocator->Locate( SSE_POOL));

      // was picking successful
      if( !picked_sse.IsDefined())
      {
        return math::MutateResult< SSEPool>( util::ShPtr< SSEPool>(), *this);
      }

      const biol::SSType picked_type( picked_sse->GetType());

      // length of sse to split
      const size_t seq_length( picked_sse->GetSize());

      // need at least coil residues plus two to split
      if( seq_length < m_CoilLength)
      {
        return math::MutateResult< SSEPool>( util::ShPtr< SSEPool>(), *this);
      }

      // determine all n residue averages
      storage::List< storage::VectorND< 3, math::RunningAverageSD< double> > > mean_list;
      if( m_CoilLength && picked_sse->GetSize() <= m_CoilLength)
      {
        // iterate over the sequence
        for
        (
          biol::AASequence::const_iterator aa_itr( picked_sse->Begin()), aa_itr_end( picked_sse->End() - m_CoilLength + 1);
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {
          storage::VectorND< 3, math::RunningAverageSD< double> > stats;

          for
          (
            biol::AASequence::const_iterator aa_itr_local( aa_itr), aa_itr_local_end( aa_itr + m_CoilLength);
            aa_itr_local != aa_itr_local_end;
            ++aa_itr_local
          )
          {
            const linal::Vector3D current_prediction( ( *aa_itr_local)->GetData()->GetSSPrediction( m_Method)->GetThreeStatePrediction());
            // add pseudo probability
            stats( 0) +=  current_prediction( 0);
            stats( 1) +=  current_prediction( 1);
            stats( 2) +=  current_prediction( 2);
          }
          mean_list.PushBack( stats);
        }
      }

      // find the lowest mean starting from pos 1
      storage::List< storage::VectorND< 3, math::RunningAverageSD< double> > >::const_iterator itr( mean_list.Begin()), itr_end( mean_list.End());
      math::RunningAverageSD< double> lowest( itr->operator ()( picked_type));
      size_t split_pos( 0);
      for( size_t i( 0); itr != itr_end; ++itr, ++i)
      {
        if( lowest.GetAverage() > itr->operator ()( picked_type).GetAverage())
        {
          lowest = itr->operator ()( picked_type);
          split_pos = i;
        }
      }

      util::ShPtr< SSE> sp_sse_left( new SSE( picked_sse->SubSequence( 0, split_pos), picked_type));
      util::ShPtr< SSE> sp_sse_coil( new SSE( picked_sse->SubSequence( split_pos, m_CoilLength), biol::GetSSTypes().COIL));
      util::ShPtr< SSE> sp_sse_right( new SSE( picked_sse->SubSequence( split_pos + m_CoilLength, seq_length - split_pos - m_CoilLength), picked_type));

      BCL_MessageDbg
      (
        "inserted coil into sse: " + picked_sse->GetIdentification() + "\nresulting in\n" +
        sp_sse_left->GetIdentification() + '\n' + sp_sse_coil->GetIdentification() + '\n' + sp_sse_right->GetIdentification()
      );

      // replace with new sses
      util::SiPtrList< const SSE>::iterator new_end
      (
        std::remove_if
        (
          pool_sses.Begin(), pool_sses.End(),
          SSECompare( *picked_sse)
        )
      );

      util::ShPtr< SSEPool> sp_pool( new SSEPool( util::SiPtrList< const SSE>( pool_sses.Begin(), new_end), false));
      if( sp_sse_left->GetSize() > 0)
      {
        sp_pool->Insert( sp_sse_left);
      }
      sp_pool->Insert( sp_sse_coil);
      if( sp_sse_right->GetSize() > 0)
      {
        sp_pool->Insert( sp_sse_right);
      }

      // end
      return math::MutateResult< SSEPool>( sp_pool, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEPoolInsertCoilIntoSSE::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Method    , ISTREAM);
      io::Serialize::Read( m_SSELocator, ISTREAM);
      io::Serialize::Read( m_CoilLength, ISTREAM);
      io::Serialize::Read( m_Scheme    , ISTREAM);

      BCL_Assert( m_CoilLength > 0, "coil length needs to be at least 1, but is: " + util::Format()( m_CoilLength));

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSEPoolInsertCoilIntoSSE::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Method    , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SSELocator, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_CoilLength, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme    , OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
