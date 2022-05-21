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
#include "fold/bcl_fold_mutate_domain_sse_split.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateDomainSSESplit::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateDomainSSESplit())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateDomainSSESplit::MutateDomainSSESplit() :
      m_Scheme(),
      m_StandardDeviation(),
      m_RandomNumberGenerator( random::GetGlobalRandom())
    {
    }

    //! @brief default constructor
    //! @param SCHEME the scheme for this mutate
    //! @param STANDARD_DEVIATION standard deviation of gaussian distribution that biases splitting towards the middle
    //! @param RNG the random number generator that should be used
    MutateDomainSSESplit::MutateDomainSSESplit
    (
      const std::string &SCHEME, const double STANDARD_DEVIATION, const random::DistributionInterface &RNG
    ) :
      m_Scheme( SCHEME),
      m_StandardDeviation( STANDARD_DEVIATION),
      m_RandomNumberGenerator( RNG)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateDomainSSESplit
    MutateDomainSSESplit *MutateDomainSSESplit::Clone() const
    {
      return new MutateDomainSSESplit( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateDomainSSESplit::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this mutate
    //! @return the scheme for this mutate
    const std::string &MutateDomainSSESplit::GetScheme() const
    {
      return m_Scheme;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that takes a domain
    //! @param DOMAIN domain which will be mutated
    //! @return MutateResult with the mutated domain
    math::MutateResult< assemble::Domain>
    MutateDomainSSESplit::operator()( const assemble::Domain &THIS_DOMAIN) const
    {
      // static empty domain
      static util::ShPtr< assemble::Domain> s_empty_domain;

      // get the sses from this domain
      const util::SiPtrVector< const assemble::SSE> domain_sses( THIS_DOMAIN.GetSSEs());

      if( domain_sses.IsEmpty())
      {
        return math::MutateResult< assemble::Domain>( s_empty_domain, *this);
      }

      // make copy of the domain
      util::ShPtr< assemble::Domain> new_domain( THIS_DOMAIN.Clone());

      // iterate through the sses in order to split them
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
          sse_itr( domain_sses.Begin()), sse_itr_end( domain_sses.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // the current sse that will be split
        const assemble::SSE &current_sse( **sse_itr);

        // split the current sse
        storage::VectorND< 2, util::ShPtr< assemble::SSE> > splitted_sses( SplitSSE( current_sse));

        // true if the sse could not be split
        if( !splitted_sses.First().IsDefined() || !splitted_sses.Second().IsDefined())
        {
          // go to next sse in the domain
          continue;
        }

        // remove the current sse from the domain
        BCL_Assert( new_domain->Remove( current_sse), "could not remove " + current_sse.GetIdentification());

        // insert the two new sses into the domain
        BCL_Assert
        (
          new_domain->Insert( splitted_sses.First()), "could not remove " + splitted_sses.First()->GetIdentification()
        );
        BCL_Assert
        (
          new_domain->Insert( splitted_sses.Second()),
          "could not remove " + splitted_sses.Second()->GetIdentification()
        );
      }

      return math::MutateResult< assemble::Domain>( new_domain, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateDomainSSESplit::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Scheme, ISTREAM);
      io::Serialize::Read( m_StandardDeviation, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateDomainSSESplit::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_StandardDeviation, OSTREAM, INDENT) << "\n";

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief splits an sse at a random position determined by the random number generator
    //! @param CURRENT_SSE the sse which will be split
    //! @return vectorND 2 with 2 shared pointers to sses which were created from the split sse
    storage::VectorND< 2, util::ShPtr< assemble::SSE> >
    MutateDomainSSESplit::SplitSSE( const assemble::SSE &CURRENT_SSE) const
    {
      static const size_t min_split_size( 1);
      // the length of the sse
      const size_t sse_length( CURRENT_SSE.GetSize());

      // true if the sse length is not sufficient to split
      if( sse_length <= min_split_size)
      {
        BCL_MessageDbg
        (
          "cannot split an sse with less than " + util::Format()( min_split_size) +
          " residues " + CURRENT_SSE.GetIdentification()
        );

        return storage::VectorND< 2, util::ShPtr< assemble::SSE> >();
      }

      // get the seq id to split at
      size_t split_point( GetSplitPoint( sse_length));

      // true if the split point is outside the sequence range of the sse
      if( split_point < min_split_size)
      {
        split_point = min_split_size;
      }
      else if( split_point > sse_length - min_split_size)
      {
        split_point = sse_length - min_split_size;
      }

      // get the subsequences
      const biol::AASequence sequence_a( CURRENT_SSE.SubSequence( 0, split_point));
      const biol::AASequence sequence_b( CURRENT_SSE.SubSequence( split_point, sse_length - split_point));

      BCL_MessageDbg( "first  split part sequence is " +  sequence_a.GetSequenceIdentification());
      BCL_MessageDbg( "second split part sequence is " +  sequence_b.GetSequenceIdentification());

      // make the two new sses
      const util::ShPtr< assemble::SSE> new_sse_a( new assemble::SSE( sequence_a, CURRENT_SSE.GetType()));
      const util::ShPtr< assemble::SSE> new_sse_b( new assemble::SSE( sequence_b, CURRENT_SSE.GetType()));

      // will be returned
      const storage::VectorND< 2, util::ShPtr< assemble::SSE> > sses( new_sse_a, new_sse_b);

      // return the sses
      return sses;
    }

    size_t MutateDomainSSESplit::GetSplitPoint( const size_t SSE_LENGTH) const
    {
      // the middle of the sequence of this sse
      const double random( m_RandomNumberGenerator.Random( double( 0), double( 1)));

      BCL_MessageDbg( "random number for splitting is " + util::Format()( random));
      BCL_MessageDbg( "SSE_LENGTH for splitting is " + util::Format()( SSE_LENGTH));

      double sum( 0);

      // u
      const double mid_point( double( SSE_LENGTH) / 2.0);
      BCL_MessageDbg( "mid_point is " + util::Format()( mid_point));

      // s
      const double range( double( SSE_LENGTH) / 2.0);
      BCL_MessageDbg( "range is " + util::Format()( range));

      size_t cut_point( util::GetUndefinedDouble());

      for( size_t pos( 1); pos < SSE_LENGTH; ++pos)
      {
        const double current_probability
        (
          ( 1.0 / ( 2.0 * range)) * ( 1.0 + cos( ( ( pos - mid_point) / range) * math::g_Pi))
        );
        sum += current_probability;
        BCL_MessageDbg( "current_probability is " + util::Format()( current_probability));
        BCL_MessageDbg( "sum is " + util::Format()( sum));

        if( sum > random)
        {
          cut_point = pos;
          BCL_MessageDbg( "cut_point is " + util::Format()( cut_point));
          break;
        }
      }

      BCL_Assert( util::IsDefined( cut_point), "cut point is undefined");

      return cut_point;
    }

  } // namespace fold
} // namespace bcl
