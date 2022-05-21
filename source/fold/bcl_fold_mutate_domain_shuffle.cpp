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
#include "fold/bcl_fold_mutate_domain_shuffle.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MutateDomainShuffle::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateDomainShuffle())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a max number swaps and scheme
    //! @param MAX_NUMBER_SWAPS Maximum number of swaps within one mutate
    //! @param BEND whether to bend the SSE after swapping to match the original phi/psi angles
    //! @param SCHEME Scheme to be used
    MutateDomainShuffle::MutateDomainShuffle
    (
      const size_t MAX_NUMBER_SWAPS,
      const bool BEND,
      const std::string &SCHEME
    ) :
      m_MaxNumberSwaps( MAX_NUMBER_SWAPS),
      m_Bend( BEND),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateDomainShuffle
    MutateDomainShuffle *MutateDomainShuffle::Clone() const
    {
      return new MutateDomainShuffle( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateDomainShuffle::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MutateDomainShuffle::GetAlias() const
    {
      static const std::string s_name( "MutateDomainShuffle");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutateDomainShuffle::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Shuffles locations and orientations of SSEs in domains.");
      serializer.AddInitializer
      (
        "maximum steps",
        "maximum number of steps",
        io::Serialization::GetAgent( &m_MaxNumberSwaps)
      );
      serializer.AddInitializer
      (
        "bend",
        "whether to adjust swapped SSEs to the dihedral angles of the previous SSE",
        io::Serialization::GetAgent( &m_Bend)
      );

      return serializer;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that takes a Domain and return a mutated Domain
    //! @param THIS_DOMAIN Domain which will be mutated
    //! @return MutateResult with the mutated Domain
    math::MutateResult< assemble::Domain> MutateDomainShuffle::operator()( const assemble::Domain &THIS_DOMAIN) const
    {
      // static empty domain
      static util::ShPtr< assemble::Domain> s_empty_domain;

      // get the sses from the Domain
      util::SiPtrVector< const assemble::SSE> sses( THIS_DOMAIN.GetSSEs());

      // if there are less than 2 SSEs
      if( sses.GetSize() < 2)
      {
        // warn user and return failed MutateResult
        BCL_MessageVrb
        (
          "The given domain has only " + util::Format()( sses.GetSize()) + " SSEs therefore skipping this move"
        );
        return math::MutateResult< assemble::Domain>( s_empty_domain, *this);
      }

      // determine total number of swaps
      const size_t number_swaps
      (
        random::GetGlobalRandom().Random< size_t>( 1, std::min( m_MaxNumberSwaps, THIS_DOMAIN.GetNumberSSEs() / 2))
      );

      BCL_MessageVrb( "Number of swaps to apply " + util::Format()( number_swaps));

      // new sses
      util::ShPtr< assemble::SSE> new_sses;

      // construct a new Domain
      util::ShPtr< assemble::Domain> new_domain( THIS_DOMAIN.Clone());

      // iterate over swaps
      for( size_t nr_swap( 0); nr_swap < number_swaps; ++nr_swap)
      {
        // pick a strand
        util::SiPtrVector< const assemble::SSE>::iterator itr_a
        (
          random::GetGlobalRandom().Iterator( sses.Begin(), sses.End(), sses.GetSize())
        );

        // make a copy of the SSE
        util::ShPtr< assemble::SSE> sse_a( ( *itr_a)->Clone());

        // remove the first SSE
        sses.RemoveElement( itr_a);

        // pick a second strand
        util::SiPtrVector< const assemble::SSE>::iterator itr_b
        (
          random::GetGlobalRandom().Iterator( sses.Begin(), sses.End(), sses.GetSize())
        );

        // make a copy of the SSE
        util::ShPtr< assemble::SSE> sse_b( ( *itr_b)->Clone());

        // remove the second SSE
        sses.RemoveElement( itr_b);

        BCL_MessageVrb
        (
          "Swap #" + util::Format()( nr_swap) + " " + sse_a->GetIdentification() + " with " + sse_b->GetIdentification()
        );

        // if the SSEs need to be bent and the SSEs are the same type
        if( m_Bend && sse_a->GetType() == sse_b->GetType())
        {
          // make a copy of the first SSE
          const util::ShPtr< assemble::SSE> sse_a_copy( sse_a->Clone());

          // fit the sses
          sse_a->FitToSSE( *sse_b);
          sse_b->FitToSSE( *sse_a_copy);
        }
        // don't bend the SSEs, keep their current phi/psi angles
        else
        {
          // prepare transformation from a to b and b to a
          math::TransformationMatrix3D transform_ab( math::Inverse( sse_a->GetOrientation()));
          transform_ab( sse_b->GetOrientation());
          math::TransformationMatrix3D transform_ba( math::Inverse( sse_b->GetOrientation()));
          transform_ba( sse_a->GetOrientation());

          // apply the transformations
          sse_a->Transform( transform_ab);
          sse_b->Transform( transform_ba);
        }

        // do the replace within the domain
        new_domain->Replace( sse_a);
        new_domain->Replace( sse_b);
      }

      // end
      return math::MutateResult< assemble::Domain>( new_domain, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateDomainShuffle::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_MaxNumberSwaps, ISTREAM);
      io::Serialize::Read( m_Bend, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateDomainShuffle::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_MaxNumberSwaps, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Bend, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
