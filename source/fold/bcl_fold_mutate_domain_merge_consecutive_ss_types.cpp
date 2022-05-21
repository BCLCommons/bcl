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
#include "assemble/bcl_assemble_domain.h"
#include "fold/bcl_fold_mutate_domain_merge_consecutive_ss_types.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateDomainMergeConsecutiveSSTypes::MutateDomainMergeConsecutiveSSTypes() :
      m_Scheme( GetStaticClassName< MutateDomainMergeConsecutiveSSTypes>()),
      m_SSType()
    {
    }

    //! @brief default constructor
    MutateDomainMergeConsecutiveSSTypes::MutateDomainMergeConsecutiveSSTypes( const biol::SSType &SSTYPE) :
      m_Scheme( GetStaticClassName< MutateDomainMergeConsecutiveSSTypes>()),
      m_SSType( SSTYPE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateDomainMergeConsecutiveSSTypes
    MutateDomainMergeConsecutiveSSTypes *MutateDomainMergeConsecutiveSSTypes::Clone() const
    {
      return new MutateDomainMergeConsecutiveSSTypes( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateDomainMergeConsecutiveSSTypes::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns scheme
    //! @return the scheme of the object
    const std::string &MutateDomainMergeConsecutiveSSTypes::GetScheme() const
    {
      return m_Scheme;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////
    math::MutateResult< bcl::assemble::Domain> MutateDomainMergeConsecutiveSSTypes::operator()
    (
      const assemble::Domain &DOMAINS
    ) const
    {
      util::ShPtr< assemble::Domain> new_domain( DOMAINS.Clone());
      util::SiPtrVector< const assemble::SSE> domain_sses( DOMAINS.GetSSEs());
      storage::Set< util::SiPtr< const assemble::SSE>, assemble::SSELessThanNoOverlap> sses( domain_sses.Begin(), domain_sses.End());

      // iterate over all SSEs
      for
      (
        storage::Set< util::SiPtr< const assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
        itr( sses.Begin()), next( sses.Begin()), itr_end( sses.End()); next != itr_end;
      )
      {
        ++next;
        if( next == itr_end)
        {
          break;
        }

        // check if the current SSEs type matches the next ones
        if( ( *itr)->GetType() == ( *next)->GetType())
        {
          util::ShPtr< assemble::SSE> new_sse( ( *itr)->Clone());

          // keep iterating and assembling the new sse while the type is still the same
          while( ( *next)->GetType() == ( *itr)->GetType() && ( *itr)->GetType() == m_SSType && next != itr_end)
          {
            new_sse->AppendSequence( ( **next), false);

            // remove the sse just appended to the new sse out of the old domain
            new_domain->Remove( **next);
            ++next;
          }

          // remove the current position out of the old domain
          new_domain->Remove( **itr);

          // insert the new sse into the domain to fill up the gap left by the merged ones
          new_domain->Insert( new_sse);
        }

        // reset the iterator to the current position right after the last merged SSE
        itr = next;
      }

      return math::MutateResult< assemble::Domain>( new_domain, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateDomainMergeConsecutiveSSTypes::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_SSType, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateDomainMergeConsecutiveSSTypes::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {

      io::Serialize::Write( m_SSType, OSTREAM) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM);
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace fold
} // namespace bcl
