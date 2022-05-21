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
#include "assemble/bcl_assemble_pick_sses_random.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> PickSSEsRandom::s_Instance
    (
      GetObjectInstances().AddInstance( new PickSSEsRandom())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! @detail sets the SSE types to be considered to helix and strand and the number of SSE to pick to one
    PickSSEsRandom::PickSSEsRandom() :
      m_SSTypes(),
      m_NumSSEsToPick( 1)
    {
      m_SSTypes.Insert( biol::GetSSTypes().HELIX);
      m_SSTypes.Insert( biol::GetSSTypes().STRAND);
    }

    //! @brief constructor from a set of SSTypes and number of SSEs to pick
    //! @param SS_TYPES types of the SSEs which shall be picked
    //! @param NUM_SSES number of SSEs to pick randomly
    PickSSEsRandom::PickSSEsRandom( const storage::Set< biol::SSType> &SS_TYPES, const size_t NUM_SSES) :
      m_SSTypes( SS_TYPES),
      m_NumSSEsToPick( NUM_SSES)
    {
    }

    //! @brief returns a pointer to a new PickSSEsRandom
    //! @return pointer to a new PickSSEsRandom
    PickSSEsRandom *PickSSEsRandom::Clone() const
    {
      return new PickSSEsRandom( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &PickSSEsRandom::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief randomly picks the given number of SSEs of the given type from the list
    //! @param SSE_LIST SSE list to pick the SSEs from
    //! @return list with the picked SSEs
    util::SiPtrList< const SSE> PickSSEsRandom::Pick( const util::SiPtrList< const SSE> &SSE_LIST) const
    {
      // collect the SSEs of the wanted type from the list
      util::SiPtrList< const SSE> sses;
      for
      (
        util::SiPtrList< const SSE>::const_iterator it( SSE_LIST.Begin()), it_end( SSE_LIST.End());
        it != it_end;
        ++it
      )
      {
        if( m_SSTypes.Contains( ( **it).GetType()))
        {
          sses.PushBack( *it);
        }
      }

      // if the resulting list is shorter than the number of SSEs to pick exit
      BCL_Assert
      (
        m_NumSSEsToPick <= sses.GetSize(),
        "tried to pick " + util::Format()( m_NumSSEsToPick) + " SSEs from a list that only contains " +
        util::Format()( sses.GetSize()) + " SSEs"
      );

      // randomly pick the given number SSEs from the list containing the given SSE types
      util::SiPtrList< const SSE> picked_sses;
      for( size_t num_picked_sses( 0); num_picked_sses < m_NumSSEsToPick; ++num_picked_sses)
      {
        util::SiPtrList< const SSE>::iterator random_iterator
        (
          random::GetGlobalRandom().Iterator( sses.Begin(), sses.End(), sses.GetSize())
        );
        picked_sses.PushBack( *random_iterator);
        sses.Remove( random_iterator);
      }

      return picked_sses;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read object from input stream
    //! @param ISTREAM input stream to read object from
    //! @return input stream which was read from
    std::istream &PickSSEsRandom::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SSTypes, ISTREAM);
      io::Serialize::Read( m_NumSSEsToPick, ISTREAM);

      return ISTREAM;
    }

    //! @brief write object into  output stream
    //! @param OSTREAM output stream to write object into
    //! @param INDENT number of indentations to separate members
    //! @return output stream object was written into
    std::ostream &PickSSEsRandom::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SSTypes, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NumSSEsToPick, OSTREAM, INDENT);

      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
