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
#include "restraint/bcl_restraint_accessibility_profile_assignment.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AccessibilityProfileAssignment::s_Instance
    (
      GetObjectInstances().AddInstance( new AccessibilityProfileAssignment())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AccessibilityProfileAssignment::AccessibilityProfileAssignment() :
      m_SSEAssignments(),
      m_NonSSEAssignments()
    {
    }

    //! @brief constructor from member variables
    //! @param SSE_ASSIGNMENTS map connecting SSEs to the individual accessibility profiles which is stored as list
    //! @param NON_SSE_ASSIGNMENTS holds any AccessibilityAAAssignments that are not associated with an SSE
    AccessibilityProfileAssignment::AccessibilityProfileAssignment
    (
      const storage::Map
      <
        util::SiPtr< const assemble::SSE>, storage::List< AccessibilityAAAssignment>, assemble::SSELessThanNoOverlap
      > &SSE_ASSIGNMENTS,
      const storage::List< AccessibilityAAAssignment> &NON_SSE_ASSIGNMENTS
    ) :
      m_SSEAssignments( SSE_ASSIGNMENTS),
      m_NonSSEAssignments( NON_SSE_ASSIGNMENTS)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AccessibilityProfileAssignment
    AccessibilityProfileAssignment *AccessibilityProfileAssignment::Clone() const
    {
      return new AccessibilityProfileAssignment( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AccessibilityProfileAssignment::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief data access to map connecting SSEs to the individual accessibility profiles which is stored as list
    //! @return map which is m_SSEAssignments
    const storage::Map
    <
      util::SiPtr< const assemble::SSE>, storage::List< AccessibilityAAAssignment>, assemble::SSELessThanNoOverlap
    >
    &AccessibilityProfileAssignment::GetSSEAssignments() const
    {
      return m_SSEAssignments;
    }

    //! @brief gives the total number of assignments within SSES
    //! @return size_t which is the number of assignments within sses
    size_t AccessibilityProfileAssignment::GetTotalNumberOfSSEAssignments() const
    {
      size_t num_assigns( 0);

      // iterate over the sse assignment map
      for
      (
        storage::Map
        <
          util::SiPtr< const assemble::SSE>, storage::List< AccessibilityAAAssignment>, assemble::SSELessThanNoOverlap
        >::const_iterator itr( m_SSEAssignments.Begin()), itr_end( m_SSEAssignments.End());
        itr != itr_end;
        ++itr
      )
      {
        // add the current number of assignments in the current sse
        num_assigns += itr->second.GetSize();
      }

      return num_assigns;
    }

    //! @brief data access to any AccessibilityAAAssignments that are not associated with an SSE
    //! @return list which is m_NonSSEAssignments
    const storage::List< AccessibilityAAAssignment> &AccessibilityProfileAssignment::GetNonSSEAssignments() const
    {
      return m_NonSSEAssignments;
    }

    //! @brief counts the number of residues with accessibility data of given environment type in all sses
    //! @param ENVIRONMENT the environment that will be counted
    //! @return size_t count of the number of residues with accessibility data of given environment type in all sses
    storage::Map< util::SiPtr< const biol::AABase>, double> AccessibilityProfileAssignment::GetSSEAccessibilities
    (
      const assemble::SSE &SSE, const AccessibilityAA::EnvironmentType &ENVIRONMENT
    ) const
    {
      // try to find the sse in the map
      storage::Map< util::SiPtr< const assemble::SSE>, storage::List< AccessibilityAAAssignment> >::const_iterator
        itr_sse( m_SSEAssignments.Find( util::ToSiPtr( SSE)));

      // true if the sse could not be found
      if( itr_sse == m_SSEAssignments.End())
      {
        return storage::Map< util::SiPtr< const biol::AABase>, double>();
      }

      // will hold the residues associated with the given sse and their accessibilty for the desired environment
      storage::Map< util::SiPtr< const biol::AABase>, double> data;

      // iterate through list of accessibility assignments associated with this sse
      for
      (
        storage::List< AccessibilityAAAssignment>::const_iterator
          data_itr( itr_sse->second.Begin()), data_itr_end( itr_sse->second.End());
        data_itr != data_itr_end; ++data_itr
      )
      {
        // get desired accessibility type by environment
        storage::Pair< bool, double> environ_type( data_itr->GetAccessibilityByEnvironment( ENVIRONMENT));

        // data for environment doesn't exist
        if( !environ_type.First())
        {
          continue;
        }

        // add the residue and its accessibility to the data map
        data.Insert
        (
          std::pair< util::SiPtr< const biol::AABase>, double>
          (
            data_itr->GetAABase(), environ_type.Second()
          )
        );
      }

      return data;
    }

    //! @brief counts the number of residues with accessibility data of given environment type in all sses
    //! @param ENVIRONMENT the environment that will be counted
    //! @return size_t count of the number of residues with accessibility data of given environment type in all sses
    size_t AccessibilityProfileAssignment::GetNumberResiduesInSSEsWithEnvironmentType
    (
      const AccessibilityAA::EnvironmentType &ENVIRONMENT
    ) const
    {
      // will count the number of residues in sses with the given ENVIRONMENT
      size_t count( 0);

      // iterate through the data of sses and their associated accessibility assignments
      for
      (
        storage::Map
        <
          util::SiPtr< const assemble::SSE>, storage::List< AccessibilityAAAssignment>, assemble::SSELessThanNoOverlap
        >::const_iterator
          sse_itr( m_SSEAssignments.Begin()), sse_itr_end( m_SSEAssignments.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // get the list of accessibilities associated with the current sse and add the size to count
        count += GetSSEAccessibilities( *sse_itr->first, ENVIRONMENT).GetSize();
      }

      BCL_MessageDbg( "count is " + util::Format()( count));
      return count;
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AccessibilityProfileAssignment::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SSEAssignments, ISTREAM);
      io::Serialize::Read( m_NonSSEAssignments, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AccessibilityProfileAssignment::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SSEAssignments, OSTREAM, INDENT);
      io::Serialize::Write( m_NonSSEAssignments, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint

} // namespace bcl
