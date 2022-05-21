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

#ifndef BCL_RESTRAINT_ACCESSIBILITY_PROFILE_ASSIGNMENT_H_
#define BCL_RESTRAINT_ACCESSIBILITY_PROFILE_ASSIGNMENT_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_restraint_accessibility_aa_assignment.h"
#include "assemble/bcl_assemble_sse.h"
#include "assemble/bcl_assemble_sse_compare.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_si_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AccessibilityProfileAssignment
    //! @brief represents a profile of accessibility data measured for over secondary structure elements
    //! @details For multiple sses, each sse can has a list of associated accessibility. Measurements that are not
    //!          associated with an SSE are also contained.
    //!
    //! @see @link example_restraint_accessibility_profile_assignment.cpp @endlink
    //! @author alexanns
    //! @date Apr 7, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AccessibilityProfileAssignment :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! map connecting SSEs to the individual accessibility profiles which is stored as list
      storage::Map
      <
        util::SiPtr< const assemble::SSE>, storage::List< AccessibilityAAAssignment>, assemble::SSELessThanNoOverlap
      > m_SSEAssignments;

      //! holds any AccessibilityAAAssignments that are not associated with an SSE
      storage::List< AccessibilityAAAssignment> m_NonSSEAssignments;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AccessibilityProfileAssignment();

      //! @brief constructor from member variables
      //! @param SSE_ASSIGNMENTS map connecting SSEs to the individual accessibility profiles which is stored as list
      //! @param NON_SSE_ASSIGNMENTS holds any AccessibilityAAAssignments that are not associated with an SSE
      AccessibilityProfileAssignment
      (
        const storage::Map
        <
          util::SiPtr< const assemble::SSE>, storage::List< AccessibilityAAAssignment>, assemble::SSELessThanNoOverlap
        > &SSE_ASSIGNMENTS,
        const storage::List< AccessibilityAAAssignment> &NON_SSE_ASSIGNMENTS
      );

      //! @brief Clone function
      //! @return pointer to new AccessibilityProfileAssignment
      AccessibilityProfileAssignment *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief data access to map connecting SSEs to the individual accessibility profiles which is stored as list
      //! @return map which is m_SSEAssignments
      const storage::Map
      <
        util::SiPtr< const assemble::SSE>, storage::List< AccessibilityAAAssignment>, assemble::SSELessThanNoOverlap
      >
      &GetSSEAssignments() const;

      //! @brief gives the total number of assignments within SSES
      //! @return size_t which is the number of assignments within sses
      size_t GetTotalNumberOfSSEAssignments() const;

      //! @brief data access to any AccessibilityAAAssignments that are not associated with an SSE
      //! @return list which is m_NonSSEAssignments
      const storage::List< AccessibilityAAAssignment> &GetNonSSEAssignments() const;

      //! @brief counts the number of residues with accessibility data of given environment type in all sses
      //! @param ENVIRONMENT the environment that will be counted
      //! @return size_t count of the number of residues with accessibility data of given environment type in all sses
      size_t GetNumberResiduesInSSEsWithEnvironmentType( const AccessibilityAA::EnvironmentType &ENVIRONMENT) const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief gives accessibilities of a specific environment type for an sse
      //! @param SSE from which accessibilities will be gotten for residues with data
      //! @param ENVIRONMENT the accessibility environment type of interest
      //! @return map of aabase to double which is the value of the experimental accessibility value
      storage::Map< util::SiPtr< const biol::AABase>, double> GetSSEAccessibilities
      (
        const assemble::SSE &SSE, const AccessibilityAA::EnvironmentType &ENVIRONMENT
      ) const;

    }; // class AccessibilityProfileAssignment

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_ACCESSIBILITY_PROFILE_ASSIGNMENT_H_
