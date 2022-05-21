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
#include "restraint/bcl_restraint_accessibility_profile.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom.h"
#include "restraint/bcl_restraint_accessibility_profile_assignment.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AccessibilityProfile::s_Instance
    (
      GetObjectInstances().AddInstance( new AccessibilityProfile())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AccessibilityProfile::AccessibilityProfile() :
      m_Accessibilities()
    {
    }

    //! @brief constructor taking member variable
    //! @param ACCESSIBILITIES the list of accessibilities indicating a series of accessibilit measurements
    AccessibilityProfile::AccessibilityProfile( const storage::List< AccessibilityAA> &ACCESSIBILITIES) :
      m_Accessibilities( ACCESSIBILITIES)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AccessibilityProfile
    AccessibilityProfile *AccessibilityProfile::Clone() const
    {
      return new AccessibilityProfile( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AccessibilityProfile::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function for creating an assignment of the profile to residues in a protein model
    //! @param PROTEIN_MODEL model from which the assignment will be created
    //! @return AccessibilityProfileAssignment which assigns the profile to the residues in PROTEIN_MODEL
    AccessibilityProfileAssignment
    AccessibilityProfile::GenerateAssignment( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // list to hold accessibility assignments which are not for residues in sses
      storage::List< AccessibilityAAAssignment> non_sse_assignments;

      // map for storing sses and the associated accessibility measurements
      storage::Map
      <
        util::SiPtr< const assemble::SSE>, storage::List< AccessibilityAAAssignment>, assemble::SSELessThanNoOverlap
      > sse_assignments;

      // iterate through the AccessibilityAAs to generate assignments
      for
      (
        storage::List< AccessibilityAA>::const_iterator
          itr( m_Accessibilities.Begin()), itr_end( m_Accessibilities.End()); itr != itr_end; ++itr
      )
      {
        // generate the current assignment
        const AccessibilityAAAssignment assignment( itr->GenerateAssignment( PROTEIN_MODEL));

        // true if the assignment is not valid
        if
        (
          !assignment.GetAABase().IsDefined() || !util::IsDefined( assignment.GetExposureValue())
        )
        {
          // add assignment to non sse assignments list
          non_sse_assignments.PushBack( assignment);

          // go to next accessibility restraint
          continue;
        }

        // get the center amino acid - the current amino acid of interest that has an accessibility measurement
        const biol::AABase &aa( *assignment.GetAABase());

        // locate the sse the amino acid is in
        util::SiPtr< const assemble::SSE> sse
        (
          assemble::LocatorAA( aa.GetChainID(), aa.GetSeqID()).LocateSSE( PROTEIN_MODEL)
        );

        // make sure at this point the sse should be found
        BCL_Assert( sse.IsDefined(), "could not find sse in protein model");

        // make sure coordinates are defined
        BCL_Assert
        (
          aa.GetFirstSidechainAtom().GetCoordinates().IsDefined(),
          "GetFirstSidechainAtom coords not defined for " + aa.GetIdentification() + "\n\n" +
          util::Format()( aa) + "\n\nfirst side chain atom is " + aa.GetFirstSidechainAtom().GetType().GetName()
        );

        // add the assignment to the list of assignments for the current sse
        sse_assignments[ sse].PushBack( assignment);
      }

      // return AccessibilityProfileAssignment
      return AccessibilityProfileAssignment( sse_assignments, non_sse_assignments);
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
    std::istream &AccessibilityProfile::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Accessibilities, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AccessibilityProfile::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Accessibilities, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint

} // namespace bcl
