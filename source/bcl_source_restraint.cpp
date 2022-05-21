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
#include "restraint/bcl_restraint_accessibility_aa_assignment.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"

// external includes - sorted alphabetically
#include <algorithm>

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AccessibilityAAAssignment::s_Instance
    (
      GetObjectInstances().AddInstance( new AccessibilityAAAssignment())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AccessibilityAAAssignment::AccessibilityAAAssignment() :
      m_AminoAcid(),
      m_ExposureValue(),
      m_Accessibility()
    {
    }

    //! @brief constructor from member variables
    //! @param AMINO_ACID the amino acid that this assignment is for
    //! @param EXPOSURE_VALUE the calculated exposure
    //! @param EXPOSURE_CALCULATOR the method to use for calculating exposure of residues from structure
    //! @param ACCESSIBILITY the experimentally measured accessibilities
    AccessibilityAAAssignment::AccessibilityAAAssignment
    (
      const util::SiPtr< const biol::AABase> AMINO_ACID,
      const double EXPOSURE_VALUE,
      const storage::Map< AccessibilityAA::EnvironmentEnum, double> &ACCESSIBILITY
    ) :
      m_AminoAcid( AMINO_ACID),
      m_ExposureValue( EXPOSURE_VALUE),
      m_Accessibility( ACCESSIBILITY)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AccessibilityAAAssignment
    AccessibilityAAAssignment *AccessibilityAAAssignment::Clone() const
    {
      return new AccessibilityAAAssignment( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AccessibilityAAAssignment::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the amino acid this assignment is for
    //! @return the amino acid this assignment is for
    const util::SiPtr< const biol::AABase> &AccessibilityAAAssignment::GetAABase() const
    {
      return m_AminoAcid;
    }

    //! @brief get the calculated exposure value of this residue
    //! @return the exposure value calculated for this residue
    const double AccessibilityAAAssignment::GetExposureValue() const
    {
      return m_ExposureValue;
    }

    //! @brief data access to the experimental accessibilities
    //! @return set to experimental accessibilities which is m_Accessibility
    const storage::Map< AccessibilityAA::EnvironmentEnum, double> &AccessibilityAAAssignment::GetAccessibility() const
    {
      return m_Accessibility;
    }

    //! @brief provides the desired accessibility based on a given environment type
    //! @return pair of bool and Accessibility where the bool indicates if that environment type exists and
    //!         the Accessibility is the experimentally measured value
    storage::Pair< bool, double>
    AccessibilityAAAssignment::GetAccessibilityByEnvironment( const AccessibilityAA::EnvironmentType &ENVIRONMENT) const
    {
      storage::Map< AccessibilityAA::EnvironmentEnum, double>::const_iterator environ_type_itr
      (
        m_Accessibility.Find( ENVIRONMENT)
      );
      // data for environment doesn't exist
      if( environ_type_itr == m_Accessibility.End())
      {
        return storage::Pair< bool, double>( false, util::GetUndefinedDouble());
      }

      return storage::Pair< bool, double>( true, environ_type_itr->second);
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AccessibilityAAAssignment::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_AminoAcid, ISTREAM);
      io::Serialize::Read( m_ExposureValue, ISTREAM);
      io::Serialize::Read( m_Accessibility, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AccessibilityAAAssignment::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_AminoAcid, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ExposureValue, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Accessibility, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint

} // namespace bcl
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
#include "restraint/bcl_restraint_accessibility_aa.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "restraint/bcl_restraint_accessibility_aa_assignment.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AccessibilityAA::s_Instance
    (
      GetObjectInstances().AddInstance( new AccessibilityAA())
    );

    //! @brief EnvironmentType as string
    //! @param TYPE the type
    //! @return the string for TYPE
    const std::string &AccessibilityAA::GetEnvironmentName( const EnvironmentType &TYPE)
    {
      static const std::string s_names[] =
      {
        "Oxygen",
        "NiEDDA",
        GetStaticClassName< EnvironmentType>()
      };
      return s_names[ TYPE];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default consructor
    AccessibilityAA::AccessibilityAA() :
      m_Accessibility(),
      m_AminoAcid(),
      m_ExposureCalculator()
    {
    }

    //! @brief construct from all the member variables
    //! @param ACCESSIBILITY double which is the accessibility
    //! @param AMINO_ACID the type of accessibility measurement made
    //! @param EXPOSURE_CALCULATOR the method to use for calculating exposure of residues from structure
    AccessibilityAA::AccessibilityAA
    (
      const storage::Map< EnvironmentEnum, double> &ACCESSIBILITY,
      const util::ShPtr< assemble::LocatorAA> &AMINO_ACID,
      const util::ShPtr< assemble::AAExposureInterface> &EXPOSURE_CALCULATOR
    ) :
      m_Accessibility( ACCESSIBILITY),
      m_AminoAcid( AMINO_ACID),
      m_ExposureCalculator( EXPOSURE_CALCULATOR)
    {
    }

    //! @brief virtual copy constructor
    AccessibilityAA *AccessibilityAA::Clone() const
    {
      return new AccessibilityAA( *this);
    }

    //! @brief virtual destructor
    AccessibilityAA::~AccessibilityAA()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AccessibilityAA::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief provides the desired accessibility based on a given environment type
    //! @return pair of bool and Accessibility where the bool indicates if that environment type exists and
    //!         the Accessibility is the experimentally measured value
    storage::Pair< bool, double>
    AccessibilityAA::GetAccessibilityByEnvironment( const AccessibilityAA::EnvironmentType &ENVIRONMENT) const
    {
      storage::Map< EnvironmentEnum, double>::const_iterator environ_type_itr
      (
        m_Accessibility.Find( ENVIRONMENT)
      );
      // data for environment doesn't exist
      if( environ_type_itr == m_Accessibility.End())
      {
        return storage::Pair< bool, double>( false, util::GetUndefinedDouble());
      }

      return storage::Pair< bool, double>( true, environ_type_itr->second);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief GenerateNeighborList generates an AANeighborList from a protein model for the desired residue
    //! @param PROTEIN_MODEL is used to create the assemble::AANeighborList
    //! @return returns assemble::AANeighborList for the located residue from protein model
    assemble::AANeighborList AccessibilityAA::GenerateNeighborList
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // create const SiPtr to const AABase "amino_acid" and initialize with the residue located by "m_AminoAcid"
      const util::SiPtr< const biol::AABase> amino_acid( m_AminoAcid->Locate( PROTEIN_MODEL));

      // if the residue is not found then return assignment with the empty SiPtr where the AAbase should be
      if( !amino_acid.IsDefined())
      {
        return assemble::AANeighborList();
      }

      // create AllAANeighborList "neighbor_list" which will be filled with the neighbors for "amino_acid"
      const assemble::AANeighborList neighbor_list
      (
        *amino_acid,
        PROTEIN_MODEL.GetAminoAcids(),
        m_ExposureCalculator->GetDistanceCutoff(),
        m_ExposureCalculator->GetMinimalSequenceSeparation(),
        true
      );

      return neighbor_list;
    }

    //! @brief GenerateAssignment generates an assignment for the located residue from a protein model
    //! @param PROTEIN_MODEL is used to create the AccessibilityAAAssignment
    //! @return returns AccessibilityAAAssignment for the located residue from protein model
    AccessibilityAAAssignment AccessibilityAA::GenerateAssignment
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // get neighbor list
      const assemble::AANeighborList neighbor_list( GenerateNeighborList( PROTEIN_MODEL));

      // true if the residue could not be found
      if
      (
        neighbor_list.GetCenterAminoAcid() ==
          util::SiPtr< const biol::AABase>( assemble::AANeighborList::GetDefaultCenterAA())
      )
      {
        // return an assignment with an undefined residue (default neighbor list aa) and exposure
        return AccessibilityAAAssignment
        (
          util::SiPtr< const biol::AABase>(),
          util::GetUndefinedDouble(),
          m_Accessibility
        );
      }

      // calculate the exposure of the residue of interest
      const double exposure( m_ExposureCalculator->operator()( neighbor_list));

      // return AccessibilityAAAssignment
      return AccessibilityAAAssignment
      (
        util::SiPtr< const biol::AABase>( neighbor_list.GetCenterAminoAcid()),
        exposure,
        m_Accessibility
      );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read distance from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AccessibilityAA::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Accessibility, ISTREAM);
      io::Serialize::Read( m_AminoAcid, ISTREAM);
      io::Serialize::Read( m_ExposureCalculator, ISTREAM);
      return ISTREAM;
    }

    //! @brief write distance to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &AccessibilityAA::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Accessibility, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AminoAcid, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ExposureCalculator, OSTREAM, INDENT);

      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_analyze_accessibility_change.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_count.h"
#include "biol/bcl_biol_aa_classes.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serialization_via_static_functions.h"
#include "math/bcl_math_running_average_sd.h"
#include "math/bcl_math_statistics.h"
#include "restraint/bcl_restraint_accessibility_aa_assignment.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AnalyzeAccessibilityChange::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new AnalyzeAccessibilityChange())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AnalyzeAccessibilityChange::AnalyzeAccessibilityChange() :
      m_OutFilePostFix( ".AnalyzeAccessibilityChange"),
      m_StartEnsemble(),
      m_ExposureMethod( util::CloneToShPtr( assemble::AANeighborCount())),
      m_ExperimentalExposures(),
      m_MeanMinCutoff(),
      m_MeanMaxCutoff(),
      m_ZScoreMinCutoff(),
      m_ZScoreMaxCutoff(),
      m_PymolOuputFilename(),
      m_EnsembleRepresentativeIndex(),
      m_EnsembleRepresentativeFromStartEnsemble(),
      m_GradientMin(),
      m_GradientMax(),
      m_DirectRelation(),
      m_SetNCRange(),
      m_NCRangeMin(),
      m_NCRangeMax()
    {
    }

    //! @brief Clone function
    //! @return pointer to new AnalyzeAccessibilityChange
    AnalyzeAccessibilityChange *AnalyzeAccessibilityChange::Clone() const
    {
      return new AnalyzeAccessibilityChange( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AnalyzeAccessibilityChange::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &AnalyzeAccessibilityChange::GetOutFilePostfix() const
    {
      return m_OutFilePostFix;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &AnalyzeAccessibilityChange::GetAlias() const
    {
      static const std::string s_Name( "AccessibilityChange");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string AnalyzeAccessibilityChange::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // to hold, over all profiles, the average of the average exposure difference calculated across all the residues
      // in the profile
      math::RunningAverageSD< double> average_mean_diff_calculated;
      // ranges for normalizing the calculated exposures in to the range of the experimental data
      const storage::VectorND< 2, math::Range< double> > mean_range_exp_range
      (
        GetCalculatedAndExperimentalRanges( ENSEMBLE, average_mean_diff_calculated)
      );
      const math::Range< double> &mean_range( mean_range_exp_range.First());
      const math::Range< double> &exp_range( mean_range_exp_range.Second());

      const double all_mean_change_mean( average_mean_diff_calculated.GetAverage()),
                   all_mean_change_stddev( average_mean_diff_calculated.GetStandardDeviation());

      BCL_MessageDbg( "all_mean_change_mean "   + util::Format()( all_mean_change_mean));
      BCL_MessageDbg( "all_mean_change_stddev " + util::Format()( all_mean_change_stddev));

      // to hold the analysis text
      std::string analysis;

      // open write to pymol output script
      io::OFStream write;
      io::File::MustOpenOFStream( write, m_PymolOuputFilename);

      // iterator to the structure for visualization in pymol, by default comes from the end ensemble
      assemble::ProteinEnsemble::const_iterator represent_itr( ENSEMBLE.Begin()), represent_itr_end( ENSEMBLE.End());

      // true if the representative structure should come from the start ensemble
      if( m_EnsembleRepresentativeFromStartEnsemble)
      {
        // set iterator to start ensemble
        represent_itr     = m_StartEnsemble.Begin();
        represent_itr_end = m_StartEnsemble.End();
      }

      // move the iterator to the correct structure
      storage::AdvanceIterator( represent_itr, represent_itr_end, m_EnsembleRepresentativeIndex);

      // get the name of the pdb the structure came from
      util::ShPtr< util::Wrapper< std::string> > pdb_filename
      (
        ( *represent_itr)->GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile)
      );

      // write intitial pymol commands
      {
        std::string pdb( "dummy.pdb");

        if( pdb_filename.IsDefined())
        {
          pdb = std::string( pdb_filename->GetData());
        }

        write <<  "load " << pdb << ",access_model\n";
        write <<  "cmd.show_as(\"cartoon\"   ,\"access_model\")\n";
        write <<  "color grey70, access_model\n";
        write << "alter access_model, q=" + util::Format()( 0.0) << '\n';
      }

      // string for selection of residues with experimental information
      std::string experimental_residues( "select exp_resi, ");

      // vectors for experimental and accessibility values for calculating the spearman correlation
      storage::Vector< double> exp_vals;
      storage::Vector< double> access_vals;

      // iterate through the data
      for
      (
        storage::List< AccessibilityProfile>::const_iterator profile_itr( m_ExperimentalExposures.Begin()),
        profile_itr_end( m_ExperimentalExposures.End());
        profile_itr != profile_itr_end;
        ++profile_itr
      )
      {
        // to hold, for a profile, the average exposure difference calculated across all the residues in the profile
        math::RunningAverageSD< double> average_profile_diff;

        // to hold the experimental data associated with this profile
        double exp_data;

        // iterate through the profile
        for
        (
          storage::List< AccessibilityAA>::const_iterator aa_itr( profile_itr->GetAccessibilities().Begin()),
          aa_itr_end( profile_itr->GetAccessibilities().End());
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {
          // reference to current data
          const AccessibilityAA &current_data( *aa_itr);

          // iterate over the start ensemble
          for
          (
            assemble::ProteinEnsemble::const_iterator
              start_ensemble_itr( m_StartEnsemble.Begin()), start_ensemble_itr_end( m_StartEnsemble.End());
            start_ensemble_itr != start_ensemble_itr_end;
            ++start_ensemble_itr
          )
          {
            // the accessibility assignment for the starting ensemble
            const AccessibilityAAAssignment assignment_start( current_data.GenerateAssignment( **start_ensemble_itr));

            // iterate over the end ensemble
            for
            (
              assemble::ProteinEnsemble::const_iterator
               end_ensemble_itr( ENSEMBLE.Begin()), end_ensemble_itr_end( ENSEMBLE.End());
              end_ensemble_itr != end_ensemble_itr_end;
              ++end_ensemble_itr
            )
            {
              // the accessibility assignment for the ending ensemble
              AccessibilityAAAssignment assignment_end( current_data.GenerateAssignment( **end_ensemble_itr));

              // right now the data is just stored as oxygen accessibility
              // this is reset with each iteration but all exp data is the same across the entire profile
              exp_data = assignment_end.GetAccessibilityByEnvironment
              (
                AccessibilityAA::e_Oxygen
              ).Second();

              // add the current exposure difference
              average_profile_diff += assignment_end.GetExposureValue() - assignment_start.GetExposureValue();
            }
          } // start ensemble
        } // profile aas

        // get the profile and locators for the beginning and ending residues
        const storage::List< AccessibilityAA> &profile( profile_itr->GetAccessibilities());
        // casts from LocatorInterface to LocatorAA
        util::ShPtr< assemble::LocatorAA> start_locator( profile.FirstElement().GetAA());
        util::ShPtr< assemble::LocatorAA> end_locator( profile.LastElement().GetAA());

        // statistical values related to the accessibility profile
        const double stddev( average_profile_diff.GetStandardDeviation());
        const double mean( average_profile_diff.GetAverage());
        double zscore( stddev == 0 ? util::GetUndefinedDouble() : mean / stddev);

        // true if the profile meets the user specified cutoffs
        if
        (
          math::Absolute( mean) >= m_MeanMinCutoff && math::Absolute( mean) <= m_MeanMaxCutoff &&
          (
            (
              util::IsDefined( zscore) && math::Absolute( zscore) >= m_ZScoreMinCutoff &&
              math::Absolute( zscore) <= m_ZScoreMaxCutoff
            )
            ||
            ( !util::IsDefined( zscore))
          )
        )
        {
          double overall_zscore( ( mean - all_mean_change_mean) / all_mean_change_stddev);
          if( !m_DirectRelation)
          {
            overall_zscore = ( -mean - all_mean_change_mean) / all_mean_change_stddev;
          }

          // analysis file
          analysis += "|" + start_locator->GetIdentification() + "| |";
          analysis += end_locator->GetIdentification() + "|\texp_data:\t";
          analysis += util::Format()( exp_data) + "\tmodel_mean_change:\t";
          analysis += util::Format()( mean) + "\tmodel_stdev:\t";
          analysis += util::Format()( stddev) + "\tzscore:\t";
          analysis += util::Format()( zscore) + "\tre-ranged_mean:\t";
          analysis += util::Format()( exp_range.Rescale( mean, mean_range));
          analysis += util::Format()( "\toverall_zscore_exp_related\t");
          analysis += util::Format()( overall_zscore) + "\n";

          // for spearman correlation
          access_vals.PushBack( mean);
          exp_vals.PushBack( exp_data);
          BCL_MessageDbg( "mean " + util::Format()( mean) + " expval " + util::Format()( exp_data));

          // pymol file
          // iterate over the residue locators for the current profile
          for
          (
            storage::List< AccessibilityAA>::const_iterator aa_itr( profile_itr->GetAccessibilities().Begin()),
            aa_itr_end( profile_itr->GetAccessibilities().End());
            aa_itr != aa_itr_end;
            ++aa_itr
          )
          {
            // need to have a "+" if not the first residue for the selection
            if( aa_itr != profile_itr->GetAccessibilities().Begin())
            {
              experimental_residues += "+";
            }

            // reference to data
            const AccessibilityAA &current_data( *aa_itr);

            // try to cast the locator to a LocatorAA
            util::ShPtr< assemble::LocatorAA> current_locator( current_data.GetAA());

            // message and continue if the locator is not defined
            if( !current_locator.IsDefined())
            {
              BCL_MessageCrt
              (
                "locator " + util::Format()( current_data.GetAA()) +
                " could not be cast to a LocatorAA"
              );
              continue;
            }

            // set b factor to calculated exposure
//            write << "alter access_model and chain " << current_locator->GetLocatorChain().GetChainID()
//                  << " and resi " + util::Format()( current_locator->GetAAID()) + ", b=" +
//                     util::Format()( exp_range.Rescale( mean, mean_range))
//                  << '\n';

            // show sticks
//            write << "cmd.show( \"sticks\", \"access_model and chain "
//                  << current_locator->GetLocatorChain().GetChainID()
//                  << " and resi " + util::Format()( current_locator->GetAAID()) << "\")\n";

            // set q value to experimental data
            write << "alter access_model and chain " << current_locator->GetLocatorChain().GetChainID()
                  << " and resi " + util::Format()( current_locator->GetAAID()) + " and name CA, q=" +
                     util::Format()( exp_data - overall_zscore)
                  << '\n';

            // add residue to experimental residues selection
            experimental_residues +=
            (
                std::string( "(access_model and chain ")
              + util::Format()( current_locator->GetLocatorChain().GetChainID())
              + std::string( " and resi ") + util::Format()( current_locator->GetAAID()) + " )"
            );
          } // pymol file : iterate over the residue locators for the current profile
        } // if residue meets criteria
      } // profiles

      // calculate spearman correlation
      const double correlation
      (
        math::Statistics::CorrelationSpearman( access_vals.Begin(), access_vals.End(), exp_vals.Begin(), exp_vals.End())
      );

      BCL_MessageStd( "CorrelationSpearman " + util::Format()( correlation));

      // write experimental residues selection
      write << experimental_residues << "\n";

      // write the ramp information
//      write <<  "ramp_new grad, access_model, [" << m_GradientMin << ", " << ( m_GradientMax + m_GradientMin) / 2.0
//            << ", "
//            << m_GradientMax << "],[red,green,blue]\n";

      // write the spectrum information
      write <<  "spectrum q, rainbow_rev, exp_resi, minimum=" << m_GradientMin << "," << "maximum=" << m_GradientMax
            << "\n";

      // write spectrum for experimental values
//      std::string experimental_gradient( m_DirectRelation ? "blue_white_red" : "red_white_blue");
//      write <<  "spectrum q, " << experimental_gradient << ", exp_resi and name CA, minimum=" << m_GradientMin << ","
//            << "maximum=" << m_GradientMax << "\n";

      // hide hydrogens and close pymol script file
      write << "cmd.hide(\"(all and hydro)\")\n";
      io::File::CloseClearFStream( write);

      // return the analysis text
      return analysis;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AnalyzeAccessibilityChange::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_OutFilePostFix,                          ISTREAM);
      io::Serialize::Read( m_StartEnsemble,                           ISTREAM);
      io::Serialize::Read( m_ExposureMethod,                          ISTREAM);
      io::Serialize::Read( m_ExperimentalExposures,                   ISTREAM);
      io::Serialize::Read( m_MeanMinCutoff,                           ISTREAM);
      io::Serialize::Read( m_MeanMaxCutoff,                           ISTREAM);
      io::Serialize::Read( m_ZScoreMinCutoff,                         ISTREAM);
      io::Serialize::Read( m_ZScoreMaxCutoff,                         ISTREAM);
      io::Serialize::Read( m_PymolOuputFilename,                      ISTREAM);
      io::Serialize::Read( m_EnsembleRepresentativeIndex,             ISTREAM);
      io::Serialize::Read( m_EnsembleRepresentativeFromStartEnsemble, ISTREAM);
      io::Serialize::Read( m_GradientMin,                             ISTREAM);
      io::Serialize::Read( m_GradientMax,                             ISTREAM);
      io::Serialize::Read( m_DirectRelation,                          ISTREAM);
      io::Serialize::Read( m_SetNCRange,                              ISTREAM);
      io::Serialize::Read( m_NCRangeMin,                              ISTREAM);
      io::Serialize::Read( m_NCRangeMax,                              ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AnalyzeAccessibilityChange::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_OutFilePostFix,                          OSTREAM, INDENT);
      io::Serialize::Write( m_StartEnsemble,                           OSTREAM, INDENT);
      io::Serialize::Write( m_ExposureMethod,                          OSTREAM, INDENT);
      io::Serialize::Write( m_ExperimentalExposures,                   OSTREAM, INDENT);
      io::Serialize::Write( m_MeanMinCutoff,                           OSTREAM, INDENT);
      io::Serialize::Write( m_MeanMaxCutoff,                           OSTREAM, INDENT);
      io::Serialize::Write( m_ZScoreMinCutoff,                         OSTREAM, INDENT);
      io::Serialize::Write( m_ZScoreMaxCutoff,                         OSTREAM, INDENT);
      io::Serialize::Write( m_PymolOuputFilename,                      OSTREAM, INDENT);
      io::Serialize::Write( m_EnsembleRepresentativeIndex,             OSTREAM, INDENT);
      io::Serialize::Write( m_EnsembleRepresentativeFromStartEnsemble, OSTREAM, INDENT);
      io::Serialize::Write( m_GradientMin,                             OSTREAM, INDENT);
      io::Serialize::Write( m_GradientMax,                             OSTREAM, INDENT);
      io::Serialize::Write( m_DirectRelation,                          OSTREAM, INDENT);
      io::Serialize::Write( m_SetNCRange,                              OSTREAM, INDENT);
      io::Serialize::Write( m_NCRangeMin,                              OSTREAM, INDENT);
      io::Serialize::Write( m_NCRangeMax,                              OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AnalyzeAccessibilityChange::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates accessibility changes as (EndEnsembleAccessibility) - (BeginEnsembleAccessibility). Takes two "
        "ensembles of models. Calculates the mean, stddev, zscore of exposure change between them"
        " for residues that have experimental data. The experimental data is provided by an input file. Cutoffs "
        "(min and max) can be given for a residue to be outputted as of interest. Outputs a text file with "
        "these numbers and the information about the residues involved. Also, outputs a pymol script file "
        "which colors the backbone of residues according to the experimental accessibility change and colors the "
        "side chains of residues according to the calculated accessibility change."
      );

      parameters.AddInitializer
      (
        "filename_postfix",
        "the postfix that will be appended to the output file in order to identify this analysis",
        io::Serialization::GetAgent( &m_OutFilePostFix),
        ".AnalyzeAccessibilityChange"
      );
      parameters.AddInitializer
      (
        "start_ensemble_filename",
        "the name of the file containing the list of pdbs for the starting ensemble",
        util::OwnPtr< io::SerializationInterface>
        (
          new io::SerializationViaStaticFunctions< assemble::ProteinEnsemble>
          (
            command::ParameterCheckFileExistence(),
            ( &EnsembleAsFilename),
            ( &EnsembleFromFilename),
            &m_StartEnsemble
          )
        ),
        "start_ensemble_pdbs.ls"
      );

      parameters.AddInitializer
      (
        "experimental_data_filename",
        "the name of the file containing the list of experimental data\nShould have the format\n"
        "'A' 33  'A' 38  2\n'A' 41 'A' 52  2\n'A' 53  'A' 81  1\n"
        "Where -2 = large accessibility decrease; -1 = decrease; 0 = no change; 1 = increase; 2 = large "
        "accessibility increase",
        util::OwnPtr< io::SerializationInterface>
        (
          new io::SerializationViaStaticFunctions< storage::List< AccessibilityProfile> >
          (
            command::ParameterCheckFileExistence(),
            ( &ExposureDataAsFilename),
            ( &ExposureDataFromFilename),
            &m_ExperimentalExposures
          )
        ),
        "exposure_data.cst"
      );

      parameters.AddInitializer
      (
        "mean_min_cutoff",
        "The absolute value (inclusive) of the minimum mean calculated accessibility change for a residue to be "
        "considered at all",
        io::Serialization::GetAgent( &m_MeanMinCutoff),
        "0"
      );

      parameters.AddInitializer
      (
        "mean_max_cutoff",
        "The absolute value (inclusive) of the maximum mean calculated accessibility change for a residue to be "
        "considered at all",
        io::Serialization::GetAgent( &m_MeanMaxCutoff),
        "0"
      );

      parameters.AddInitializer
      (
        "zscore_min_cutoff",
        "The absolute value (inclusive) of the minimum zscore calculated for an accessibility change for a "
        "residue to be considered at all",
        io::Serialization::GetAgent( &m_ZScoreMinCutoff),
        "0"
      );

      parameters.AddInitializer
      (
        "zscore_max_cutoff",
        "The absolute value (inclusive) of the maximum zscore calculated for an accessibility change for a "
        "residue to be considered at all",
        io::Serialization::GetAgent( &m_ZScoreMaxCutoff),
        "0"
      );

      parameters.AddInitializer
      (
        "pymol_output_filename",
        "The filename of the pymol script that will be outputted showing the accessibility changes",
        io::Serialization::GetAgent( &m_PymolOuputFilename),
        "accessibilities.pml"
      );

      parameters.AddInitializer
      (
        "ensemble_representative_index",
        "The index of the model in the desired ensemble list that should be used as the representative in the "
        "pymol script. First model = index 0",
        io::Serialization::GetAgent( &m_EnsembleRepresentativeIndex),
        "0"
      );

      parameters.AddInitializer
      (
        "ensemble_representative_from_start_ensemble",
        "The model that should be used as the representative in the pymol script should come from the start "
        "ensemble. 1 = true; 0=false. If false representative will be taken from the end ensemble",
        io::Serialization::GetAgent( &m_EnsembleRepresentativeFromStartEnsemble),
        "0"
      );

      parameters.AddInitializer
      (
        "gradient_min",
        "The minimum value for the color gradient in pymol",
        io::Serialization::GetAgent( &m_GradientMin),
        "0"
      );

      parameters.AddInitializer
      (
        "gradient_max",
        "The maximum value for the color gradient in pymol",
        io::Serialization::GetAgent( &m_GradientMax),
        "0"
      );

      parameters.AddInitializer
      (
        "direct_relation",
        "Indicates the experimental data and exposures calculated are directly related. If true, a larger "
        "experimental value indicates the neighbor count should be larger (experiment measures buriedness). If"
        "false, a larger experimental value indicates the neigbor count should be smaller (experiment measures"
        "exposure).1=true;0=false",
        io::Serialization::GetAgent( &m_DirectRelation),
        "0"
      );

      parameters.AddInitializer
      (
        "set_nc_range",
        "Indicates that the range calculated from the NC's should not be determined automatically. Automatically"
        "the range is calculated as between the negative and postive value of the most extremem absolute value."
        "If this flag is set to true, the next parameter will be used as the range of NCs used to transform"
        "their values into the range of the experimental values. 1=true;0=false",
        io::Serialization::GetAgent( &m_SetNCRange),
        "0"
      );

      parameters.AddInitializer
      (
        "nc_range_min",
        "The minimimum value for the range of NCs (see parameter set_nc_range)",
        io::Serialization::GetAgent( &m_NCRangeMin),
        "0"
      );

      parameters.AddInitializer
      (
        "nc_range_max",
        "The maximum value for the range of NCs (see parameter set_nc_range)",
        io::Serialization::GetAgent( &m_NCRangeMax),
        "0"
      );

      return parameters;
    }

    //! @brief returns dummy name for ensemble
    //! @param ENSEMBLE for which a name will be created
    //! @return string which is the dummy name of ensemble
    std::string AnalyzeAccessibilityChange::EnsembleAsFilename
    (
      const assemble::ProteinEnsemble &ENSEMBLE
    )
    {
      return "start_ensemble_pdbs.ls";
    }

    //! @brief create ensemble from filename
    //! @param ENSEMBLE ensemble to setup
    //! @param NAME string name of file which will be used to create the ensemble
    //! @param ERR_STREAM stream to write out errors to
    //! @return ensemble created from the filename
    bool AnalyzeAccessibilityChange::EnsembleFromFilename( assemble::ProteinEnsemble &ENSEMBLE, const std::string &NAME, std::ostream &ERR_STREAM)
    {
      ENSEMBLE = assemble::ProteinEnsemble( NAME, 0, biol::GetAAClasses().e_AAComplete);
      return true;
    }

    //! @brief returns dummy name for exposure data file
    //! @return string which could be the name of the file the data comes from
    std::string AnalyzeAccessibilityChange::ExposureDataAsFilename
    (
      const storage::List< AccessibilityProfile> &DATA
    )
    {
      return "exposure_data.cst";
    }

    //! @brief reads in exposure data from a file given the filename
    //! @param PROFILES a list of accessibility profiles to set
    //! @param NAME the name of the file the data will be read from
    //! @param ERR_STREAM the stream any error will be written to
    //! @return true on success
    bool AnalyzeAccessibilityChange::ExposureDataFromFilename
    (
      storage::List< AccessibilityProfile> &PROFILES,
      const std::string &NAME,
      std::ostream &ERR_STREAM
    )
    {
      PROFILES.Reset();
      // open the data file
      io::IFStream read;
      io::File::MustOpenIFStream( read, NAME);

      // read in data file
      while( !read.eof() && read.peek() != std::istream::traits_type::eof())
      {
        // read in the chain and residue numbers
        char chain_a, chain_b;
        int aa_id_a, aa_id_b;
        io::Serialize::Read( chain_a, read);
        io::Serialize::Read( aa_id_a, read);
        io::Serialize::Read( chain_b, read);
        io::Serialize::Read( aa_id_b, read);

        // read in the exposure value
        double exposure;
        io::Serialize::Read( exposure, read);

        // create profile
        storage::List< AccessibilityAA> current_profile;
        // iterate through the residue range to make the profile
        for( int current_aa_id( aa_id_a); current_aa_id <= aa_id_b; ++current_aa_id)
        {
          // accessibility restraint
          AccessibilityAA restraint
          (
            storage::Map< AccessibilityAA::EnvironmentEnum, double>
            (
              storage::Map< AccessibilityAA::EnvironmentEnum, double>::Create
              (
                std::make_pair( AccessibilityAA::EnvironmentEnum( AccessibilityAA::e_Oxygen), exposure)
              )
            ),
            util::CloneToShPtr( assemble::LocatorAA( chain_a, current_aa_id, true)),
            util::ShPtr< assemble::AAExposureInterface>( new assemble::AANeighborCount())
          );

          // add current restraint
          current_profile.PushBack( restraint);
        }

        // create profile and add it to data
        const AccessibilityProfile profile( current_profile);
        PROFILES.PushBack( profile);

        // read end of line character
        char tmp;
        read.get( tmp);
      }

      BCL_MessageStd( "read in " + util::Format()( PROFILES.GetSize()) + " profiles");

      // return the data
      return true;
    }

    //! @brief determines the range of the calculated exposures and experimental exposures
    //! @param MEAN_DIFF_STATS to hold, over all profiles, the average of the average exposure difference calculated
    //!                        across all the residues in the profiles
    //! @return vector nd with two ranges, one for calculate the other for experimental exposures
    storage::VectorND< 2, math::Range< double> > AnalyzeAccessibilityChange::GetCalculatedAndExperimentalRanges
    (
      const assemble::ProteinEnsemble &ENSEMBLE, math::RunningAverageSD< double> &MEAN_DIFF_STATS
    ) const
    {
      // vectors to hold the values of calculated and experimental exposures
      storage::Vector< double> profile_mean;
      storage::Vector< double> profile_expval;

      // to hold, over all profiles, the average of the average exposure difference calculated across all the residues
      // in the profile
      MEAN_DIFF_STATS.Reset();

      // iterate through the data
      for
      (
        storage::List< AccessibilityProfile>::const_iterator profile_itr( m_ExperimentalExposures.Begin()),
        profile_itr_end( m_ExperimentalExposures.End());
        profile_itr != profile_itr_end;
        ++profile_itr
      )
      {
        // to hold, for a profile, the average exposure difference calculated across all the residues in the profile
        math::RunningAverageSD< double> average_profile_diff;

        // to hold the experimental data associated with this profile
        double exp_data;

        // iterate through the profile
        for
        (
          storage::List< AccessibilityAA>::const_iterator aa_itr( profile_itr->GetAccessibilities().Begin()),
          aa_itr_end( profile_itr->GetAccessibilities().End());
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {
          // reference to current data
          const AccessibilityAA &current_data( *aa_itr);

          // iterate over the start ensemble
          for
          (
            assemble::ProteinEnsemble::const_iterator
              start_ensemble_itr( m_StartEnsemble.Begin()), start_ensemble_itr_end( m_StartEnsemble.End());
            start_ensemble_itr != start_ensemble_itr_end;
            ++start_ensemble_itr
          )
          {
            // the accessibility assignment for the starting ensemble
            const AccessibilityAAAssignment assignment_start( current_data.GenerateAssignment( **start_ensemble_itr));

            // iterate over the end ensemble
            for
            (
              assemble::ProteinEnsemble::const_iterator
               end_ensemble_itr( ENSEMBLE.Begin()), end_ensemble_itr_end( ENSEMBLE.End());
              end_ensemble_itr != end_ensemble_itr_end;
              ++end_ensemble_itr
            )
            {
              // the accessibility assignment for the ending ensemble
              const AccessibilityAAAssignment assignment_end( current_data.GenerateAssignment( **end_ensemble_itr));

              // right now the data is just stored as oxygen accessibility
              // this is reset with each iteration but all exp data is the same across the entire profile
              exp_data = assignment_end.GetAccessibilityByEnvironment
              (
                AccessibilityAA::e_Oxygen
              ).Second();

              // add in the current exposure difference to the statistics calculator
              average_profile_diff += assignment_end.GetExposureValue() - assignment_start.GetExposureValue();
            }
          } // start ensemble
        } // profile aas

        // get the mean and add it and the experimental data to their vectors
        const double mean( average_profile_diff.GetAverage());
        profile_mean.PushBack( mean);
        profile_expval.PushBack( exp_data);

        // add the mean change for the current profile into the overall statistics mean
        if( !m_DirectRelation)
        {
          MEAN_DIFF_STATS += -mean;
        }
        else
        {
          MEAN_DIFF_STATS += mean;
        }
      } // profiles

      // get range of means
      const storage::Vector< double> &mean_values( profile_mean);
      const double mean_min( math::Absolute( math::Statistics::MinimumValue( mean_values.Begin(), mean_values.End())));
      const double mean_max( math::Absolute( math::Statistics::MaximumValue( mean_values.Begin(), mean_values.End())));
      const double mean_range_value( mean_min > mean_max ? mean_min : mean_max);
      math::Range< double> mean_range( -mean_range_value, mean_range_value);
      if( m_SetNCRange)
      {
        mean_range = math::Range< double>( m_NCRangeMin, m_NCRangeMax);
      }

      // get range of exp_data
      const storage::Vector< double> &exp_values( profile_expval);
      const double exp_min( math::Absolute( math::Statistics::MinimumValue( exp_values.Begin(), exp_values.End())));
      const double exp_max( math::Absolute( math::Statistics::MaximumValue( exp_values.Begin(), exp_values.End())));
      const double exp_range_value( exp_min > exp_max ? exp_min : exp_max);
      const math::Range< double> exp_range( -exp_range_value, exp_range_value);
      BCL_MessageStd
      (
        "calculated accessibilities will be transformed from range\n"
        + util::Format()( mean_range) + "\ninto range\n" + util::Format()( exp_range)
      );

      // holds the calculated and experimental ranges
      storage::VectorND< 2, math::Range< double> > ranges( mean_range, exp_range);

      // return ranges
      return ranges;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_analyze_atom_distance_heatmap.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "math/bcl_math_gnuplot_heatmap.h"
#include "math/bcl_math_histogram.h"
#include "restraint/bcl_restraint_epr_distance_data.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AnalyzeAtomDistanceHeatmap::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new AnalyzeAtomDistanceHeatmap())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AnalyzeAtomDistanceHeatmap::AnalyzeAtomDistanceHeatmap() :
      m_OutFilePostFix( ".AtomDistanceHeatmap"),
      m_ProteinModelDataType( EPRDistanceData::GetDefaultHandler()),
      m_HistogramMinimum( 10),
      m_HistogramBinSize( 2),
      m_HistogramNumberOfBins( 25),
      m_ShowAtomDistance( true)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AnalyzeAtomDistanceHeatmap
    AnalyzeAtomDistanceHeatmap *AnalyzeAtomDistanceHeatmap::Clone() const
    {
      return new AnalyzeAtomDistanceHeatmap( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AnalyzeAtomDistanceHeatmap::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &AnalyzeAtomDistanceHeatmap::GetOutFilePostfix() const
    {
      return m_OutFilePostFix;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &AnalyzeAtomDistanceHeatmap::GetAlias() const
    {
      static const std::string s_Name( "AtomDistanceHeatmap");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string AnalyzeAtomDistanceHeatmap::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // string to hold analysis
      std::string analysis;

      // get the protein model data from one of the models of the ensemble and make sure it could be cast
      util::ShPtrVector< AtomDistance> atom_distances( m_ProteinModelDataType->ReadRestraintsFromFile());

      // get the histograms of distances
      const storage::Vector< math::Histogram> score_histograms( GetDistanceHistograms( atom_distances, ENSEMBLE));

      // get the names of the restraints that will be used as tics
      const storage::Vector< std::string> restraint_tics( GetRestraintNameTics( atom_distances));

      // get the gnuplot script object
      math::GnuplotHeatmap heatmap( GetHeatMap( score_histograms, restraint_tics));

      // true if the atom distance restraints should be displayed
      if( m_ShowAtomDistance)
      {
        const storage::Vector< storage::VectorND< 2, storage::VectorND< 2, double> > > boxes
        (
          GetDistanceBoxes( atom_distances)
        );
        const double x_range( score_histograms.GetSize());
        const double x_min( 0);
        const double y_range( score_histograms.FirstElement().GetBoundaries().Second() - score_histograms.FirstElement().GetBoundaries().First());
        const double y_min( score_histograms.FirstElement().GetBoundaries().First());
        heatmap.SetBoxes( boxes, x_range, y_range, x_min, y_min);
      }

      // write heat map to string stream
      std::stringstream stream;
      heatmap.WriteScript( stream);

      // return the string
      return stream.str();
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AnalyzeAtomDistanceHeatmap::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Makes a heat map showing the frequency with which a distance occurs for each atom distance restraint."
        "Can optionally also show in the heat map the restraint distance with upper and lower bounds."
      );

      parameters.AddInitializer
      (
        "filename_postfix",
        "the postfix that will be appended to the output file in order to identify this analysis",
        io::Serialization::GetAgent( &m_OutFilePostFix),
        ".AtomDistanceHeatmap"
      );

      parameters.AddInitializer
      (
        "restraint",
        "the type of atom distance related restraint needed for analysis",
        io::Serialization::GetAgent( &m_ProteinModelDataType),
        EPRDistanceData::GetDefaultHandler().GetString()
      );

      parameters.AddInitializer
      (
        "histogram_minimum",
        "the minimal value representing the left boundary of the score histogram",
        io::Serialization::GetAgent( &m_HistogramMinimum),
        "-1.0"
      );

      parameters.AddInitializer
      (
        "histogram_binsize",
        "the width of one bin of the score histograms",
        io::Serialization::GetAgent( &m_HistogramBinSize),
        "0.1"
      );

      parameters.AddInitializer
      (
        "histogram_num_bins",
        "the number of bins in the score histograms",
        io::Serialization::GetAgent( &m_HistogramNumberOfBins),
        "10"
      );

      parameters.AddInitializer
      (
        "show_atom_distance",
        "one if the atom distance restraint should be indicated on the heatmap with upper and lower bounds - 0 otherwise",
        io::Serialization::GetAgent( &m_ShowAtomDistance),
        "1"
      );

      return parameters;
    }

    //! @brief creates histograms of how frequently a distance is observed in models for each restraint
    //! @param DATA the restraints that will be used to get distances
    //! @param ENSEMBLE the ensemble of models that will be used to get distances
    //! @return vector of histograms - one histogram for each restraint
    storage::Vector< math::Histogram> AnalyzeAtomDistanceHeatmap::GetDistanceHistograms
    (
      const util::ShPtrVector< AtomDistance> &DATA, const assemble::ProteinEnsemble &ENSEMBLE
    ) const
    {
      // hold the frequency of distance observed for each atom distance restraint
      storage::Vector< math::Histogram> distance_histograms
      (
        DATA.GetSize(),
        math::Histogram( m_HistogramMinimum, m_HistogramBinSize, m_HistogramNumberOfBins)
      );

      // iterator for the vector of distance histograms
      storage::Vector< math::Histogram>::iterator histogram_itr
      (
        distance_histograms.Begin()), histogram_itr_end( distance_histograms.End()
      );

      // iterate through the atom distance restraints
      for
      (
        util::ShPtrVector< AtomDistance>::const_iterator restraint_itr( DATA.Begin()), restraint_itr_end( DATA.End());
        restraint_itr != restraint_itr_end && histogram_itr != histogram_itr_end;
        ++restraint_itr, ++histogram_itr
      )
      {
        const DataPairwise &data_pair( ( *restraint_itr)->GetData());

        // get the distance from the ensemble fro the current data pair
        const storage::Vector< double> distances( ENSEMBLE.GetDistances( data_pair));

        // add the distances to the current histogram
        histogram_itr->CalculateHistogram( distances);

        histogram_itr->Normalize();
      }

      return distance_histograms;
    }

    //! @brief creates the tics that will be used in the heatmap - converts each restraint into a string to use as tic
    //! @param DATA the restraints that will be used to create tics
    //! @return vector of strings which are the tics representing each restraint
    storage::Vector< std::string> AnalyzeAtomDistanceHeatmap::GetRestraintNameTics
    (
      const util::ShPtrVector< AtomDistance> &DATA
    )
    {
      storage::Vector< std::string> tics;

      // iterate through the data to get the tics from the restraints
      for
      (
        util::ShPtrVector< AtomDistance>::const_iterator data_itr( DATA.Begin()), data_itr_end( DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        const assemble::LocatorAtomCoordinatesInterface &atom_locator_a( *( *data_itr)->GetData().First());
        const assemble::LocatorAtomCoordinatesInterface &atom_locator_b( *( *data_itr)->GetData().Second());
        tics.PushBack( assemble::LocatorAtomCoordinatesInterface::GetNameFromPair( atom_locator_a, atom_locator_b));
      }

      return tics;
    }

    //! @brief creates gnuplot heat map object from the distance histograms and the restraint names as tics
    //! @param HISTOGRAMS the histograms that will be used to make the heat map
    //! @param TICS the names of the restraints
    //! @return heat map object which represents the distribution of distances for each restraint
    math::GnuplotHeatmap AnalyzeAtomDistanceHeatmap::GetHeatMap
    (
      const storage::Vector< math::Histogram> &HISTOGRAMS, const storage::Vector< std::string> &TICS
    )
    {
      const util::SiPtrVector< const math::Histogram> histograms
      (
        util::ConvertToConstSiPtrVector( HISTOGRAMS)
      );

      math::GnuplotHeatmap heatmap;

      heatmap.SetFromHistograms( histograms, true, true);
      heatmap.SetTitleAndLabel( "Title", "restraint", "distance", "Frequency (fraction of models)");
      heatmap.SetTicsX( TICS, true, 1);
      heatmap.SetRotationXTics( 90);
      heatmap.SetFilename( "AnalyzeAtomDistanceHeatmap.gnuplot");
      heatmap.SetFont( "/usr/share/fonts/dejavu-lgc/DejaVuLGCSansMono.ttf", 10);

      return heatmap;
    }

    //! @brief gives the coordinates of boxes representing the atom distance restraint data
    //!        each atom distance restraint gives two boxes, one for the upper portion of the restraint, and
    //!        one for the lower portion of the restraint. The actual distance is given by the junction of the two
    //!        boxes
    //! @param DATA the restraints that will be used to create the boxes
    //! @return vector of [(x,y),(x,y)] coordinates - the lower left corner and upper right corner of the box,
    //!         respectively
    storage::Vector< storage::VectorND< 2, storage::VectorND< 2, double> > >
    AnalyzeAtomDistanceHeatmap::GetDistanceBoxes( const util::ShPtrVector< AtomDistance> &DATA)
    {
      storage::Vector< storage::VectorND< 2, storage::VectorND< 2, double> > > box_coords;

      // iterate through the atom distance restraint data
      for
      (
        util::ShPtrVector< AtomDistance>::const_iterator restraint_itr( DATA.Begin()), restraint_itr_end( DATA.End());
        restraint_itr != restraint_itr_end;
        ++restraint_itr
      )
      {
        const double distance( ( *restraint_itr)->GetDistance()->GetDistance());
        const double upper_bound( ( *restraint_itr)->GetUpperBound());
        const double lower_bound( ( *restraint_itr)->GetLowerBound());
        const double lower_left_x( restraint_itr - DATA.Begin());
        const double upper_right_x( restraint_itr - DATA.Begin() + 1);

        // lower bound box
        {
          const double lower_left_y( lower_bound);
          const double upper_right_y( distance);

          storage::VectorND< 2, storage::VectorND< 2, double> > corner_coords;
          corner_coords.First() = storage::VectorND< 2, double>( lower_left_x, lower_left_y);
          corner_coords.Second() = storage::VectorND< 2, double>( upper_right_x, upper_right_y);
          box_coords.PushBack( corner_coords);
        }

        // upper bound box
        {
          const double lower_left_y( distance);
          const double upper_right_y( upper_bound);

          storage::VectorND< 2, storage::VectorND< 2, double> > corner_coords;
          corner_coords.First() = storage::VectorND< 2, double>( lower_left_x, lower_left_y);
          corner_coords.Second() = storage::VectorND< 2, double>( upper_right_x, upper_right_y);
          box_coords.PushBack( corner_coords);
        }
      }

      return box_coords;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_analyze_atom_distance_mean_sd.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "restraint/bcl_restraint_epr_distance_data.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AnalyzeAtomDistanceMeanSD::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new AnalyzeAtomDistanceMeanSD())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AnalyzeAtomDistanceMeanSD::AnalyzeAtomDistanceMeanSD() :
      m_PrintRestaintDistance( false),
      m_OutFilePostFix( ".AnalyzeAtomDistanceMeanSD"),
      m_ProteinModelDataType( EPRDistanceData::GetDefaultHandler()),
      m_PrintAllModelDistances( false)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AnalyzeAtomDistanceMeanSD
    AnalyzeAtomDistanceMeanSD *AnalyzeAtomDistanceMeanSD::Clone() const
    {
      return new AnalyzeAtomDistanceMeanSD( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AnalyzeAtomDistanceMeanSD::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &AnalyzeAtomDistanceMeanSD::GetOutFilePostfix() const
    {
      return m_OutFilePostFix;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &AnalyzeAtomDistanceMeanSD::GetAlias() const
    {
      static const std::string s_Name( "AtomDistanceMeanSD");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string AnalyzeAtomDistanceMeanSD::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // string to hold analysis
      std::string analysis;

      util::ShPtrVector< AtomDistance> data( m_ProteinModelDataType->ReadRestraintsFromFile());

      static const size_t s_width( 17);
      {
        // method to format the entries in the output string
        static const util::Format s_format( util::Format().W( s_width));
        // add the first line of the analysis
        analysis += ( GetRestraintHeader( data, ENSEMBLE.GetNameFormatter(), s_format) + "\n");
      }

      static const util::Format s_format( util::Format().W( s_width).FFP( 3));
      // add the lines for the distances mean and std dev
      analysis += GetMeanSDAnalysis( data, ENSEMBLE, s_format, m_PrintAllModelDistances);

      // true if the restraint information should be given
      if( m_PrintRestaintDistance)
      {
        // add the restraint information to the output
        analysis += ( "\n" + GetRestraintInformation( data, ENSEMBLE.GetNameFormatter(), s_format));
      }

      return analysis;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AnalyzeAtomDistanceMeanSD::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_PrintRestaintDistance, ISTREAM);
      io::Serialize::Read( m_OutFilePostFix, ISTREAM);
      io::Serialize::Read( m_ProteinModelDataType, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AnalyzeAtomDistanceMeanSD::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_PrintRestaintDistance, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_OutFilePostFix, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_ProteinModelDataType, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AnalyzeAtomDistanceMeanSD::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        " calculates the mean and std deviation between residues of atom distance restraints in a protein ensemble"
        " For each atom distance of a restraint indicated by the desired protein model data type, the mean"
        " and standard deviation of the distance between the corresponding atoms is determined. Optionally,"
        " the information about the restraint experimental distance can be also output so that the ensemble"
        " can be compared to the experiment."
      );

      parameters.AddInitializer
      (
        "print_restraint_distance",
        "one if the restraint distance should be also be printed for comparison - 0 otherwise",
        io::Serialization::GetAgent( &m_PrintRestaintDistance),
        "1"
      );

      parameters.AddInitializer
      (
        "filename_postfix",
        "the postfix that will be appended to the output file in order to identify this analysis",
        io::Serialization::GetAgent( &m_OutFilePostFix),
        ".AnalyzeAtomDistanceMeanSD"
      );

      parameters.AddInitializer
      (
        "restraint",
        "the type of atom distance related restraint needed for analysis",
        io::Serialization::GetAgent( &m_ProteinModelDataType),
        EPRDistanceData::GetDefaultHandler().GetString()
      );

      parameters.AddInitializer
      (
        "print_all_models_distances",
        "one if the distances in each of the models for each of the restraints should be printed - 0 otherwise",
        io::Serialization::GetAgent( &m_PrintAllModelDistances),
        "1"
      );

      return parameters;
    }

    //! @brief prints top line of output which is the each of the restraints
    //! @param DATA the list of atom distance objects which are the restraints
    //! @param LINE_NAME_FORMAT the format to use to format the name of the line i.e. the first column
    //! @param FORMAT the format object to format the restraints names
    //! @return string which has the list of restraints in a single line
    std::string AnalyzeAtomDistanceMeanSD::GetRestraintHeader
    (
      const util::ShPtrVector< AtomDistance> &DATA, const util::Format &LINE_NAME_FORMAT, const util::Format &FORMAT
    )
    {

      std::string header_string( LINE_NAME_FORMAT( "AtomPairs"));

      // iterate through the restraints
      for
      (
        util::ShPtrVector< AtomDistance>::const_iterator data_itr( DATA.Begin()), data_itr_end( DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        const assemble::LocatorAtomCoordinatesInterface &atom_locator_a( *( *data_itr)->GetData().First());
        const assemble::LocatorAtomCoordinatesInterface &atom_locator_b( *( *data_itr)->GetData().Second());

        header_string += FORMAT
          (
            assemble::LocatorAtomCoordinatesInterface::GetNameFromPair( atom_locator_a, atom_locator_b)
          );
      }

      return header_string;
    }

    //! @brief gives the mean and standard deviations for each restraint in formatted string
    //! @param DATA the list of atom distance objects which are the restraints
    //! @PARAM ENSEMBLE the ensemble from which distance mean and std devs will be calculated
    //! @param FORMAT the format object to format the means and std devs
    //! @return string which has the mean and standard deviations for each restraint - means on one line, sd on next
    std::string AnalyzeAtomDistanceMeanSD::GetMeanSDAnalysis
    (
      const util::ShPtrVector< AtomDistance> &DATA, const assemble::ProteinEnsemble &ENSEMBLE,
      const util::Format &FORMAT, const bool PRINT_ALL_MODEL_DISTANCES
    )
    {

      storage::Vector< std::string> protein_names( ENSEMBLE.GetPDBNames());
      const util::Format name_formatter( ENSEMBLE.GetNameFormatter());
      std::string mean_string( name_formatter( "Mean"));
      std::string sd_string( name_formatter( "StdDev"));

      // iterate through the restraint data to get the mean and std dev of each in the ensemble
      for
      (
        util::ShPtrVector< AtomDistance>::const_iterator data_itr( DATA.Begin()), data_itr_end( DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        math::RunningAverageSD< double> mean_sd;

        storage::Vector< std::string>::iterator
          name_itr( protein_names.Begin()), name_itr_end( protein_names.End());
        // iterate through the ensemble to build up the statistics
        for
        (
          assemble::ProteinEnsemble::const_iterator ensemble_itr( ENSEMBLE.Begin()), ensemble_itr_end( ENSEMBLE.End());
          ensemble_itr != ensemble_itr_end && name_itr != name_itr_end;
          ++ensemble_itr, ++name_itr
        )
        {
          const double current_distance( ( *data_itr)->GetData().EuclidianDistance( **ensemble_itr));
          //static const util::Format s_format( util::Format().W( s_width).FFP( 3));
          ( *name_itr) += FORMAT( current_distance);

          mean_sd += current_distance;

        }

        //const math::RunningAverageSD< double> mean_sd( ENSEMBLE.GetDistanceStatistics( ( *data_itr)->GetData()));

        mean_string += FORMAT( mean_sd.GetAverage());
        sd_string += FORMAT( mean_sd.GetStandardDeviation());
      }

      // holds the mean and std dev lines
      std::string analysis;

      if( PRINT_ALL_MODEL_DISTANCES)
      {
        for
        (
          storage::Vector< std::string>::const_iterator
            name_itr( protein_names.Begin()), name_itr_end( protein_names.End());
          name_itr != name_itr_end;
          ++name_itr
        )
        {
          analysis += ( *name_itr) + "\n";
        }
      }

       analysis += mean_string + "\n" + sd_string;

      return analysis;
    }

    //! @brief prints top line of output which is the each of the restraints
    //! @param DATA the list of atom distance objects which are the restraints
    //! @param FORMAT the format object to format the restraints
    //! @return string which has the list of restraints in a single line
    std::string AnalyzeAtomDistanceMeanSD::GetRestraintInformation
    (
      const util::ShPtrVector< AtomDistance> &DATA, const util::Format &LINE_NAME_FORMAT, const util::Format &FORMAT
    )
    {
      std::string distance_string( LINE_NAME_FORMAT( "CstDist"));
      std::string lower_bound_string( LINE_NAME_FORMAT( "CstLowBnd"));
      std::string upper_bound_string( LINE_NAME_FORMAT( "CstUpBnd"));

      // iterate through the restraint data to get the mean and std dev of each in the ensemble
      for
      (
        util::ShPtrVector< AtomDistance>::const_iterator data_itr( DATA.Begin()), data_itr_end( DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        distance_string    += FORMAT( ( *data_itr)->GetDistance()->GetDistance());
        lower_bound_string += FORMAT( ( *data_itr)->GetLowerBound());
        upper_bound_string += FORMAT( ( *data_itr)->GetUpperBound());
      }

      // holds the mean and std dev lines
      const std::string analysis( distance_string + "\n" + lower_bound_string + "\n" + upper_bound_string);

      return analysis;
    }

  } // namespace restraint

} // namespace bcl
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
#include "restraint/bcl_restraint_analyze_atom_distance_pymol.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "io/bcl_io_directory_entry.h"
#include "restraint/bcl_restraint_atom_distance_assignment.h"
#include "restraint/bcl_restraint_epr_distance_data.h"
#include "util/bcl_util_color_gradient.h"
#include "util/bcl_util_colors.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AnalyzeAtomDistancePymol::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new AnalyzeAtomDistancePymol())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AnalyzeAtomDistancePymol::AnalyzeAtomDistancePymol() :
      m_OutFilePostFix( ".AnalyzeAtomDistancePymol.py"),
      m_Score(),
      m_ColorGradient
      (
        util::ColorGradient
        (
          math::Range< double>( -1, 0),
          storage::Vector< util::Color>::Create( util::GetColors().e_Blue, util::GetColors().e_Red)
        )
      ),
      m_LongDistanceWideLine( false),
      m_EnsembleRepresentative( 0),
      m_UseCA( false),
      m_RestraintType( EPRDistanceData::GetDefaultHandler())
    {
    }

    //! @brief Clone function
    //! @return pointer to new AnalyzeAtomDistancePymol
    AnalyzeAtomDistancePymol *AnalyzeAtomDistancePymol::Clone() const
    {
      return new AnalyzeAtomDistancePymol( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AnalyzeAtomDistancePymol::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &AnalyzeAtomDistancePymol::GetOutFilePostfix() const
    {
      return m_OutFilePostFix;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &AnalyzeAtomDistancePymol::GetAlias() const
    {
      static const std::string s_name( "AtomDistancePymol");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string AnalyzeAtomDistancePymol::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // to hold the analysis string
      std::string analysis;

      // get the necessary python includes for the script
      analysis += GetScriptHeader();

      // get command for loading pdb model of interest
      analysis += GetLoadPDBCommand( ENSEMBLE);

      // get the commands for showing the distances
      analysis += GetDistanceLines( ENSEMBLE);

      // return the analysis string
      return analysis;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AnalyzeAtomDistancePymol::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_OutFilePostFix, ISTREAM);
      io::Serialize::Read( m_Score, ISTREAM);
      io::Serialize::Read( m_ColorGradient, ISTREAM);
      io::Serialize::Read( m_LongDistanceWideLine, ISTREAM);
      io::Serialize::Read( m_EnsembleRepresentative, ISTREAM);
      io::Serialize::Read( m_UseCA, ISTREAM);
      io::Serialize::Read( m_RestraintType, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AnalyzeAtomDistancePymol::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_OutFilePostFix, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_Score, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_ColorGradient, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_LongDistanceWideLine, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_EnsembleRepresentative, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_UseCA, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_RestraintType, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AnalyzeAtomDistancePymol::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Writes pymol script to show atom distance restraints in pymol."
        "Can optionally color the restraint lines by gradient according to their score if a score object is"
        "provided. Can optionally indicate whether the distance in the protein model is too long or short by"
        "making the restraint line wide in the former and thinner in the latter."
      );

      parameters.AddInitializer
      (
        "filename_postfix",
        "the postfix that will be appended to the output file in order to identify this analysis",
        io::Serialization::GetAgent( &m_OutFilePostFix),
        ".AnalyzeAtomPymol.py"
      );

      parameters.AddInitializer
      (
        "restraint_type",
        "the atom distance restraint type where the score will come from",
        io::Serialization::GetAgent( &m_Score),
        EPRDistanceData().GetAlias()
      );

      parameters.AddInitializer
      (
        "color_gradient",
        "the color gradient that should be used to color the restraint lines",
        io::Serialization::GetAgent( &m_ColorGradient)
      );

      parameters.AddInitializer
      (
        "distance_affects_line_width",
        "one if the width of the line should be wide if the distance is too long and thinner if distance is too short",
        io::Serialization::GetAgent( &m_LongDistanceWideLine),
        "1"
      );

      parameters.AddInitializer
      (
        "ensemble_representative",
        "the model in the ensemble that should be loaded in pymol as the representative of the ensemble - start at 0",
        io::Serialization::GetAgent( &m_EnsembleRepresentative),
        "0"
      );

      parameters.AddInitializer
      (
        "use_ca",
        "one to use CA atoms for drawing the lines, zero to use the locator atom types",
        io::Serialization::GetAgent( &m_UseCA),
        "1"
      );

      parameters.AddInitializer
      (
        "restraint",
        "the type of atom distance related restraint needed for analysis",
        io::Serialization::GetAgent( &m_RestraintType),
        EPRDistanceData::GetDefaultHandler().GetString()
      );

      return parameters;
    }

    //! @brief gives the text that is necessary at the top of the script file for it to work
    //! @return string which has the necessary text
    std::string AnalyzeAtomDistancePymol::GetScriptHeader()
    {
      return "from pymol import cmd\nimport string\nfrom pymol.cgo import *\nfrom pymol.vfont import plain\n";
    }

    //! @brief gives the command to load the desired representative model from an ensemble
    //! @param ENSEMBLE the ensemble from which the desired protein model will be loaded
    //! @return string the string which is the command needed to load the desired protein model from the ensemble
    std::string AnalyzeAtomDistancePymol::GetLoadPDBCommand( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure the desired model index within range of ensemble size
      BCL_Assert
      (
        m_EnsembleRepresentative < ENSEMBLE.GetSize(),
        "ensemble representative number is " + util::Format()( m_EnsembleRepresentative)
        + " but the size of the ensemble is " + util::Format()( ENSEMBLE.GetSize()) + ". First model in ensemble is "
        "model 0."
      );

      // get the desired protein from the ensemble
      assemble::ProteinEnsemble::const_iterator ensemble_itr( ENSEMBLE.Begin());
      const assemble::ProteinModel &model
      (
        **storage::AdvanceIterator( ensemble_itr, ENSEMBLE.End(), m_EnsembleRepresentative)
      );

      // get the protein model data from the model and make sure it could be cast
      const util::ShPtr< util::Wrapper< std::string> > data
      (
        model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile)
      );
      BCL_Assert( data.IsDefined(), "could not get pdb name from protein model");

      // the pdb filename
      const io::DirectoryEntry pdb_filename( data->GetData());

      // pymol python command
      const std::string command( "cmd.load( \"" + pdb_filename.GetFullName() + "\")\n");

      // return load command
      return command;
    }

    //! @brief gives the commands necessary to display the distance lines in the desired manner in pymol
    //! @param ENSEMBLE the ensemble from which distances are going to be shown
    //! @return string which has the commands necessary to display the distance lines in the desired manner in pymol
    std::string AnalyzeAtomDistancePymol::GetDistanceLines( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // hold the commands
      std::string analysis;

      // get the protein model data from one of the models of the ensemble and make sure it could be cast
      util::ShPtrVector< AtomDistance> data( m_RestraintType->ReadRestraintsFromFile());

      // get the scoring statistics
      const storage::Vector< math::RunningAverageSD< double> > score_stats( GetScoreStatistics( data, ENSEMBLE));

      // get the distance statistics
      const storage::Vector< math::RunningAverageSD< double> > distance_stats( GetDistanceStatistics( data, ENSEMBLE));

      // iterate through the restraint data and statistics
      util::ShPtrVector< AtomDistance>::const_iterator restraint_itr( data.Begin()), restraint_itr_end( data.End());
      for
      (
        storage::Vector< math::RunningAverageSD< double> >::const_iterator
          score_stat_itr( score_stats.Begin()), score_stat_itr_end( score_stats.End()),
          distance_stat_itr( distance_stats.Begin()), distance_stat_itr_end( distance_stats.End());
        score_stat_itr != score_stat_itr_end && restraint_itr != restraint_itr_end &&
          distance_stat_itr != distance_stat_itr_end;
        ++score_stat_itr, ++restraint_itr, ++distance_stat_itr
      )
      {
        // get the command for the current restraint distance line
        const std::string command
        (
          GetLineCommand( **restraint_itr, *score_stat_itr, *distance_stat_itr)
        );

        // add current distance command to analysis
        analysis += command;
      }

      // return analysis string
      return analysis;
    }

    //! @brief collects the score statistics of each atom distance restraint across the ensemble
    //! @param DATA the list of atom distance objects which are the restraints
    //! @param ENSEMBLE the ensemble from which distance mean and std devs will be calculated
    //! @return vector of statistics object holding the mean and std dev of scores for the restraints in the ensemble
    storage::Vector< math::RunningAverageSD< double> > AnalyzeAtomDistancePymol::GetScoreStatistics
    (
      const util::ShPtrVector< AtomDistance> &DATA, const assemble::ProteinEnsemble &ENSEMBLE
    ) const
    {
      // to hold statistics of scores
      storage::Vector< math::RunningAverageSD< double> > statistics;

      // the scorer of individual restraints
      const math::FunctionInterfaceSerializable< AtomDistanceAssignment, double> &scorer( *m_Score);

      // iterate through the restraint data
      for
      (
        util::ShPtrVector< AtomDistance>::const_iterator
          restraint_itr( DATA.Begin()), restraint_itr_end( DATA.End());
        restraint_itr != restraint_itr_end;
        ++restraint_itr
      )
      {
        const AtomDistance &atom_distance( **restraint_itr);

        // to hold the score statistics of the current restraint
        math::RunningAverageSD< double> score_stats;

        // iterate through the ensemble
        for
        (
          assemble::ProteinEnsemble::const_iterator
            ensemble_itr( ENSEMBLE.Begin()), ensemble_itr_end( ENSEMBLE.End());
          ensemble_itr != ensemble_itr_end;
          ++ensemble_itr
        )
        {
          // generate the assignment and get the score
          const AtomDistanceAssignment assignment( atom_distance.GenerateAssignment( **ensemble_itr));
          const double current_score( scorer( assignment));

          // add the score to the statistics object
          score_stats += current_score;
        }

        // add the statistics for the current restraint
        statistics.PushBack( score_stats);
      }

      return statistics;
    }

    //! @brief collects the distance statistics of each atom distance restraint across the ensemble
    //! @param DATA the list of atom distance objects which are the restraints
    //! @param ENSEMBLE the ensemble from which distance mean and std devs will be calculated
    //! @return vector of statistics object holding the mean and std dev of distances for restraints in the ensemble
    storage::Vector< math::RunningAverageSD< double> > AnalyzeAtomDistancePymol::GetDistanceStatistics
    (
      const util::ShPtrVector< AtomDistance> &DATA, const assemble::ProteinEnsemble &ENSEMBLE
    ) const
    {
      // to hold statistics of distances
      storage::Vector< math::RunningAverageSD< double> > statistics;

      // iterate through the restraint data to get the mean and std dev of each in the ensemble
      for
      (
        util::ShPtrVector< AtomDistance>::const_iterator data_itr( DATA.Begin()), data_itr_end( DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // calculate the mean and standard deviation of distance in ensemble
        const math::RunningAverageSD< double> mean_sd( ENSEMBLE.GetDistanceStatistics( ( *data_itr)->GetData()));

        // add to vector of statistics
        statistics.PushBack( mean_sd);
      }

      // return the vector of statistics
      return statistics;
    }

    //! @brief gives the command for showing a restraint as a line connecting atoms in a protein model
    //! @param RESTRAINT the restraint that will be shown
    //! @param SCORE_STATS the score statistics for the current restraint
    //! @param DISTANCE_STATS the distance statistics for the current restraint
    //! @return string which has the commands to show the restraint in pymol as line connecting restraint atoms
    std::string AnalyzeAtomDistancePymol::GetLineCommand
    (
      const AtomDistance &RESTRAINT, const math::RunningAverageSD< double> &SCORE_STATS,
      const math::RunningAverageSD< double> &DISTANCE_STATS
    ) const
    {
      // hold the python commands
      std::string command;

      // name of the restraint in pymol
      const std::string name
      (
        assemble::LocatorAtomCoordinatesInterface::GetNameFromPair
        (
          *RESTRAINT.GetData().First(), *RESTRAINT.GetData().Second()
        )
      );

      // add commands
      command += GetDistanceLineCommand( RESTRAINT.GetData(), name) + "\n";

      // true if distance line should be colored by score
      if( m_Score.IsDefined())
      {
        BCL_Assert( m_ColorGradient.IsDefined(), "the scoring object is defined but color gradient object is not");
        command += GetColorCommand( SCORE_STATS, name) + "\n";
      }

      // true if the width of the distance line should be altered - wider for distances that are longer than restraint
      if( m_LongDistanceWideLine)
      {
        command += GetLineWidthCommand( *RESTRAINT.GetDistance(), DISTANCE_STATS, name) + "\n";
      }

      // set dash width to zero so lines are solid
      command += "cmd.set( \"dash_gap\", 0, \"" + name + "\")\n";

      return command;
    }

    //! @brief gives the command to make a distance object in pymol
    //! @param DATA the pair of points the distance will be between
    //! @param NAME the name of the created distance object
    //! @return string which is the command to make a distance object in pymol between the points in data
    std::string
    AnalyzeAtomDistancePymol::GetDistanceLineCommand( const DataPairwise &DATA, const std::string &NAME) const
    {
      // pymol command for first selection
      const std::string atom_a( GetAtomSelection( *DATA.First()));

      // pymol command for second selection
      const std::string atom_b( GetAtomSelection( *DATA.Second()));

      // command to make a distance line
      const std::string command( "cmd.distance( \"" + NAME + "\", \"" + atom_a + "\",\"" + atom_b + "\")");

      // return command
      return command;
    }

    //! @brief gives the string to indicate an atom within pymol
    //! @param LOCATOR locator which will be indicated
    //! @return string that can be used to indicate an atom within pymol
    std::string
    AnalyzeAtomDistancePymol::GetAtomSelection( const assemble::LocatorAtomCoordinatesInterface &LOCATOR) const
    {
      // assume ca is used
      std::string atom_name( "CA");

      // get the atom type from the locator
      const biol::AtomType atom_type( LOCATOR.GetAtomType());

      // true if the atom type is defined and don't want to force use of CA
      if( atom_type.IsDefined() && !m_UseCA)
      {
        // atom name is the type from the locator
        atom_name = atom_type.GetName();
      }

      // pymol syntax for a selection
      std::string selection
      (
        "chain " + util::Format()( LOCATOR.GetChainID()) + " and resi " + util::Format()( LOCATOR.GetSeqID()) +
        " and name " + atom_name
      );

      // return selection string
      return selection;
    }

    //! @brief gives the command to set the line of the distance restraint to the correct color
    //! @param STATISTICS the score statistics that will determine the current color
    //! @param NAME the name of the current retraint selection
    //! @return string with the commands to set the restraint line to the correct color
    std::string AnalyzeAtomDistancePymol::GetColorCommand
    (
      const math::RunningAverageSD< double> &STATISTICS, const std::string &NAME
    ) const
    {
      // to hold the command string
      std::string command;

      // defined color name
      command += "color = \"color_" + NAME + "\"\n";

      // calculate color value
      const linal::Vector3D color( m_ColorGradient->operator ()( STATISTICS.GetAverage()));

      // set color value
      command += "cmd.set_color( color, [ " + util::Format()( color.X()) +
        "," + util::Format()( color.Y()) + "," + util::Format()( color.Z()) + "])\n";

      // assign color to line based on gradient
      command += "cmd.set( \"dash_color\", color,\"" + NAME + "\")\n";

      // return the command
      return command;
    }

    //! @brief gives the commands to set the with of the restraint line
    //! @param DISTANCE the distance object giving the restraint distance
    //! @param STATISTICS the statistics object giving the average model distance
    //! @param NAME the name of the current retraint selection
    std::string AnalyzeAtomDistancePymol::GetLineWidthCommand
    (
      const Distance &DISTANCE, const math::RunningAverageSD< double> &STATISTICS, const std::string &NAME
    ) const
    {
      // initial width of distance line
      double dash_width( 2);

      // messages
      BCL_MessageStd( "DISTANCE.GetDistance() " + util::Format()( DISTANCE.GetDistance()));
      BCL_MessageStd( "STATISTICS.GetAverage() " + util::Format()( STATISTICS.GetAverage()));

      // true if the average distance of the models is larger than the experimental distance
      if( m_LongDistanceWideLine && DISTANCE.GetDistance() < STATISTICS.GetAverage())
      {
        dash_width = 8;
      }

      // command for setting line width in pymol
      const std::string command( "cmd.set( \"dash_width\"," + util::Format()( dash_width) + ", \"" + NAME + "\")");

      // return command
      return command;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_analyze_atom_distance_score.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "restraint/bcl_restraint_analyze_atom_distance_mean_sd.h"
#include "restraint/bcl_restraint_epr_distance_data.h"
#include "score/bcl_score_restraint_distance_spin_label.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AnalyzeAtomDistanceScore::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new AnalyzeAtomDistanceScore())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AnalyzeAtomDistanceScore::AnalyzeAtomDistanceScore() :
      m_Score(), // EPRDistanceData().GetAlias()
      m_PrintAllModelScores( true),
      m_OutFilePostFix( ".AnalyzeAtomDistanceScore"),
      m_PrintStatistics( true),
      m_RestraintType( EPRDistanceData::GetDefaultHandler())
    {
    }

    //! @brief Clone function
    //! @return pointer to new AnalyzeAtomDistanceScore
    AnalyzeAtomDistanceScore *AnalyzeAtomDistanceScore::Clone() const
    {
      return new AnalyzeAtomDistanceScore( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AnalyzeAtomDistanceScore::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &AnalyzeAtomDistanceScore::GetOutFilePostfix() const
    {
      return m_OutFilePostFix;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &AnalyzeAtomDistanceScore::GetAlias() const
    {
      static const std::string s_Name( "AtomDistanceScore");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string AnalyzeAtomDistanceScore::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // string to hold analysis
      std::string analysis;

      // get the protein model data from one of the models of the ensemble and make sure it could be cast
      util::ShPtrVector< AtomDistance> data( m_RestraintType->ReadRestraintsFromFile());

      // format for the line names
      static const size_t s_width( 17);
      util::Format line_name_format( util::Format().W( s_width));

      // true if scores for each model will be printed -  need to get length of pdb filename
      if( m_PrintAllModelScores)
      {
        // get the name of the first model in order to determine its length
        util::ShPtr< util::Wrapper< std::string> > pdb_name
        (
          ( *ENSEMBLE.Begin())->GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile)
        );
        BCL_Assert( pdb_name.IsDefined(), "could not get pdb name from protein model");

        line_name_format = util::Format().W( pdb_name->length());
      }

      static const util::Format s_format( util::Format().W( s_width).FFP( 3));
      // add the first line of the analysis
      {
        analysis += AnalyzeAtomDistanceMeanSD::GetRestraintHeader( data, line_name_format, s_format);
        analysis += ( s_format( "sum") + s_format( "total_score") + "\n");
      }

      // add the score analysis lines
      analysis += GetScoreAnalysis( data, ENSEMBLE, line_name_format, s_format);

      return analysis;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AnalyzeAtomDistanceScore::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Score, ISTREAM);
      io::Serialize::Read( m_PrintAllModelScores, ISTREAM);
      io::Serialize::Read( m_OutFilePostFix, ISTREAM);
      io::Serialize::Read( m_PrintStatistics, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AnalyzeAtomDistanceScore::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Score, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_PrintAllModelScores, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_OutFilePostFix, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_PrintStatistics, OSTREAM, INDENT) ;

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AnalyzeAtomDistanceScore::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates atom distance scores for ensemble of models and gives mean and sd of score of each restraint."
        "Can also optionally write the individual restraint score for every model in the ensemble"
      );

      parameters.AddInitializer
      (
        "score",
        "the atom distance score",
        io::Serialization::GetAgent( &m_Score),
        score::RestraintDistanceSpinLabel::GetDefaultScheme()
      );

      parameters.AddInitializer
      (
        "print_all_scores",
        "one if scores should be printed for every model - 0 otherwise",
        io::Serialization::GetAgent( &m_PrintAllModelScores),
        "1"
      );

      parameters.AddInitializer
      (
        "filename_postfix",
        "the postfix that will be appended to the output file in order to identify this analysis",
        io::Serialization::GetAgent( &m_OutFilePostFix),
        ".AnalyzeAtomDistanceScore"
      );

      parameters.AddInitializer
      (
        "print_statistics",
        "one if the mean, std dev, min, and max statistics should be printed for each atom distance restraint - 0 otherwise",
        io::Serialization::GetAgent( &m_PrintStatistics),
        "1"
      );

      parameters.AddInitializer
      (
        "restraint_type",
        "actual type of restraint that holds the data",
        io::Serialization::GetAgent( &m_RestraintType),
        EPRDistanceData().GetAlias()
      );
      return parameters;
    }

    //! @brief gives the string that is the scoring analysis
    //! @param DATA the list of atom distance objects which are the restraints
    //! @param ENSEMBLE the ensemble from which distance mean and std devs will be calculated
    //! @param LINE_NAME_FORMAT the format to use to format the name of the line i.e. the first column
    //! @param FORMAT the format object to format the scores
    std::string AnalyzeAtomDistanceScore::GetScoreAnalysis
    (
      const util::ShPtrVector< AtomDistance> &DATA, const assemble::ProteinEnsemble &ENSEMBLE,
      const util::Format &LINE_NAME_FORMAT, const util::Format &FORMAT
    ) const
    {
      std::string analysis_string;

      // get statistics of atom distance restraint scores and add formatted analysis for each model if desired
      storage::Pair
      <
        storage::Vector< math::RunningAverageSD< double> >, storage::Vector< math::RunningMinMax< double> >
      > statistics( GetScoreStatistics( DATA, ENSEMBLE, LINE_NAME_FORMAT, FORMAT, analysis_string));

      if( m_PrintStatistics)
      {
        analysis_string +=
        (
          "\n" + GetScoreStatisticsFormattedAnalysis
          (
            statistics.First(), statistics.Second(), LINE_NAME_FORMAT, FORMAT
          )
        );
      }

      return analysis_string;
    }

    //! @brief gives statistics for each atom distance restraint and add individual model analysis if necessary
    //! @param DATA the list of atom distance objects which are the restraints
    //! @PARAM ENSEMBLE the ensemble from which distance mean and std devs will be calculated
    //! @param LINE_NAME_FORMAT the format to use to format the name of the line i.e. the first column
    //! @param FORMAT the format object to format the scores
    //! @param ANALYSIS_STRING string holding formatted analysis
    //! @return pair of vectors with a statistic object for each atom distance restraint
    storage::Pair
    <
      storage::Vector< math::RunningAverageSD< double> >, storage::Vector< math::RunningMinMax< double> >
    > AnalyzeAtomDistanceScore::GetScoreStatistics
    (
      const util::ShPtrVector< AtomDistance> &DATA, const assemble::ProteinEnsemble &ENSEMBLE,
      const util::Format &LINE_NAME_FORMAT, const util::Format &FORMAT, std::string &ANALYSIS_STRING
    ) const
    {
      // takes into account adding columns for the total score and the simple sum of the individual restraint scores
      static const size_t s_additional_columns( 2);

      // to hold mean and sd of scores
      storage::Vector< math::RunningAverageSD< double> > mean_sd
      (
        DATA.GetSize() + s_additional_columns, math::RunningAverageSD< double>()
      );

      // to hold min and max of scores
      storage::Vector< math::RunningMinMax< double> > min_max
      (
        DATA.GetSize() + s_additional_columns, math::RunningMinMax< double>()
      );

      // score to score the individual restraints
      const math::FunctionInterfaceSerializable< AtomDistanceAssignment, double> &individual_score( *m_Score->GetScore());

      // iterate through the ensemble
      for
      (
        assemble::ProteinEnsemble::const_iterator protein_itr( ENSEMBLE.Begin()), protein_itr_end( ENSEMBLE.End());
        protein_itr != protein_itr_end;
        ++protein_itr
      )
      {
        const assemble::ProteinModel &model( **protein_itr);

        // iterators for the statistics vectors
        storage::Vector< math::RunningAverageSD< double> >::iterator
          mean_sd_itr( mean_sd.Begin()), mean_sd_itr_end( mean_sd.End());
        storage::Vector< math::RunningMinMax< double> >::iterator
          min_max_itr( min_max.Begin()), min_max_itr_end( min_max.End());

        // true if scores for each model will be printed -  need to print of pdb filename
        if( m_PrintAllModelScores)
        {
          // get the name of the first model in order to determine its length
          util::ShPtr< util::Wrapper< std::string> > pdb_name
          (
            ( *protein_itr)->GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile)
          );
          BCL_Assert( pdb_name.IsDefined(), "could not get atom distance restraints from protein model");
          if( protein_itr != ENSEMBLE.Begin())
          {
            ANALYSIS_STRING += "\n";
          }
          ANALYSIS_STRING += LINE_NAME_FORMAT( std::string( *pdb_name));
        }

        // to hold sum of restraint scores for this protein
        double score_sum( 0);

        // iterate through the atom distance restraints
        for
        (
          util::ShPtrVector< AtomDistance>::const_iterator restraint_itr( DATA.Begin()), restraint_itr_end( DATA.End());
          restraint_itr != restraint_itr_end && mean_sd_itr != mean_sd_itr_end && min_max_itr != min_max_itr_end;
          ++restraint_itr, ++mean_sd_itr, ++min_max_itr
        )
        {
          // get the assignment
          const AtomDistanceAssignment assignment( ( *restraint_itr)->GenerateAssignment( model));

          // get the current score
          const double score( individual_score( assignment));

          score_sum += score;

          // add the score to the two statistics objects
          *mean_sd_itr += score;
          *min_max_itr += score;

          // true if scores for each model will be printed -  need to print current score
          if( m_PrintAllModelScores)
          {
            ANALYSIS_STRING += FORMAT( score);
          }
        } // end iterate through restraints

        // add the score sum and total score to the two statistics objects
        const double score_total( m_Score->operator()( model));
        *mean_sd_itr += score_sum;
        *++mean_sd_itr += score_total;
        *min_max_itr += score_sum;
        *++min_max_itr += score_total;

        // true if scores for each model will be printed -  need to print score sum and total score
        if( m_PrintAllModelScores)
        {
          ANALYSIS_STRING += FORMAT( score_sum);
          ANALYSIS_STRING += FORMAT( score_total);
        }
      } // end iterate through ensemble

      return storage::Pair
        <
          storage::Vector< math::RunningAverageSD< double> >, storage::Vector< math::RunningMinMax< double> >
        >( mean_sd, min_max);
    }

    //! @brief gives formatted analysis of the statistics of the atom distance restraints
    //! @param MEAN_SD the mean and standard deviation statistics of the atom distance restraints
    //! @param MIN_MAX the min and max statistics of the atom distance restraints
    //! @param LINE_NAME_FORMAT the format to use to format the name of the line i.e. the first column
    //! @param FORMAT the format object to format the scores
    //! @return string which has the formatted statistics analysis of atom restraint scores
    std::string AnalyzeAtomDistanceScore::GetScoreStatisticsFormattedAnalysis
    (
      const storage::Vector< math::RunningAverageSD< double> > &MEAN_SD,
      const storage::Vector< math::RunningMinMax< double> > &MIN_MAX,
      const util::Format &LINE_NAME_FORMAT, const util::Format &FORMAT
    ) const
    {
      std::string analysis;
      analysis += ( GetMeanSDLines( MEAN_SD, LINE_NAME_FORMAT, FORMAT) + "\n");
      analysis += GetMinMaxLines( MIN_MAX, LINE_NAME_FORMAT, FORMAT);

      return analysis;
    }

    //! @brief gives formatted analysis of the mean and sd statistics of the atom distance restraints
    //! @param MEAN_SD the mean and standard deviation statistics of the atom distance restraints
    //! @param LINE_NAME_FORMAT the format to use to format the name of the line i.e. the first column
    //! @param FORMAT the format object to format the scores
    //! @return string which has the formatted mean and sd statistics analysis of atom restraint scores
    std::string AnalyzeAtomDistanceScore::GetMeanSDLines
    (
      const storage::Vector< math::RunningAverageSD< double> > &MEAN_SD,
      const util::Format &LINE_NAME_FORMAT, const util::Format &FORMAT
    )
    {
      std::string mean_line( LINE_NAME_FORMAT( "Mean"));
      std::string sdev_line( LINE_NAME_FORMAT( "StdDev"));

      // iterate through the mean and std dev statistics
      for
      (
        storage::Vector< math::RunningAverageSD< double> >::const_iterator
          stat_itr( MEAN_SD.Begin()), stat_itr_end( MEAN_SD.End());
        stat_itr != stat_itr_end; ++stat_itr
      )
      {
        mean_line += FORMAT( stat_itr->GetAverage());
        sdev_line += FORMAT( stat_itr->GetStandardDeviation());
      }

      std::string analysis( mean_line + "\n" + sdev_line);

      return analysis;
    }

    //! @brief gives formatted analysis of the min and max statistics of the atom distance restraints
    //! @param MIN_MAX the min and max statistics of the atom distance restraints
    //! @param LINE_NAME_FORMAT the format to use to format the name of the line i.e. the first column
    //! @param FORMAT the format object to format the scores
    //! @return string which has the formatted min and max statistics analysis of atom restraint scores
    std::string AnalyzeAtomDistanceScore::GetMinMaxLines
    (
      const storage::Vector< math::RunningMinMax< double> > &MIN_MAX,
      const util::Format &LINE_NAME_FORMAT, const util::Format &FORMAT
    )
    {
      std::string min_line( LINE_NAME_FORMAT( "Min"));
      std::string max_line( LINE_NAME_FORMAT( "Max"));

      // iterate through the mean and std dev statistics
      for
      (
        storage::Vector< math::RunningMinMax< double> >::const_iterator
          stat_itr( MIN_MAX.Begin()), stat_itr_end( MIN_MAX.End());
        stat_itr != stat_itr_end; ++stat_itr
      )
      {
        min_line += FORMAT( stat_itr->GetMin());
        max_line += FORMAT( stat_itr->GetMax());
      }

      std::string analysis( min_line + "\n" + max_line);

      return analysis;
    }

  } // namespace restraint

} // namespace bcl
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
#include "restraint/bcl_restraint_analyze_atom_distance_score_heatmap.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "math/bcl_math_gnuplot_heatmap.h"
#include "restraint/bcl_restraint_epr_distance_data.h"
#include "score/bcl_score_restraint_distance_spin_label.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AnalyzeAtomDistanceScoreHeatmap::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new AnalyzeAtomDistanceScoreHeatmap())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AnalyzeAtomDistanceScoreHeatmap::AnalyzeAtomDistanceScoreHeatmap() :
      m_Score(),
      m_OutFilePostFix( ".AtomDistanceScoreHeatmap"),
      m_HistogramMinimum( -1.0),
      m_HistogramBinSize( 0.1),
      m_HistogramNumberOfBins( 10),
      m_RestraintType( EPRDistanceData::GetDefaultHandler())
    {
    }

    //! @brief Clone function
    //! @return pointer to new AnalyzeAtomDistanceScoreHeatmap
    AnalyzeAtomDistanceScoreHeatmap *AnalyzeAtomDistanceScoreHeatmap::Clone() const
    {
      return new AnalyzeAtomDistanceScoreHeatmap( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AnalyzeAtomDistanceScoreHeatmap::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &AnalyzeAtomDistanceScoreHeatmap::GetOutFilePostfix() const
    {
      return m_OutFilePostFix;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &AnalyzeAtomDistanceScoreHeatmap::GetAlias() const
    {
      static const std::string s_Name( "AtomDistanceScoreHeatmap");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string AnalyzeAtomDistanceScoreHeatmap::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // get the protein model data from one of the models of the ensemble and make sure it could be cast
      util::ShPtrVector< AtomDistance> data( m_RestraintType->ReadRestraintsFromFile());

      const storage::Vector< math::Histogram> score_histograms( GetScoreHistograms( data, ENSEMBLE));

      const storage::Vector< std::string> restraint_tics( GetRestraintNameTics( data));

      math::GnuplotHeatmap heatmap( GetHeatMap( score_histograms, restraint_tics));

      // write heat map to string stream
      std::stringstream stream;
      heatmap.WriteScript( stream);

      // return the string
      return stream.str();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AnalyzeAtomDistanceScoreHeatmap::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Score, ISTREAM);
      io::Serialize::Read( m_OutFilePostFix, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AnalyzeAtomDistanceScoreHeatmap::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Score, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_OutFilePostFix, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AnalyzeAtomDistanceScoreHeatmap::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "creates heat map showing how frequently a given score is achieved by each atom distance restraint"
      );

      parameters.AddInitializer
      (
        "restraint_score_type",
        "the atom distance restraint type where the score will come from",
        io::Serialization::GetAgent( &m_Score),
        score::RestraintDistanceSpinLabel::GetDefaultScheme()
      );

      parameters.AddInitializer
      (
        "restraint_type",
        "the restraint type; used to load the actual restraints",
        io::Serialization::GetAgent( &m_RestraintType),
        EPRDistanceData::GetDefaultHandler().GetString()
      );

      parameters.AddInitializer
      (
        "filename_postfix",
        "the postfix that will be appended to the output file in order to identify this analysis",
        io::Serialization::GetAgent( &m_OutFilePostFix),
        ".AnalyzeAtomDistanceScoreHeatmap"
      );

      parameters.AddInitializer
      (
        "histogram_minimum",
        "the minimal value representing the left boundary of the score histogram",
        io::Serialization::GetAgent( &m_HistogramMinimum),
        "-1.0"
      );

      parameters.AddInitializer
      (
        "histogram_binsize",
        "the width of one bin of the score histograms",
        io::Serialization::GetAgent( &m_HistogramBinSize),
        "0.1"
      );

      parameters.AddInitializer
      (
        "histogram_num_bins",
        "the number of bins in the score histograms",
        io::Serialization::GetAgent( &m_HistogramNumberOfBins),
        "10"
      );

      return parameters;
    }

    //! @brief creates histograms of how frequently a score is achieved by each restraint
    //! @param DATA the restraints that will be scored
    //! @param ENSEMBLE the ensemble of models that will be scored according to the atom distance restraints
    //! @return vector of histograms - one histogram for each restraint + one for average sum
    storage::Vector< math::Histogram> AnalyzeAtomDistanceScoreHeatmap::GetScoreHistograms
    (
      const util::ShPtrVector< AtomDistance> &DATA, const assemble::ProteinEnsemble &ENSEMBLE
    ) const
    {
      // takes into account adding columns for the average score
      static const size_t s_additional_columns( 1);

      // hold the score data for each restraint and score sum and total score
      storage::Vector< math::Histogram> score_histograms
      (
        DATA.GetSize() + s_additional_columns,
        math::Histogram( m_HistogramMinimum, m_HistogramBinSize, m_HistogramNumberOfBins)
      );

      // score to score the individual restraints
      const math::FunctionInterfaceSerializable< AtomDistanceAssignment, double> &individual_score( *m_Score->GetScore());

      const double data_size( DATA.GetSize());

      // iterate through the ensemble
      for
      (
        assemble::ProteinEnsemble::const_iterator protein_itr( ENSEMBLE.Begin()), protein_itr_end( ENSEMBLE.End());
        protein_itr != protein_itr_end;
        ++protein_itr
      )
      {
        const assemble::ProteinModel &model( **protein_itr);

        // iterators for the statistics vectors
        storage::Vector< math::Histogram>::iterator
          score_histograms_itr( score_histograms.Begin()), score_histograms_itr_end( score_histograms.End());

        // to hold sum of restraint scores for this protein
        double score_sum( 0);

        // iterate through the atom distance restraints
        for
        (
          util::ShPtrVector< AtomDistance>::const_iterator restraint_itr( DATA.Begin()), restraint_itr_end( DATA.End());
          restraint_itr != restraint_itr_end && score_histograms_itr != score_histograms_itr_end;
          ++restraint_itr, ++score_histograms_itr
        )
        {
          // get the assignment
          const AtomDistanceAssignment assignment( ( *restraint_itr)->GenerateAssignment( model));

          // get the current score
          const double score( individual_score( assignment));

          score_sum += score;

          // add the score to the two statistics objects
          score_histograms_itr->PushBack( score);
        } // end iterate through restraints

        // add the score sum and total score to the two histograms objects
        score_histograms_itr->PushBack( score_sum / data_size);
      } // end iterate through ensemble

      // normalize all the histograms
      for
      (
        storage::Vector< math::Histogram>::iterator
          histogram_itr( score_histograms.Begin()), histogram_itr_end( score_histograms.End());
        histogram_itr != histogram_itr_end;
        ++histogram_itr
      )
      {
        histogram_itr->Normalize();
        BCL_MessageDbg( "histogram " + util::Format()( histogram_itr_end - histogram_itr) + "\n" + util::Format()( *histogram_itr));
      }

      return score_histograms;
    }

    //! @brief creates the tics that will be used in the heatmap - converts each restraint into a string to use as tic
    //! @param DATA the restraints that will be used to create tics
    //! @return vector of strings which are the tics representing each restraint
    storage::Vector< std::string> AnalyzeAtomDistanceScoreHeatmap::GetRestraintNameTics
    (
      const util::ShPtrVector< AtomDistance> &DATA
    )
    {
      storage::Vector< std::string> tics;

      // iterate through the data to get the tics from the restraints
      for
      (
        util::ShPtrVector< AtomDistance>::const_iterator data_itr( DATA.Begin()), data_itr_end( DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        const assemble::LocatorAtomCoordinatesInterface &atom_locator_a( *( *data_itr)->GetData().First());
        const assemble::LocatorAtomCoordinatesInterface &atom_locator_b( *( *data_itr)->GetData().Second());
        tics.PushBack( assemble::LocatorAtomCoordinatesInterface::GetNameFromPair( atom_locator_a, atom_locator_b));
      }

      tics.PushBack( "mean");

      return tics;
    }

    //! @brief creates gnuplot heat map object from the score histograms and the restraint names as tics
    //! @param SCORE_HISTOGRAMS the histograms that will be used to make the heat map
    //! @param TICS the names of the restraints
    //! @return heat map object which represents the distribution of scores for each restraint
    math::GnuplotHeatmap AnalyzeAtomDistanceScoreHeatmap::GetHeatMap
    (
      const storage::Vector< math::Histogram> &SCORE_HISTOGRAMS, const storage::Vector< std::string> &TICS
    )
    {
      const util::SiPtrVector< const math::Histogram> histograms
      (
        util::ConvertToConstSiPtrVector( SCORE_HISTOGRAMS)
      );

      math::GnuplotHeatmap heatmap;

      heatmap.SetFromHistograms( histograms, true, true);
      heatmap.SetTitleAndLabel( "Title", "restraint", "score", "Frequency (fraction of models)");
      heatmap.SetTicsX( TICS, true, 1);
      heatmap.SetRotationXTics( 90);
      heatmap.SetFilename( "AnalyzeAtomDistanceScoreHeatmap.gnuplot");
      heatmap.SetFont( "/usr/share/fonts/dejavu-lgc/DejaVuLGCSansMono.ttf", 10);

      return heatmap;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_analyze_coordinate_distance_distribution.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_gnuplot_heatmap.h"
#include "math/bcl_math_gnuplot_multiplot.h"
#include "math/bcl_math_histogram.h"
#include "math/bcl_math_running_average_sd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AnalyzeCoordinateDistanceDistribution::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new AnalyzeCoordinateDistanceDistribution())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AnalyzeCoordinateDistanceDistribution::AnalyzeCoordinateDistanceDistribution() :
      m_OutFilePostFix( ".CoordinateDistanceDistribution"),
      m_CoordinateA( assemble::LocatorAtom( 'A', 1, biol::GetAtomTypes().CA)),
      m_CoordinateB( assemble::LocatorAtom( 'A', 1, biol::GetAtomTypes().CA)),
      m_ComparisonFunction(),
      m_HistogramMinimum( 0),
      m_HistogramBinSize( 2),
      m_HistogramNumberOfBins( 30),
      m_Title( "Title"),
      m_PixelX   ( 600),
      m_PixelY   ( 400),
      m_Font     ( "Arial"),
      m_FontSize ( 12),
      m_GreyScale( false),
      m_PlotRatio( util::GetUndefinedDouble()),
      m_MeanStdDevOutFile( "mean_stddev.txt"),
      m_BinsPerTic( 1),
      m_CenterTics( false),
      m_MinZ( util::GetUndefinedDouble()),
      m_MaxZ( util::GetUndefinedDouble()),
      m_Normalize( true)
      {
      }

    //! @brief Clone function
    //! @return pointer to new AnalyzeCoordinateDistanceDistribution
    AnalyzeCoordinateDistanceDistribution *AnalyzeCoordinateDistanceDistribution::Clone() const
    {
      return new AnalyzeCoordinateDistanceDistribution( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AnalyzeCoordinateDistanceDistribution::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &AnalyzeCoordinateDistanceDistribution::GetOutFilePostfix() const
    {
      return m_OutFilePostFix;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &AnalyzeCoordinateDistanceDistribution::GetAlias() const
    {
      static const std::string s_Name( "CoordinateDistanceDistribution");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string AnalyzeCoordinateDistanceDistribution::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // intialize
      math::GnuplotMultiplot multiplot;
      multiplot.SetRowsCols( 2, 1);
      multiplot.SetFont( m_Font, m_FontSize);
      multiplot.SetPixelAndRatio( m_PixelX, m_PixelY, util::GetUndefinedDouble());

      // ensemble histogram plot
      multiplot.Insert( GetHeatMap( GetDistanceHistogram( ENSEMBLE)));

      // comparison function histogram plot
      multiplot.Insert( GetReferenceHeatMap());

      // write heat map to string stream
      std::stringstream stream;
      multiplot.WriteScript( stream);
      std::string analysis( stream.str());

      // return analysis string of the gnuplot script
      return analysis;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AnalyzeCoordinateDistanceDistribution::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Creates heat map of showing frequency with which a distance between two coordinates is observed. The heat "
        "map of frequency with which a distance is observed for two coordinates located in a protein ensemble can "
        "be compared with a second heat map. The second heat map is given by a math function. This allows the "
        "distribution of distances coming from the ensemble to be compared with a second distribution which is "
        "perhaps the expected distribution. An example use of this is to see the distribution of distances coming "
        "from two residues in the ensemble whose distance has been measured by EPR. This ensemble distribution of "
        "distances can then be compared to the EPR distribution described by a gaussian function defined by the "
        "mean and standard deviation of the distance probability distribution measured by EPR."
      );

      parameters.AddInitializer
      (
        "filename_postfix",
        "the postfix that will be appended to the output file in order to identify this analysis",
        io::Serialization::GetAgent( &m_OutFilePostFix),
        ".ChiAnglePairDistribution"
      );

      parameters.AddInitializer
      (
        "coord_a",
        "the first coordinate that will be used in the distance",
        io::Serialization::GetAgent( &m_CoordinateA)
      );

      parameters.AddInitializer
      (
        "coord_b",
        "the second coordinate that will be used in the distance",
        io::Serialization::GetAgent( &m_CoordinateB)
      );

      parameters.AddInitializer
      (
        "function",
        "the function that the coordinate heat map will be compared against",
        io::Serialization::GetAgent( &m_ComparisonFunction)
      );

      parameters.AddInitializer
      (
        "histogram_minimum",
        "the minimal value representing the left boundary of the score histogram",
        io::Serialization::GetAgent( &m_HistogramMinimum),
        "-1.0"
      );

      parameters.AddInitializer
      (
        "histogram_binsize",
        "the width of one bin of the score histograms",
        io::Serialization::GetAgent( &m_HistogramBinSize),
        "0.1"
      );

      parameters.AddInitializer
      (
        "histogram_num_bins",
        "the number of bins in the score histograms",
        io::Serialization::GetAgent( &m_HistogramNumberOfBins),
        "10"
      );

      parameters.AddInitializer
      (
        "title",
        "the title that will label the resulting heat map",
        io::Serialization::GetAgent( &m_Title),
        "Title"
      );

      parameters.AddInitializer
      (
        "pixel_x",
        "The size of the plot png in the x direction",
        io::Serialization::GetAgent( &m_PixelX),
        "600"
      );

      parameters.AddInitializer
      (
        "pixel_y",
        "The size of the plot png in the y direction",
        io::Serialization::GetAgent( &m_PixelY),
        "400"
      );

      parameters.AddInitializer
      (
        "font",
        "The font to be used for text in the plot",
        io::Serialization::GetAgent( &m_Font),
        "Arial"
      );

      parameters.AddInitializer
      (
        "font_size",
        "The size of the font for the plot",
        io::Serialization::GetAgent( &m_FontSize),
        "12"
      );

      parameters.AddInitializer
      (
        "grey_scale",
        "boolean true the color palette for the gradient should be grey scale. 1=true;0=false",
        io::Serialization::GetAgent( &m_GreyScale),
        "0"
      );

      parameters.AddInitializer
      (
        "plot_ratio",
        "How much of the png area the plot will cover top to bottom and left to write."
        " (this does not include any labels or tics). So if you say 0.5, plot will go from 0 to 0.5 from bottom"
        " to top and left to right",
        io::Serialization::GetAgent( &m_PlotRatio),
        util::Format()( util::GetUndefinedDouble())
      );

      parameters.AddInitializer
      (
        "mean_stddev_out_file",
        "Full path and filename for holding the mean and standard deviation of the distribution",
        io::Serialization::GetAgent( &m_MeanStdDevOutFile),
        "mean_stddev.txt"
      );

      parameters.AddInitializer
      (
        "bins_per_tic",
        "how many bins per tic in the plot. 1 means every bin is labeled, 2 means every other, etc.",
        io::Serialization::GetAgent( &m_BinsPerTic),
        "1"
      );

      parameters.AddInitializer
      (
        "center_tics",
        "boolean if true, tics will be centered on the bins, if false, will be on edges of bins. 1=true;0=false",
        io::Serialization::GetAgent( &m_CenterTics),
        "0"
      );

      parameters.AddInitializer
      (
        "min_z",
        "the minimum value used in pymol for the z axis (i.e. color). Values below this will just get the minimum value color.",
        io::Serialization::GetAgent( &m_MinZ),
        util::Format()( util::GetUndefinedDouble())
      );
      parameters.AddInitializer
      (
        "max_z",
        "the maximum value used in pymol for the z axis (i.e. color). Values above this will just get the maximum value color.",
        io::Serialization::GetAgent( &m_MaxZ),
        util::Format()( util::GetUndefinedDouble())
      );

      parameters.AddInitializer
      (
        "normalize_histogram",
        "boolean true indicates that the distance distribution histogram will be normalized. 1=true;0=false",
        io::Serialization::GetAgent( &m_Normalize),
        "1"
      );

      return parameters;
    }

    //! @brief gives the histogram of the distance distribution of the protein ensemble
    //! @param ENSEMBLE the ensemble for which the distance distribution will be calculated
    //! @return histogram that has the distance distribution calculated from the provided protein ensemble
    math::Histogram AnalyzeCoordinateDistanceDistribution::GetDistanceHistogram
    (
      const assemble::ProteinEnsemble &ENSEMBLE
    ) const
    {
      // initialize the histogram distance distribution
      math::Histogram distance_distribution( m_HistogramMinimum, m_HistogramBinSize, m_HistogramNumberOfBins);

      // to keep track of the mean and standard deviation of the distances in the ensemble
      math::RunningAverageSD< double> mean_stddev;

      // iterate through the ensemble to get the distribution of distances
      for
      (
        assemble::ProteinEnsemble::const_iterator
          ensemble_itr( ENSEMBLE.Begin()), ensemble_itr_end( ENSEMBLE.End());
        ensemble_itr != ensemble_itr_end;
        ++ensemble_itr
      )
      {
        // locate the two coordinates and make sure they are defined
        const linal::Vector3D coord_a( m_CoordinateA->Locate( **ensemble_itr));
        BCL_Assert( coord_a.IsDefined(), "coord_a could not be located with " + util::Format()( *m_CoordinateA));
        const linal::Vector3D coord_b( m_CoordinateB->Locate( **ensemble_itr));
        BCL_Assert( coord_b.IsDefined(), "coord_b could not be located with " + util::Format()( *m_CoordinateB));

        // calculate the distance between the two coordinates and make sure it is defined
        const double distance( linal::Distance( coord_a, coord_b));
        BCL_Assert( util::IsDefined( distance), "distance is not defined");

        // add the distance into the distance distribution and the statistics object
        distance_distribution.PushBack( distance);
        mean_stddev += distance;
      }

      // write the mean and standard deviation of the observed distances to the desired output file
      io::OFStream write;
      io::File::MustOpenOFStream( write, m_MeanStdDevOutFile);
      write << m_Title << "\tmean " << mean_stddev.GetAverage() << "\tstddev " << mean_stddev.GetStandardDeviation();
      io::File::CloseClearFStream( write);

      // true if desired to normalize
      if( m_Normalize)
      {
        // normalize distribution
        distance_distribution.Normalize();
      }

      // return the histogram of the distances
      return distance_distribution;
    }

    //! @brief gives the heat map corresponding to the provided histogram made from the protein ensemble
    //! @param HISTGRAM the histogram for which the heat map will be made
    //! @return heat map representing the distance distribution of the protein ensemble
    util::ShPtr< math::GnuplotHeatmap>
    AnalyzeCoordinateDistanceDistribution::GetHeatMap( const math::Histogram &HISTOGRAM) const
    {
      // make new heatmap
      util::ShPtr< math::GnuplotHeatmap> gnuplot( new math::GnuplotHeatmap());

      // set the heatmap from the provided histogram
      gnuplot->SetFromHistogram( HISTOGRAM, false, false);

      // set options on the heat map
      gnuplot->SetShowColorBox( false);
      gnuplot->SetTitleAndLabel( "", "", "", "");
      gnuplot->SetMargins( 0.5, 0.5 - m_PlotRatio, 0.99, 0.01);
      gnuplot->SetWritePreHeader( false);
      gnuplot->SetWriteHeader( true);
      gnuplot->SetPalette( math::GnuplotHeatmap::e_GreyScale);
      gnuplot->SetRotationXTics( 90);
      gnuplot->SetNoMirrorTics( true);

      // set the tics of the heat map according to the preferences of the user
      {
        const bool center_tic_x( m_CenterTics);
        const linal::Vector< double> binning_x
        (
          linal::FillVector< double>
          (
            HISTOGRAM.GetNumberOfBins() + size_t( !center_tic_x),
            HISTOGRAM.GetBoundaries().First() + ( center_tic_x ? 0.5 * HISTOGRAM.GetBinSize() : 0.0),
            HISTOGRAM.GetBinSize()
          )
        );
        BCL_Assert
        (
          gnuplot->SetTicsX
          (
            math::GnuplotHeatmap::TicsFromBinning( binning_x, m_BinsPerTic, util::Format().W( 4)),
            center_tic_x,
            m_BinsPerTic
          ), "unable to set tics x"
        );
      }

      // set the min max z values
      gnuplot->SetMinMaxZ( m_MinZ, m_MaxZ);

      // return the heat map
      return gnuplot;
    }

    //! @brief gives the heat map corresponding to the reference function
    //! @return shptr to a heat map that represents the reference function
    util::ShPtr< math::GnuplotHeatmap> AnalyzeCoordinateDistanceDistribution::GetReferenceHeatMap() const
    {
      // make new heatmap
      util::ShPtr< math::GnuplotHeatmap> gnuplot( new math::GnuplotHeatmap());

      // set the heat map from the reference function
      gnuplot->SetFromFunction
      (
        *m_ComparisonFunction, m_HistogramNumberOfBins,
        m_HistogramMinimum + ( m_HistogramBinSize / 2.0), m_HistogramBinSize, false, false, m_Normalize
      );

      // set options on the heat map
      gnuplot->SetMargins( 0.5 + m_PlotRatio, 0.5, 0.99, 0.01);
      gnuplot->SetNoMirrorTics( true);
      gnuplot->SetWritePreHeader( false);
      gnuplot->SetWriteHeader( true);
      gnuplot->SetTitleAndLabel( m_Title, "", "", "");
      gnuplot->SetPalette( math::GnuplotHeatmap::e_GreyScale);
      gnuplot->SetShowXTics( false);
      // set the min max z values
      gnuplot->SetMinMaxZ( m_MinZ, m_MaxZ);

      // return the heat map
      return gnuplot;
    }
  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_analyze_per_residue_rmsd_between_ensembles.h"

// includes from bcl - sorted alphabetically
#include "align/bcl_align_alignment_interface.h"
#include "assemble/bcl_assemble_locator_atom.h"
#include "assemble/bcl_assemble_quality.h"
#include "command/bcl_command_flag_interface.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serialization_via_static_functions.h"
#include "quality/bcl_quality_rmsd.h"
#include "restraint/bcl_restraint_analyze_accessibility_change.h"
#include "restraint/bcl_restraint_analyze_per_residue_rmsd.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AnalyzePerResidueRMSDBetweenEnsembles::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new AnalyzePerResidueRMSDBetweenEnsembles())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AnalyzePerResidueRMSDBetweenEnsembles::AnalyzePerResidueRMSDBetweenEnsembles() :
      m_OutFilePostFix( ".AnalyzePerResidueRMSDBetweenEnsembles"),
      m_SuperimposeMeasure(),
      m_TemplateModel(),
      m_QualityMeasure( quality::GetMeasures().e_RMSD),
      m_PymolOuputFilename( "AnalyzePerResidueRMSDBetweenEnsembles.py"),
      m_StartEnsemble(),
      m_NormalizeByInternalRMSD( false)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AnalyzePerResidueRMSDBetweenEnsembles
    AnalyzePerResidueRMSDBetweenEnsembles *AnalyzePerResidueRMSDBetweenEnsembles::Clone() const
    {
      return new AnalyzePerResidueRMSDBetweenEnsembles( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AnalyzePerResidueRMSDBetweenEnsembles::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &AnalyzePerResidueRMSDBetweenEnsembles::GetOutFilePostfix() const
    {
      return m_OutFilePostFix;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &AnalyzePerResidueRMSDBetweenEnsembles::GetAlias() const
    {
      static const std::string s_name( "PerResidueRMSDBetweenEnsembles");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string AnalyzePerResidueRMSDBetweenEnsembles::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // string to hold analysis
      std::string analysis;

      // open write to pymol output script
      io::OFStream write;
      io::File::MustOpenOFStream
      (
        write, GetFlagOutFilePrefix()->GetFirstParameter()->GetValue() + m_PymolOuputFilename
      );

      // write python script loading of template model
      write << "cmd.load( \""
            <<  util::ShPtr< util::Wrapper< std::string> >
                (
                  m_TemplateModel.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile)
                )->GetData()
            << "\", \"template_model\")\n";

      // map to hold for each residue the corresponding coordinates
      storage::Map< assemble::LocatorAtom, math::RunningAverageSD< double> > atom_coordinate_list;

      // for holding atom lists if want to normalize by the internal per residue rmsd of ENSEMBLE
      storage::Map< assemble::LocatorAtom, storage::Vector< linal::Vector3D> > atom_coordinate_list_norm;
      if( m_NormalizeByInternalRMSD)
      {
        // map to hold for each residue the corresponding coordinates
        atom_coordinate_list_norm = AnalyzePerResidueRMSD::GetAtomCoordinates( ENSEMBLE, m_SuperimposeMeasure, m_TemplateModel);
      }

      // iterate through the current ensemble
      for
      (
        assemble::ProteinEnsemble::const_iterator ensemble_itr( ENSEMBLE.Begin()), ensemble_itr_end( ENSEMBLE.End());
        ensemble_itr != ensemble_itr_end;
        ++ensemble_itr
      )
      {
        // copy the protein model
        assemble::ProteinModel model( **ensemble_itr);

        // atom type for superimposing
        static const biol::AtomType s_atom_type( biol::GetAtomTypes().CA);

        // atom types for superimposing
        static const storage::Set< biol::AtomType> s_atom_types( s_atom_type);

        // iterate over the start ensemble
        for
        (
          assemble::ProteinEnsemble::const_iterator
            start_ensemble_itr( m_StartEnsemble.Begin()), start_ensemble_itr_end( m_StartEnsemble.End());
            start_ensemble_itr != start_ensemble_itr_end;
          ++start_ensemble_itr
        )
        {
          // superimpose the current model onto the start structure
          assemble::Quality::SuperimposeModel( m_SuperimposeMeasure, model, **start_ensemble_itr, s_atom_types);

          // get alignments for each chain
          storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > alignments
          (
            assemble::Quality::CreateAlignmentProteinModels( model, **start_ensemble_itr)
          );

          // iterate over the chain alignments
          for
          (
            storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >::const_iterator
              chain_itr( alignments.Begin()), chain_itr_end( alignments.End());
            chain_itr != chain_itr_end;
            ++chain_itr
          )
          {
            util::ShPtrList< align::Assignment< biol::AABase> > assignments( chain_itr->second->GetAssignments());

            // iterate over the assignments
            for
            (
              util::ShPtrList< align::Assignment< biol::AABase> >::const_iterator
                assignment_itr( assignments.Begin()), assignment_itr_end( assignments.End());
              assignment_itr != assignment_itr_end;
              ++assignment_itr
            )
            {
              util::SiPtrList< const biol::AABase> members( ( *assignment_itr)->GetMembers());

              // to hold the coords in this assignment
              util::SiPtrVector< const linal::Vector3D> coords;

              // create locator atom for current residue
              assemble::LocatorAtom locator;

              // iterate over the members to get the coordinates
              for
              (
                util::SiPtrList< const biol::AABase>::const_iterator
                  member_itr( members.Begin()), member_itr_end( members.End());
                member_itr != member_itr_end;
                ++member_itr
              )
              {
                const util::SiPtr< const biol::AABase> &aa_base( *member_itr);

                // if the residue is defined
                if( aa_base.IsDefined())
                {
                  // add coordinates to list of coordinates
                  util::SiPtrVector< const linal::Vector3D> current_coord( aa_base->GetAtomCoordinates( s_atom_types));

                  // make sure one coordinate was gotten
                  BCL_Assert( current_coord.GetSize() == 1, "cord size not 1");
                  const linal::Vector3D &ca_coord( *current_coord.FirstElement());

                  // true if the coordinates of the atom are defined
                  if( ca_coord.IsDefined())
                  {
                    // add coords and set locator (locator will be set with each iteration but will be the same)
                    coords.PushBack( ca_coord);
                    locator = assemble::LocatorAtom( aa_base->GetChainID(), aa_base->GetSeqID(), s_atom_type);
                  }
                }
              }

              // should be two coordinates in the list, one for atom from start model one from end model
              if( coords.GetSize() == 2)
              {
                // get the rmsd
                const storage::Pair< storage::VectorND< 2, double>, math::RunningAverageSD< double> > rmsd_info
                (
                  quality::RMSD::RealSpaceRMSDPairwise( coords)
                );

                // for subtracting the background that is the variation within the ENSEMBLE for this residue
                double normalization_factor( 0);
                if( m_NormalizeByInternalRMSD)
                {
                  // try to find the current residue
                  storage::Map< assemble::LocatorAtom, storage::Vector< linal::Vector3D> >::const_iterator norm_itr
                  (
                    atom_coordinate_list_norm.Find( locator)
                  );

                  // true if the residue was found
                  if( norm_itr != atom_coordinate_list_norm.End())
                  {
                    // calculate the rmsd of the coordinates to one another within ENSEMBLE
                    storage::Vector< linal::Vector3D> coords( norm_itr->second);
                    const util::SiPtrVector< const linal::Vector3D> all_coords( util::ConvertToSiPtrVector( coords));
                    const storage::Pair< storage::VectorND< 2, double>, math::RunningAverageSD< double> > rmsd_info_norm
                    (
                      quality::RMSD::RealSpaceRMSDPairwise( all_coords)
                    );
                    // rmsd is the background normalization factor
                    normalization_factor = rmsd_info_norm.First().First();
                  }
                }

                // current rmsd between the ensembles is the actual rmsd between the two ensembles minus the
                // normalization factor
                const double current_rmsd( rmsd_info.First().First() - normalization_factor);

                BCL_MessageDbg
                (
                  "considering rmsd " + util::Format()( current_rmsd) +
                   " for residue " + locator.GetIdentification()
                );

                // add the rmsd to the map
                atom_coordinate_list[ locator] += current_rmsd;
              }
              else //< could be that coordinate is not defined for either of the models, then can't calculate rmsd
              {
                continue;
              }
            } // assignments
          } // alignment chain map
        } // ensemble
      } // ensemble

      // iterate over the map of rmsds for each residue
      for
      (
        storage::Map< assemble::LocatorAtom, math::RunningAverageSD< double> >::const_iterator
          map_itr( atom_coordinate_list.Begin()), map_itr_end( atom_coordinate_list.End());
        map_itr != map_itr_end;
        ++map_itr
      )
      {
        // current atom
        const assemble::LocatorAtom &atom( map_itr->first);

        // write the analysis information
        analysis += atom.GetIdentification() + " rmsd " + util::Format().W( 9).FFP( 3)( map_itr->second.GetAverage()) +
                    " rmsd_stddev " + util::Format().W( 9).FFP( 3)( map_itr->second.GetStandardDeviation()) + '\n';

        // write to the python script for visualization
        write << "cmd.alter( \"template_model and " << atom.GetPymolResidueSelection() << "\", \"b= "
              << util::Format()( map_itr->second.GetAverage()) << "\")\n";
      }

      // write the ending information to the python script for visualization
      write <<  "cmd.show_as(\"cartoon\"   ,\"template_model\")\n";
      write <<  "cmd.cartoon(\"putty\"   ,\"template_model\")\n";
      write <<  "cmd.unset(\"cartoon_smooth_loops\"   ,\"template_model\")\n";
      write <<  "cmd.unset(\"cartoon_flat_sheets\"   ,\"template_model\")\n";

      // return the analysis string
      return analysis;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AnalyzePerResidueRMSDBetweenEnsembles::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates the per residue RMSD between two sets of ensembles. CA atoms are used in the calculation. Each "
        " model in each ensemble is pairwise superimposed according to the provided superimposition method, then the "
        " per residue rmsd is calculated along the sequence. These rmsds are averaged over all pairwise RMSDs between "
        " the ensembles. Also outputs a python file that can be opened with pymol in order to vizualize the regions "
        " of the provided template structure that have large per residue RMSDs. This is depicted as thickness using "
        " the putty representation in pymol."
      );

      parameters.AddInitializer
      (
        "filename_postfix",
        "the postfix that will be appended to the output file in order to identify this analysis",
        io::Serialization::GetAgent( &m_OutFilePostFix),
        ".AnalyzeAtomDistanceScore"
      );

      parameters.AddInitializer
      (
        "superimpose_method",
        "the method that should be used for superimposing onto the template",
        io::Serialization::GetAgent( &m_SuperimposeMeasure),
        quality::GetSuperimposeMeasures().e_RMSD.GetName()
      );

      parameters.AddInitializer
      (
        "template_pdb_filename",
        "the pdb filename of the protein model the other models will be superimposed onto",
        util::OwnPtr< io::SerializationInterface>
        (
          new io::SerializationViaStaticFunctions< assemble::ProteinModel>
          (
            command::ParameterCheckFileExistence(),
            ( &AnalyzePerResidueRMSD::ProteinModelAsString),
            ( &AnalyzePerResidueRMSD::ProteinModelFromString),
            &m_TemplateModel
          )
        ),
        "template.pdb"
      );

      parameters.AddInitializer
      (
        "pymol_output_filename",
        "The filename of the pymol script that will be outputted showing the accessibility changes",
        io::Serialization::GetAgent( &m_PymolOuputFilename),
        "accessibilities.pml"
      );

      parameters.AddInitializer
      (
        "start_ensemble_filename",
        "the name of the file containing the list of pdbs for the starting ensemble",
        util::OwnPtr< io::SerializationInterface>
        (
          new io::SerializationViaStaticFunctions< assemble::ProteinEnsemble>
          (
            command::ParameterCheckFileExistence(),
            ( &AnalyzeAccessibilityChange::EnsembleAsFilename),
            ( &AnalyzeAccessibilityChange::EnsembleFromFilename),
            &m_StartEnsemble
          )
        ),
        "start_ensemble_pdbs.ls"
      );

      parameters.AddInitializer
      (
        "normalize_by_background",
        "The background is the per residue rmsd within the ending ensemble i.e. the one provided"
        " by the main ensemble flag. If this flag is set (referring to this flag, not the main ensemble flag)"
        " then the per residue rmsd of the ending ensemble will be calculated after superimposing each model"
        " onto the provided template pdb. This per residue rmsd will then be subtracted from the per residue RMSD"
        " calculated between the two provided ensembles. This is helpful for seeing differences between the "
        " ensembles that are larger than the difference within the ending ensemble. 1=true;0=false",
        io::Serialization::GetAgent( &m_NormalizeByInternalRMSD),
        "0"
      );

      return parameters;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_analyze_per_residue_rmsd.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_atom.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_quality.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serialization_via_static_functions.h"
#include "pdb/bcl_pdb_factory.h"
#include "quality/bcl_quality_rmsd.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AnalyzePerResidueRMSD::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new AnalyzePerResidueRMSD())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AnalyzePerResidueRMSD::AnalyzePerResidueRMSD() :
      m_OutFilePostFix( ".AnalyzePerResidueRMSD"),
      m_SuperimposeMeasure(),
      m_TemplateModel(),
      m_QualityMeasure( quality::GetMeasures().e_RMSD),
      m_PymolOuputFilename( "AnalyzePerResidueRMSD.py")
    {
    }

    //! @brief Clone function
    //! @return pointer to new AnalyzePerResidueRMSD
    AnalyzePerResidueRMSD *AnalyzePerResidueRMSD::Clone() const
    {
      return new AnalyzePerResidueRMSD( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AnalyzePerResidueRMSD::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &AnalyzePerResidueRMSD::GetOutFilePostfix() const
    {
      return m_OutFilePostFix;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &AnalyzePerResidueRMSD::GetAlias() const
    {
      static const std::string s_name( "PerResidueRMSD");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string AnalyzePerResidueRMSD::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // string to hold analysis
      std::string analysis;

      // open write to pymol output script
      io::OFStream write;
      io::File::MustOpenOFStream
      (
        write, GetFlagOutFilePrefix()->GetFirstParameter()->GetValue() + m_PymolOuputFilename
      );

      // write python script loading of template model
      write << "cmd.load( \""
            <<  util::ShPtr< util::Wrapper< std::string> >
                (
                  m_TemplateModel.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile)
                )->GetData()
            << "\", \"template_model\")\n";

      // map to hold for each residue the corresponding coordinates
      storage::Map< assemble::LocatorAtom, storage::Vector< linal::Vector3D> > atom_coordinate_list
      (
        GetAtomCoordinates( ENSEMBLE, m_SuperimposeMeasure, m_TemplateModel)
      );

      // iterate over the map of coordinates for each residue
      for
      (
        storage::Map< assemble::LocatorAtom, storage::Vector< linal::Vector3D> >::const_iterator
          map_itr( atom_coordinate_list.Begin()), map_itr_end( atom_coordinate_list.End());
        map_itr != map_itr_end;
        ++map_itr
      )
      {
        // current atom
        const assemble::LocatorAtom &atom( map_itr->first);

        // vector of coordinates for the current atom over the whole ensemble
        storage::Vector< linal::Vector3D> coords( map_itr->second);

        BCL_MessageDbg( " atom " + atom.GetIdentification());

        // calculate the rmsd of the coordinates to one another
        const util::SiPtrVector< const linal::Vector3D> all_coords( util::ConvertToSiPtrVector( coords));
        const storage::Pair< storage::VectorND< 2, double>, math::RunningAverageSD< double> > rmsd_info
        (
          quality::RMSD::RealSpaceRMSDPairwise( all_coords)
        );

        // if rmsd is not defined
        if( !util::IsDefined( rmsd_info.First().First()))
        {
          // skip to next
          continue;
        }

        // write the analysis information
        analysis += atom.GetIdentification() + " rmsd " + util::Format().W( 9).FFP( 3)( rmsd_info.First().First()) +
                    " rmsd_stddev " + util::Format().W( 9).FFP( 3)( rmsd_info.First().Second()) +
                    " mean_dist " + util::Format().W( 9).FFP( 3)( rmsd_info.Second().GetAverage()) +
                    " dist_stddev " + util::Format().W( 9).FFP( 3)( rmsd_info.Second().GetStandardDeviation()) + '\n';

        // write to the python script for visualization
        write << "cmd.alter( \"template_model and " << atom.GetPymolResidueSelection() << "\", \"b= "
              << util::Format()( rmsd_info.First().First()) << "\")\n";

      }

      // write the ending information to the python script for visualization
      write <<  "cmd.show_as(\"cartoon\"   ,\"template_model\")\n";
      write <<  "cmd.cartoon(\"putty\"   ,\"template_model\")\n";
      write <<  "cmd.unset(\"cartoon_smooth_loops\"   ,\"template_model\")\n";
      write <<  "cmd.unset(\"cartoon_flat_sheets\"   ,\"template_model\")\n";

      return analysis;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AnalyzePerResidueRMSD::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "For an ensemble of models, superimposes them onto a template structure using the provided superimpose method. "
        "The rmsd of each residue between the ensemble of models is then calculated. Also outputs a python script that"
        " can be used with pymol to visualize the rmsd over the structure by using the putty representation. Residues"
        " with larger rmsd will be thicker. The template model is used in the pymol script."
      );

      parameters.AddInitializer
      (
        "filename_postfix",
        "the postfix that will be appended to the output file in order to identify this analysis",
        io::Serialization::GetAgent( &m_OutFilePostFix),
        ".AnalyzeAtomDistanceScore"
      );

      parameters.AddInitializer
      (
        "superimpose_method",
        "the method that should be used for superimposing onto the template",
        io::Serialization::GetAgent( &m_SuperimposeMeasure),
        quality::GetSuperimposeMeasures().e_RMSD.GetName()
      );

      parameters.AddInitializer
      (
        "template_pdb_filename",
        "the pdb filename of the protein model the other models will be superimposed onto",
        util::OwnPtr< io::SerializationInterface>
        (
          new io::SerializationViaStaticFunctions< assemble::ProteinModel>
          (
            command::ParameterCheckFileExistence(),
            ( &ProteinModelAsString),
            ( &ProteinModelFromString),
            &m_TemplateModel
          )
        ),
        "template.pdb"
      );

      parameters.AddInitializer
      (
        "pymol_output_filename",
        "The filename of the pymol script that will be outputted showing the accessibility changes",
        io::Serialization::GetAgent( &m_PymolOuputFilename),
        "accessibilities.pml"
      );

      return parameters;
    }

    //! @brief converts a protein model into a string by returning the pdb filename stored in the protein model data
    //! @param MODEL the model that will be converted into a string
    //! @return the string of the pdb filename stored in the protein model data
    std::string AnalyzePerResidueRMSD::ProteinModelAsString( const assemble::ProteinModel &MODEL)
    {
      util::ShPtr< util::Wrapper< std::string> > pdb
      (
        MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile)
      );

      if( !pdb.IsDefined())
      {
        return "none.pdb";
      }

      return *pdb;
    }

    //! @brief converts a string into a protein model by assuming the string is a pdb filename
    //! @param MODEL protein model to setup
    //! @param NAME string which is the pdb filename
    //! @param ERR_STREAM stream to write out erros to
    bool AnalyzePerResidueRMSD::ProteinModelFromString
    (
      assemble::ProteinModel &MODEL, const std::string &NAME, std::ostream &ERR_STREAM
    )
    {
      // instantiate the pdb factory
      pdb::Factory factory;

      // instantiate protein model
      MODEL = factory.ProteinModelFromPDBFilename( NAME);

      // set the pdb filename in the protein model data
      util::ShPtr< assemble::ProteinModelData>( MODEL.GetProteinModelData())->Insert
      (
        assemble::ProteinModelData::e_PDBFile,
        util::ShPtr< util::Wrapper< std::string> >( new util::Wrapper< std::string>( NAME))
      );
      return true;
    }

    //! @brief gives the CA atom coordinates for each residue for each model in an ensemble
    //!        Each model in the ensemble is superimposed onto a templete model according to the provided
    //!        superimposition method.
    //! @param ENSEMBLE the models whose ca coordinate will be gotten for each residue
    //! @param SUPERIMPOSE_MEASURE the method of superimposition used for the models of the ensemble onto the template
    //! @param TEMPLATE_MODEL the model that the ensemble models will be superimposed on
    //! @return map with a locator for each atom and the vector of associated coordinates coming from the models in
    //!         the ensemble
    storage::Map< assemble::LocatorAtom, storage::Vector< linal::Vector3D> > AnalyzePerResidueRMSD::GetAtomCoordinates
    (
      const assemble::ProteinEnsemble &ENSEMBLE,
      const quality::SuperimposeMeasure &SUPERIMPOSE_MEASURE,
      const assemble::ProteinModel &TEMPLATE_MODEL
    )
    {

      // map to hold for each residue the corresponding coordinates
      storage::Map< assemble::LocatorAtom, storage::Vector< linal::Vector3D> > atom_coordinate_list;

      // iterate through the ensemble
      for
      (
        assemble::ProteinEnsemble::const_iterator ensemble_itr( ENSEMBLE.Begin()), ensemble_itr_end( ENSEMBLE.End());
        ensemble_itr != ensemble_itr_end;
        ++ensemble_itr
      )
      {
        // copy the protein model
        assemble::ProteinModel model( **ensemble_itr);

        // atom type for superimposing
        static const biol::AtomType s_atom_type( biol::GetAtomTypes().CA);

        // atom types for superimposing
        static const storage::Set< biol::AtomType> s_atom_types( s_atom_type);

        // superimpose the current model onto the template
        assemble::Quality::SuperimposeModel( SUPERIMPOSE_MEASURE, model, TEMPLATE_MODEL, s_atom_types);

        // get the residues from the superimposed model
        const util::SiPtrVector< const biol::AABase> aas( model.GetAminoAcids());

        // iterate through the residues of the current model
        for
        (
          util::SiPtrVector< const biol::AABase>::const_iterator aa_itr( aas.Begin()), aa_itr_end( aas.End());
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {
          // current residue
          const biol::AABase &current_aa( **aa_itr);

          // create locator atom for current residue
          const assemble::LocatorAtom locator( current_aa.GetChainID(), current_aa.GetSeqID(), s_atom_type);

          // get the coordinates for the atom types
          const util::SiPtrVector< const linal::Vector3D> coords( current_aa.GetAtomCoordinates( s_atom_types));

          // assert there is only one set of coordinates, otherwise calculating rmsd is worse
          BCL_Assert
          (
            coords.GetSize() == 1, "number of coordinates is " + util::Format()( coords.GetSize()) +
            " but probably should be 1"
          );

          if( !coords.FirstElement()->IsDefined())
          {
            continue;
          }

          // add coordinates for the atom into the map
          atom_coordinate_list[ locator].Append
          (
            util::ConvertToStorageVector< linal::Vector3D>( current_aa.GetAtomCoordinates( s_atom_types))
          );
        }
      }

      return atom_coordinate_list;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_analyze_sas.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "biol/bcl_biol_sasa_data.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_histogram.h"
#include "restraint/bcl_restraint_sas_analysis.h"
#include "restraint/bcl_restraint_sas_debye_interface.h"
#include "restraint/bcl_restraint_sas_optimization.h"
#include "restraint/bcl_restraint_sas_transformation.h"
#include "restraint/bcl_restraint_saxs_opti_result.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AnalyzeSas::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new AnalyzeSas())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AnalyzeSas::AnalyzeSas() :
        m_OutFilePostFix( ".AnalyzeSas"),
        m_ExperimentalDataFileName(),
        m_SasaDataFileName(),
        m_ComputedDataFileName(),
        m_ExcludedVolumeParameter( 1.0),
        m_HydrationShellThickness( 0.0),
        m_MaximumDimension( util::GetUndefined< double>()),
        m_OptimizeHydrationShellParameters( false),
        m_C1Min( SasOptimization::default_C1Min),
        m_C1Max( SasOptimization::default_C1Max),
        m_C2Min( SasOptimization::default_C2Min),
        m_C2Max( SasOptimization::default_C2Max),
        m_C1StepSize( SasOptimization::default_C1StepSize),
        m_C2StepSize( SasOptimization::default_C2StepSize),
        m_ScoreType( "chi"),
        m_Cpu( false),
        m_ShouldApproximateSideChains( true),
        m_ShouldApproximateLoops( false),
        m_ShouldUseSans( false),
        m_DeuteriumExchangeParameter( 0.0)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AnalyzeAtomDistanceHeatmap
    AnalyzeSas *AnalyzeSas::Clone() const
    {
      return new AnalyzeSas( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AnalyzeSas::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &AnalyzeSas::GetOutFilePostfix() const
    {
      return m_OutFilePostFix;
    }

    //! @brief gives the experimental data file name
    //! @return the filename string
    const std::string &AnalyzeSas::GetExperimentalDataFileName() const
    {
      return m_ExperimentalDataFileName;
    }

    //! @brief gives the sasa data file name
    //! @return the filename string
    const std::string &AnalyzeSas::GetSasaDataFileName() const
    {
      return m_SasaDataFileName;
    }

    //! @brief gives the computed data file name
    //! @return the filename string
    const std::string &AnalyzeSas::GetComputedDataFileName() const
    {
      return m_ComputedDataFileName;
    }

    //! @brief gives the excluded volume parameter
    //! @return the c1 parameter
    const float &AnalyzeSas::GetExcludedVolumeParameter() const
    {
      return m_ExcludedVolumeParameter;
    }

    //! @brief gives the hydration shell thickness parameter
    //! @return the c2 parameter
    const float &AnalyzeSas::GetHydrationShellParameter() const
    {
      return m_HydrationShellThickness;
    }

    //! @brief gives the maximum dimension parameter
    //! @return dmax
    const float &AnalyzeSas::GetMaximumDimension() const
    {
      return m_MaximumDimension;
    }

    //! @brief gives the Minimum value for C1 on the search grid
    //! @return C1Min
    const float &AnalyzeSas::GetC1Min() const
    {
      return m_C1Min;
    }

    //! @brief gives the Maximum value for C1 on the search grid
    //! @return C1Max
    const float &AnalyzeSas::GetC1Max() const
    {
      return m_C1Max;
    }

    //! @brief gives the Minimum value for C2 on the search grid
    //! @return C2Min
    const float &AnalyzeSas::GetC2Min() const
    {
      return m_C2Min;
    }

    //! @brief gives the Maximum value for C2 on the search grid
    //! @return C2Max
    const float &AnalyzeSas::GetC2Max() const
    {
      return m_C1Max;
    }

    //! @brief gives the StepSize for C1 on the search grid
    //! @return C1StepSize
    const float &AnalyzeSas::GetC1StepSize() const
    {
      return m_C1StepSize;
    }

    //! @brief gives the StepSize for C2 on the search grid
    //! @return C2StepSize
    const float &AnalyzeSas::GetC2StepSize() const
    {
      return m_C2StepSize;
    }

    const float &AnalyzeSas::GetDeuteriumExchangeParameter() const
    {
      return m_DeuteriumExchangeParameter;
    }

    //! @brief gives name of the score to be used to compare two saxs profiles
    //! @return Score type
    const std::string &AnalyzeSas::GetScoreTypeName() const
    {
      return m_ScoreType;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &AnalyzeSas::GetAlias() const
    {
      static const std::string s_Name( "AnalyzeSas");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string AnalyzeSas::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {

      // This will need to find a new place to reside.  This provides the ability to read a raw
      // ExperimentalandCalculated Data file and perform any transformation desired without re-computing the
      // sas profiles

      BCL_Assert( !m_ExperimentalDataFileName.empty(), "experimental data file name is not defined");

      // Read and store the Experimental Data.  The SAS curve generated from the protein model will be fit to this data
      util::ShPtr< SasScatteringData> sp_experimental_data( new SasScatteringData());

      io::IFStream read;
      io::File::MustOpenIFStream( read, m_ExperimentalDataFileName);
      sp_experimental_data->ReadFromDataFile( read);
      io::File::CloseClearFStream( read);

      double max_intensity( SasAnalysis::MaxIntensity( *sp_experimental_data));

      if( !m_ComputedDataFileName.empty())
      {
        // Read in ExperimentalandCalculated Data file
        util::ShPtr< SasExperimentalAndCalculatedData> sp_SasExperimentalAndCalculatedData
        (
          new SasExperimentalAndCalculatedData()
        );

        io::IFStream read;
        io::File::MustOpenIFStream( read, m_ComputedDataFileName);
        sp_SasExperimentalAndCalculatedData->ReadFromDataFile( read);
        io::File::CloseClearFStream( read);

        // Transform the data
        SasExperimentalAndCalculatedData transformed_data( m_Transform( *sp_SasExperimentalAndCalculatedData));

        // Score the transformed data
        double result
        (
          score::SasType( m_Transform.GetUseErrors(), m_ScoreType)( transformed_data)
        );

        const storage::Vector< SasTransformation::TransformationTypeEnum> transforms_performed
        (
          m_Transform.GetTransformationTypes()
        );

        std::string transforms;

        for
        (
          storage::Vector< SasTransformation::TransformationTypeEnum>::const_iterator
           transform_itr( transforms_performed.Begin()),
           transform_itr_end( transforms_performed.End());
           transform_itr != transform_itr_end;
          ++transform_itr
        )
        {
          transforms += " " + transform_itr->GetString();
        }

        std::string summary
        (
          "\n ScoreFunction: " + util::Format()( m_ScoreType) + "\n"
          " Max Intensity: " + util::Format()( max_intensity) + "\n"
          " Use errors: " + util::Format()( m_Transform.GetUseErrors()) + "\n"
          " Transforms: " + util::Format()( transforms) + "\n"
          " Score: " + util::Format()( result) + "\n"
          " Hardware: " + util::Format()( m_Cpu) + "\n"
        );

        BCL_MessageStd( summary);

        return util::Format()( result);
      }

      double result( 0.0);

      BCL_Assert( ENSEMBLE.GetSize() != 0, "protein ensemble is empty");

      // iterate over all models in the pdbs.ls file
      for
      (
        assemble::ProteinEnsemble::const_iterator
        model_itr( ENSEMBLE.Begin()), model_itr_end( ENSEMBLE.End());
        model_itr != model_itr_end;
        ++model_itr
      )
      {
        util::ShPtr< assemble::ProteinModel> sp_protein_model( *model_itr);
        assemble::ProteinModel protein_model( *sp_protein_model);

        // Get a shared pointer to the protein model data
        util::ShPtr< assemble::ProteinModelData> sp_protein_model_data( protein_model.GetProteinModelData());

        // If Sasa File Provided insert SASA data into Protein Model
        if( !m_SasaDataFileName.empty())
        {
          // Sasa file must be preprocessed from PDBConvert and MSMS
          util::ShPtr< biol::SasaData> sp_sasa_data( new biol::SasaData());

          io::IFStream read_sasa;
          io::File::MustOpenIFStream( read_sasa, m_SasaDataFileName);
          sp_sasa_data->ReadFromDataFile( read_sasa);
          io::File::CloseClearFStream( read_sasa);

          // Insert the Sasa Data into the Protein Model Data
          sp_protein_model_data->Insert( assemble::ProteinModelData::e_Sasa, sp_sasa_data);
        }

        //Insert the Protein Model Data into the protein model
        protein_model.SetProteinModelData( sp_protein_model_data);

        // Create Local Variable to hold C1 and C2 values.  If The values are optimized, update the variable
        float c1_value( m_ExcludedVolumeParameter), c2_value( m_HydrationShellThickness);

        if( m_OptimizeHydrationShellParameters)
        {
          BCL_Assert( !m_SasaDataFileName.empty(), "sasa data file must be defined to optimize parameters");

          // Create default object
          SasOptimization optimization_object;

          if( !m_DefaultSearchGrid)
          {
            // setup grid optmization with the score type and the error flag
            optimization_object.SetC1Min( m_C1Min);
            optimization_object.SetC1Max( m_C1Max);
            optimization_object.SetC2Min( m_C2Min);
            optimization_object.SetC2Max( m_C2Max);
            optimization_object.SetC1StepSize( m_C1StepSize);
            optimization_object.SetC2StepSize( m_C2StepSize);
          }
          optimization_object.SetScoreFunction( m_ScoreType);
          optimization_object.SetErrorFlag( m_Transform.GetUseErrors());
          optimization_object.SetTransformationTypes( m_Transform);
          optimization_object.SetExperimentalData( sp_experimental_data);
          optimization_object.SetApproximateSideChainsFlag( m_ShouldApproximateSideChains);
          optimization_object.SetHardwareType( m_Cpu);
          optimization_object.SetSasType( m_ShouldUseSans);
          optimization_object.SetDeuteriumExchangeParameter( m_DeuteriumExchangeParameter);

          // Optimize c1 and c2 values. The experimental and sasa data were previously inserted into the protein_model
          SaxsOptiResult optimal_parameters( optimization_object( protein_model));

          // Update c1 and c2 values with optimized parameters
          c1_value = optimal_parameters.GetC1();
          c2_value = optimal_parameters.GetC2();
        }

        // Setup Commandline Strings for either the opencl or non-opencl version of the code
        util::Implementation< SasDebyeInterface> saxs
        (
          SasAnalysis::SetDebyeImplementation
          (
            m_ShouldApproximateLoops,
            m_ShouldApproximateSideChains,
            c1_value,
            c2_value,
            m_Cpu,
            m_ShouldUseSans,
            m_DeuteriumExchangeParameter
          )
        );

        // Set the Experimental Data
        saxs->SetExperimentalData( sp_experimental_data);

        //Create container to hold experimental and computed SAXS profiles and store results in container
        util::ShPtr< SasExperimentalAndCalculatedData> sp_experimental_and_calculated_data
        (
          new SasExperimentalAndCalculatedData( saxs->operator()( protein_model))
        );

        // Transform the data
        SasExperimentalAndCalculatedData transformed_data( m_Transform( *sp_experimental_and_calculated_data));

        // Score the transformed data
        result = score::SasType( m_Transform.GetUseErrors(), m_ScoreType)( transformed_data);

        const storage::Vector< SasTransformation::TransformationTypeEnum> transforms_performed
        (
          m_Transform.GetTransformationTypes()
        );

        std::string transforms;

        for
        (
          storage::Vector< SasTransformation::TransformationTypeEnum>::const_iterator
            transform_itr( transforms_performed.Begin()),
            transform_itr_end( transforms_performed.End());
            transform_itr != transform_itr_end;
          ++transform_itr
        )
        {
          transforms += " " + transform_itr->GetString();
        }

        // get the name of the pdb the structure came from
        util::ShPtr< util::Wrapper< std::string> > pdb_filename
        (
           ( *model_itr)->GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile)
        );

        std::string pdb( "dummy.pdb");

        if( pdb_filename.IsDefined())
        {
          pdb = std::string( pdb_filename->GetData());
        }

        std::string summary
        (
          "\n ScoreFunction: " + util::Format()( m_ScoreType) + "\n"
          " Max_Intensity: " + util::Format()( max_intensity) + "\n"
          " Use errors: " + util::Format()( m_Transform.GetUseErrors()) + "\n"
          " C1: " + util::Format()( c1_value) + "\n"
          " C2: " + util::Format()( c2_value) + "\n"
          " Transforms: " + util::Format()( transforms) + "\n"
          " Score: " + util::Format()( result) + " " + util::Format()( pdb) + "\n"
          " Hardware: " + util::Format()( m_Cpu) + "\n"
        );

        BCL_MessageStd( summary);
      }

      return util::Format()( result);

      // return the stringSasAnalysis
      //return results;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief hook that derived classes can override, called after TryRead successfully reads a serializer containing only initializer info (no data variables)
    //! @param SERIALIZER the serializer object with initialization information
    //! @param ERR_STREAM stream to write out errors to
    //! @return true, unless there were new errors
    bool AnalyzeSas::ReadInitializerSuccessHook( const util::ObjectDataLabel &SERIALIZER, std::ostream &ERR_STREAM)
    {
      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AnalyzeSas::GetSerializer() const
    {
      io::Serializer parameters( m_Transform.GetSerializer());
      parameters.SetClassDescription
      (
        "Computes a SAXS profile from PDB or BCL model and compares the curve with Experimental SAXS Data."
        "hydration layer controlled by excluded solvent c1, and border thickness c2 variables."
        "The Solvent assessible surface area (SASA) file is read in from MSMS"
      );

      parameters.AddInitializer
      (
        "filename_postfix",
        "the postfix that will be appended to the output file in order to identify this analysis",
        io::Serialization::GetAgent( &m_OutFilePostFix),
        ".AnalyzeSas"
      );

      parameters.AddInitializer
      (
        "experimental_profile",
        "the experimental SAXS data",
        io::Serialization::GetAgentInputFilename( &m_ExperimentalDataFileName)
      );

      parameters.AddOptionalInitializer
      (
        "sasa_profile",
        "the solvent accessible surface area of the protein model",
        io::Serialization::GetAgentInputFilename( &m_SasaDataFileName)
      );

      parameters.AddOptionalInitializer
      (
        "computed_profile",
        "the computed SAXS data",
        io::Serialization::GetAgentInputFilename( &m_ComputedDataFileName)
      );

      parameters.AddInitializer
      (
        "c1",
        "the excluded volume tuning parameter",
        io::Serialization::GetAgent( &m_ExcludedVolumeParameter)
      );

      parameters.AddInitializer
      (
        "c2",
        "the hydration shell thickness",
        io::Serialization::GetAgent( &m_HydrationShellThickness)
      );

      parameters.AddOptionalInitializer
      (
        "c1_min",
        "the minimum excluded volume parameter value on the grid search",
        io::Serialization::GetAgent( &m_C1Min)
      );

      parameters.AddOptionalInitializer
      (
        "c1_max",
        "the maximum excluded volume parameter value on the grid search",
        io::Serialization::GetAgent( &m_C1Max)
      );

      parameters.AddOptionalInitializer
      (
        "c2_min",
        "the minimum hydration thickness parameter value on the grid search",
        io::Serialization::GetAgent( &m_C2Min)
      );

      parameters.AddOptionalInitializer
      (
        "c2_max",
        "the maximum hydration thickness parameter value on the grid search",
        io::Serialization::GetAgent( &m_C2Max)
      );

      parameters.AddOptionalInitializer
      (
        "c1_stepsize",
        "the excluded volume parameter stepsize on the grid search",
        io::Serialization::GetAgent( &m_C1StepSize)
      );

      parameters.AddOptionalInitializer
      (
        "c2_stepsize",
        "the hydration shell parameter stepsize on the grid search",
        io::Serialization::GetAgent( &m_C2StepSize)
      );

      parameters.AddOptionalInitializer
      (
        "dmax",
        "the maximum dimension of the particle from Gnome",
        io::Serialization::GetAgent( &m_MaximumDimension)
      );

      parameters.AddOptionalInitializer
      (
        "optimize_hydration_parameters",
        "1 if the hydration parameters should be optimized - 0 otherwise",
        io::Serialization::GetAgent( &m_OptimizeHydrationShellParameters)
      );

      parameters.AddInitializer
      (
        "default_search_grid",
        "1 to use default parameters- 0 otherwise",
        io::Serialization::GetAgent( &m_DefaultSearchGrid)
      );

      parameters.AddInitializer
      (
        "scoring_function",
        "The scoring function to use for analysis",
        io::Serialization::GetAgent( &m_ScoreType)
      );

      parameters.AddOptionalInitializer
       (
         "cpu",
         "true to force cpu implementation of debye formula - false otherwise",
         io::Serialization::GetAgent( &m_Cpu)
       );

      parameters.AddInitializer
       (
         "approximate_side_chains",
         "true sum form factors to cb atom - false otherwise",
         io::Serialization::GetAgent( &m_ShouldApproximateSideChains)
       );

      parameters.AddInitializer
       (
         "approximate_loops",
         "true approximate loops - false otherwise",
         io::Serialization::GetAgent( &m_ShouldApproximateLoops)
       );

      parameters.AddInitializer
       (
         "use_sans",
         "true use sans implementation - false use saxs implementation",
         io::Serialization::GetAgent( &m_ShouldUseSans)
       );

      parameters.AddOptionalInitializer
      (
        "deuterium_percent",
        "the percent of deuterium in the solvent",
        io::Serialization::GetAgent( &m_DeuteriumExchangeParameter)
      );

      return parameters;
    }
  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_atom_distance_assignment.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AtomDistanceAssignment::s_Instance
    (
      GetObjectInstances().AddInstance( new AtomDistanceAssignment())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AtomDistanceAssignment::AtomDistanceAssignment() :
      m_AtomA(),
      m_AtomB(),
      m_Distance()
    {
    }

    //! @brief construct from atoms and distance
    //! @param ATOM_A first atom
    //! @param ATOM_B second atom
    //! @param DISTANCE experimental distance
    AtomDistanceAssignment::AtomDistanceAssignment
    (
      const biol::Atom &ATOM_A,
      const biol::Atom &ATOM_B,
      const util::ShPtr< Distance> &DISTANCE
    ) :
      m_AtomA( ATOM_A),
      m_AtomB( ATOM_B),
      m_Distance( DISTANCE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AtomDistanceAssignment
    AtomDistanceAssignment *AtomDistanceAssignment::Clone() const
    {
      return new AtomDistanceAssignment( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AtomDistanceAssignment::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief calculates the distance between the member atoms
    //! @return the distance between the member atoms
    double AtomDistanceAssignment::CalculateAtomDistance() const
    {
      return linal::Distance( m_AtomA.GetCoordinates(), m_AtomB.GetCoordinates());
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AtomDistanceAssignment::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_AtomA, ISTREAM);
      io::Serialize::Read( m_AtomB, ISTREAM);
      io::Serialize::Read( m_Distance, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AtomDistanceAssignment::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_AtomA, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AtomB, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Distance, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
  
} // namespace bcl
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
#include "restraint/bcl_restraint_atom_distance.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_atom.h"
#include "restraint/bcl_restraint_atom_distance_assignment.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AtomDistance::s_Instance
    (
      GetObjectInstances().AddInstance( new AtomDistance())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AtomDistance::AtomDistance() :
      m_DataPair(),
      m_Distance(),
      m_Confidence( 1.0)
    {
    }

    //! @brief construct from locators and distance
    //! @param LOCATOR_A first atom locator
    //! @param LOCATOR_B second atom locator
    //! @param DISTANCE experimental distance information
    //! @param CONFIDENCE in the given restraint
    AtomDistance::AtomDistance
    (
      const assemble::LocatorAtom &LOCATOR_A,
      const assemble::LocatorAtom &LOCATOR_B,
      const util::ShPtr< Distance> &DISTANCE,
      const double &CONFIDENCE
    ) :
      m_DataPair
      (
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>( LOCATOR_A.Clone()),
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>( LOCATOR_B.Clone())
      ),
      m_Distance( DISTANCE),
      m_Confidence( CONFIDENCE)
    {
    }

    //! @brief construct from locators and distance
    //! @param LOCATOR_A first atom locator
    //! @param LOCATOR_B second atom locator
    //! @param DISTANCE experimental distance information
    //! @param CONFIDENCE in the given restraint
    AtomDistance::AtomDistance
    (
      const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &LOCATOR_A,
      const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &LOCATOR_B,
      const util::ShPtr< Distance> &DISTANCE,
      const double &CONFIDENCE
    ) :
      m_DataPair( LOCATOR_A, LOCATOR_B),
      m_Distance( DISTANCE),
      m_Confidence( CONFIDENCE)
    {
    }

    //! @brief construct from data pair and distance
    //! @param DATA_PAIR the data pair corresponding to this restraint
    //! @param DISTANCE experimental distance information
    //! @param CONFIDENCE in the given restraint
    AtomDistance::AtomDistance
    (
      const DataPairwise &DATA_PAIR,
      const util::ShPtr< Distance> &DISTANCE,
      const double &CONFIDENCE
    ) :
      m_DataPair( DATA_PAIR),
      m_Distance( DISTANCE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AtomDistance
    AtomDistance *AtomDistance::Clone() const
    {
      return new AtomDistance( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AtomDistance::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief writes the restraint information in a compact format
    //! @return the restraint information in a compact format
    std::string AtomDistance::GetIdentification() const
    {
      return m_DataPair.GetIdentification() + "\t" + m_Distance->GetIdentification();
    }

    //! @brief data access to get confidence for this restraint
    //! @return double which is the confidence in this restraint
    double AtomDistance::GetConfidence() const
    {
      return m_Confidence;
    }

    //! @brief data access to get the data pair the distance corresponds to
    //! @return restraint::DataPairwise which is m_DataPair
    const DataPairwise &AtomDistance::GetData() const
    {
      return m_DataPair;
    }

    //! @brief data access to get distance contained in the distance object
    //! @return double which is the distance contained in the distance object
    util::ShPtr< Distance> AtomDistance::GetDistance() const
    {
      return m_Distance;
    }

    //! @brief data access to get upper bound of the distance contained in the distance object
    //! @return double which is the upper bound of the distance contained in the distance object
    double AtomDistance::GetUpperBound() const
    {
      return m_Distance->UpperBound();
    }

    //! @brief data access to get lower bound of the distance contained in the distance object
    //! @return double which is the lower bound of the distance contained in the distance object
    double AtomDistance::GetLowerBound() const
    {
      return m_Distance->LowerBound();
    }

    //! @brief data access to get upper error of the distance contained in the distance object
    //! @return double which is the upper error of the distance contained in the distance object
    double AtomDistance::GetUpperError() const
    {
      return m_Distance->GetUpperError();
    }

    //! @brief data access to get lower error of the distance contained in the distance object
    //! @return double which is the lower error of the distance contained in the distance object
    double AtomDistance::GetLowerError() const
    {
      return m_Distance->GetLowerError();
    }

    //! @brief returns whether the restraint is defined
    //! @return whether the restraint is defined
    bool AtomDistance::IsDefined() const
    {
      return m_DataPair.IsSet() && m_Distance->IsDefined();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief generates the assignment from the protein model, undefined atoms are passed if they could not
    //!        be located, if the atom is a hydrogen a new atom is created at the CB position
    //!        using the H atom type
    //! @param PROTEIN_MODEL protein model to be used to generate the assignment
    //! @return assignment containing located atoms and experimental distance
    AtomDistanceAssignment AtomDistance::GenerateAssignment( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      BCL_Assert( m_DataPair.IsSet(), "data pair is not set yet");

      // create the assignment and return
      return AtomDistanceAssignment
      (
        m_DataPair.First()->LocateAtomCopy( PROTEIN_MODEL),
        m_DataPair.Second()->LocateAtomCopy( PROTEIN_MODEL),
        m_Distance
      );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AtomDistance::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_DataPair, ISTREAM);
      io::Serialize::Read( m_Distance, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AtomDistance::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_DataPair, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Distance, OSTREAM, INDENT) << '\n';

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint

} // namespace bcl
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
#include "restraint/bcl_restraint_body.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "restraint/bcl_restraint_assignment.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Body::s_Instance
    (
      GetObjectInstances().AddInstance( new Body())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Body::Body() :
      m_Bodies(),
      m_DetermineOccupancy()
    {
    }

    //! @brief construct from bodies
    //! @param BODIES the ShPtrVector of bodies which will serve as restraints
    //! @param DETERMINE_OCCUPANCY is a BinaryFunctionInterface which is used to determine if a restraint body is occupied by a SSE
    Body::Body
    (
      const util::ShPtr< util::ShPtrVector< assemble::SSEGeometryInterface> > &BODIES,
      const util::ShPtr< util::BinaryFunctionInterface< assemble::SSEGeometryInterface, assemble::SSE, bool> >
        &DETERMINE_OCCUPANCY
    ) :
      m_Bodies( BODIES),
      m_DetermineOccupancy( DETERMINE_OCCUPANCY)
    {
    }

    //! @brief virtual copy constructor
    Body *Body::Clone() const
    {
      return new Body( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Body::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief GetBody returns a const reference to "m_Bodies"
    //! @return returns a const reference to "m_Bodies"
    const util::ShPtr< util::ShPtrVector< assemble::SSEGeometryInterface> > &Body::GetBody() const
    {
      return m_Bodies;
    }

    //! @brief SetBody "m_Bodies" to a new ShPtr of bodies
    //! @param BODIES is the ShPtrVector of bodies which "m_Bodies" will be changed to
    void Body::SetBody( const util::ShPtr< util::ShPtrVector< assemble::SSEGeometryInterface> > &BODIES)
    {
      m_Bodies = BODIES;
    }

  ///////////////
  // operation //
  ///////////////

    //! @brief GetUnoccupied determines which of the bodies in "m_Bodies" is not occupied
    //! @param SSES a ShPtrVector of sses which could occupy "m_Bodies"
    //! @return returns a ShPtrVector of bodies which are not occupied by any of "BODIES"
    util::ShPtrVector< assemble::SSEGeometryInterface>
    Body::GetUnoccupied( const util::SiPtrVector< const assemble::SSE> &SSES) const
    {
      BCL_Assert( m_DetermineOccupancy.IsDefined(), "restraint::Body::GetUnoccupied m_DetermineOccupancy is not defined");

      // create ShPtrVector "unoccupied" to hold bodies of "m_Bodies" that are not occupied by "BODIES"
      util::ShPtrVector< assemble::SSEGeometryInterface> unoccupied;

      // iterate through bodies
      for
      (
        util::ShPtrVector< assemble::SSEGeometryInterface>::const_iterator
          itr_body( m_Bodies->Begin()), itr_body_end( m_Bodies->End());
        itr_body != itr_body_end;
        ++itr_body
      )
      {
        bool is_occupied( false);

        // iterate over all SSEs and check if it occupies current body
        for
        (
          util::SiPtrVector< const assemble::SSE>::const_iterator
            itr_sse( SSES.Begin()), itr_sse_end( SSES.End());
          itr_sse != itr_sse_end;
          ++itr_sse
        )
        {
          // true if the body denoted by "itr_sse" falls within the body denoted by "itr_body"
          if( m_DetermineOccupancy->operator()( **itr_body, **itr_sse))
          {
            // set "is_occupied" to true since the body denoted by "itr_body" is occupied
            is_occupied = true;

            // leave this for loop
            break;
          }
        }

        // true if the body denoted by "itr_body" is unoccupied
        if( !is_occupied)
        {
          // add the unoccupied body denoted by "itr_body" to "unoccupied"
          unoccupied.PushBack( *itr_body);
        }
      }

      BCL_MessageDbg
      (
        "restraint::Body::GetUnoccupied unoccupied size " + util::Format()( unoccupied.GetSize())
      );

      // return "unoccupied"
      return unoccupied;
    }

    //! @brief get body that is occupied by given SSE
    //! @param SSE the sse to be considered to identify the body that is occupied by it
    //! @return SHPtr to occupied body - will be undefined if there is non
    util::ShPtr< assemble::SSEGeometryInterface> Body::GetOccupiedBody( const assemble::SSE &SSE) const
    {
      // iterate through bodies to find which one is occupied by the SSE
      for
      (
        util::ShPtrVector< assemble::SSEGeometryInterface>::const_iterator
          itr_body( m_Bodies->Begin()), itr_body_end( m_Bodies->End());
        itr_body != itr_body_end;
        ++itr_body
      )
      {
        // true if the body denoted by "itr_body" is occupied by "BODY"
        if( m_DetermineOccupancy->operator()( **itr_body, SSE))
        {
          return *itr_body;
        }
      }

      return util::ShPtr< assemble::SSEGeometryInterface>();
    }

    //! @brief GenerateAssignment creates the assignment of "m_Bodies" with other assemble::SSEs
    //! @param SSES SiPtrVector of SSE which will be assigned with "m_Bodies"
    //! @return returns an Assignment which assigns "m_Bodies" with the appropriate members of "SSES"
    const SSEAssignment
    Body::GenerateAssignment( const util::SiPtrVector< const assemble::SSE> &SSES) const
    {
      // create GroupCollection "group_collection" to hold elements of BODIES which occupy "m_Bodies"
      GroupCollection< size_t, assemble::SSE> occupied_group_collection;
      BCL_Assert( m_DetermineOccupancy.IsDefined(), "restraint::Body::GenerateAssignment m_DetermineOccupancy is not defined");

      // iterate through "m_Bodies" for instance the density rods
      for
      (
        util::ShPtrVector< assemble::SSEGeometryInterface>::const_iterator itr_body( m_Bodies->Begin()), itr_body_end( m_Bodies->End());
        // util::ShPtrVector< assemble::SSEGeometryInterface>::const_iterator itr_body( m_Bodies.Begin()), itr_body_end( m_Bodies.End());
        itr_body != itr_body_end;
        ++itr_body
      )
      {
        // iterate through SSES
        for
        (
          util::SiPtrVector< const assemble::SSE>::const_iterator
            itr_sse( SSES.Begin()), itr_sse_end( SSES.End());
          itr_sse != itr_sse_end;
          ++itr_sse
        )
        {
          const bool occupants( m_DetermineOccupancy->operator()( **itr_body, **itr_sse));
          BCL_MessageDbg
          (
            "origin of restraint bodies (e.g. density rod): " + util::Format()( ( *itr_body)->GetCenter()) +
            ", origin of argument (e.g. protein model SSE): " + util::Format()( ( *itr_sse)->GetCenter()) +
            ", occupied?: " + util::Format()( occupants)
          );
          if( occupants)
          {
            BCL_MessageDbg
            (
              "start point of density rod " + util::Format()( ( *itr_body)->GetCenter() - ( *itr_body)->GetAxis( coord::GetAxes().e_Z) * ( *itr_body)->GetExtent( coord::GetAxes().e_Z)) +
              " end point of density rod " + util::Format()( ( *itr_body)->GetCenter() + ( *itr_body)->GetAxis( coord::GetAxes().e_Z) * ( *itr_body)->GetExtent( coord::GetAxes().e_Z)) +
              " radius of density rod " + util::Format()( ( *itr_body)->GetExtent( coord::GetAxes().e_X))
            );
            BCL_MessageDbg
            (
              "start point of sse " + util::Format()( ( *itr_sse)->GetCenter() - ( *itr_sse)->GetAxis( coord::GetAxes().e_Z) * ( *itr_sse)->GetExtent( coord::GetAxes().e_Z)) +
              " end point of sse " + util::Format()( ( *itr_sse)->GetCenter() + ( *itr_sse)->GetAxis( coord::GetAxes().e_Z) * ( *itr_sse)->GetExtent( coord::GetAxes().e_Z)) +
              " radius of density rod " + util::Format()( ( *itr_sse)->GetExtent( coord::GetAxes().e_X))
            );
          }
          BCL_MessageDbg( "occupancy gotten\n ");

          // true if the body denoted by "itr_sse" occupies the body denoted by "itr_body" where occupancy is determined by "m_DetermineOccupancy"
          // for example if the origins of the body are the same
          if( occupants)
          {
            //BCL_Assert( math::EqualWithinTolerance( ( *itr_body)->GetOrigin(), ( *itr_sse)->GetOrigin()), "origin of bodies is not the same");
            BCL_MessageDbg( "getting index\n ");

            const size_t index( itr_body - m_Bodies->Begin());
            BCL_MessageDbg( "index is " + util::Format()( index));

            // insert the element of "BODIES" into "occupied_group_collection"
            occupied_group_collection.Insert( index, Group< assemble::SSE>( 1, *itr_sse));

            // leave this for loop since the occupant of the body denoted by "itr_body" has been found
            break;
          }
        }
      }

//        BCL_Assert
//        (
//          occupied_group_collection.TotalDepth() == SSES.GetSize(),
//            "\n \nevery sse in the protein model should be within a restraint body, therefore the number of bodies in SSES "
//            "(from the protein models) should equal the number of occupied restraint bodies \n but the number of occupied restraint bodies is \n"
//          + util::Format()( occupied_group_collection.TotalDepth())
//          + "\n and the number of sses in the protein model is "
//          + util::Format()( SSES.GetSize())
//          + "\n\n the occupied_group_collection (corresponding to occupied density rods) is "
//          + util::Format()( occupied_group_collection)
//          + "\n while the protein model contains the following SSEs "
//          + util::Format()( SSES)
//        );

      // return Assignment constructed with "m_Bodies" and "occupied_group_collection"
      return SSEAssignment
      (
        m_Bodies, occupied_group_collection
      );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read distance from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Body::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Bodies, ISTREAM);
      io::Serialize::Read( m_DetermineOccupancy, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write distance to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &Body::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Bodies, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DetermineOccupancy, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl

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
#include "restraint/bcl_restraint_cone_model.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_collector_aa_type.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "quality/bcl_quality_rmsd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ConeModel::s_Instance
    (
      GetObjectInstances().AddInstance( new ConeModel())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ConeModel::ConeModel()
    {
    }

    //! @brief Clone function
    //! @return pointer to new ConeModel
    ConeModel *ConeModel::Clone() const
    {
      return new ConeModel( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ConeModel::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

      //! @brief calculates the maximum angle (SL->CB->SL) observed for an ensemble of models
      //! @param ENSEMBLE ensemble for which the maximum SL->CB->SL angle will be determined
      //! @return the maximum SL->CB->SL angle observed in the provided ensemble
      double ConeModel::SLCBSLMaxAngle( const assemble::ProteinEnsemble &ENSEMBLE)
      {

        // will hold the maximum SL->SB->SL angle seen for this ensemble
        double sl_cb_sl_max_angle( 0);

        // iterate through the ensemble
        for
        (
          assemble::ProteinEnsemble::const_iterator model_itr_a( ENSEMBLE.Begin()), model_itr_end( ENSEMBLE.End());
          model_itr_a != model_itr_end; ++model_itr_a
        )
        {
          // get the spin label residue from model_a
          const util::SiPtr< const biol::AABase> sl_a( GetSpinLabelResidue( **model_itr_a));

          // get the atom coordinates of the spin label side chain
          const util::SiPtrVector< const linal::Vector3D> sl_a_coords
          (
            sl_a->GetAtomCoordinates( GetSuperimposeAtomTypes())
          );

          // get the unpaired electron coordinates from spin label a
          const linal::Vector3D unpaired_electron_a_coords( GetUnpairedElectronCoordinates( *sl_a));

          // iterate through the ensemble a second time
          assemble::ProteinEnsemble::const_iterator model_itr_b( model_itr_a);
          ++model_itr_b;
          for( ; model_itr_b != model_itr_end; ++model_itr_b)
          {
            // get the spin label coordinates from spin label b copy that has been superimposed on spin label a
            const linal::Vector3D unpaired_electron_b_coords( GetSuperimposedUnpairedElectronCoordinates( **model_itr_b, sl_a_coords));

            // calculate the angle between the CB of sl_a, and the electron in sl_a and sl_b_copy
            const double proj_angle
            (
              linal::ProjAngle
              (
                sl_a->GetAtom( biol::GetAtomTypes().CB).GetCoordinates(),
                unpaired_electron_a_coords,
                unpaired_electron_b_coords
              )
            );

            BCL_MessageDbg( "current SL->CB-SL angle is " + util::Format()( proj_angle));

            // true if the current SL->CB->SL angle is larger than max angle seen so far
            if( proj_angle > sl_cb_sl_max_angle)
            {
              // set the max angle to the current projection angle
              sl_cb_sl_max_angle = proj_angle;
            }
          }
        }

        // return the maximum SL->CB->SL angle observed in the ensemble
        return sl_cb_sl_max_angle;
      }

      //! @brief calculates the angle between the effective spin label position, the CB, and the CA
      //!        SLeffective->CB->CA
      //! @param ENSEMBLE the ensemble of protein models from which the SLeffective->CB->CA
      //! @return double which is the calculated SLeffective->CB->CA angle
      double ConeModel::SLeffectiveCBCAAngle( const assemble::ProteinEnsemble &ENSEMBLE)
      {
        storage::Pair< linal::Vector3D, util::SiPtr< const biol::AABase> >
          slef_tmpltaa( GetSLEffectivePositionAndTemplateResidue( ENSEMBLE));

        // now calculate the SLeffective->CB->CA angle using the sl_a coordinates for CB and CA
        const double proj_angle
        (
          linal::ProjAngle
          (
            slef_tmpltaa.Second()->GetAtom( biol::GetAtomTypes().CB).GetCoordinates(),
            slef_tmpltaa.Second()->GetAtom( biol::GetAtomTypes().CA).GetCoordinates(),
            slef_tmpltaa.First()
          )
        );

        BCL_MessageDbg( "SLeffective->CB->CA angle is " + util::Format()( proj_angle));

        // return the calculated SLeffective->CB->CA angle
        return proj_angle;
      }

      //! @brief calculates the distance between the effective spin label position and the CB for an ensemble
      //! @param ENSEMBLE the protein ensemble for which the SL effective position will be calculated
      //! @return double which is the distance between the CB and SLeffective position
      double ConeModel::SLeffectiveCBDistance( const assemble::ProteinEnsemble &ENSEMBLE)
      {
        const storage::Pair< linal::Vector3D, util::SiPtr< const biol::AABase> > sleff_tmpltaa
        (
          GetSLEffectivePositionAndTemplateResidue( ENSEMBLE)
        );

        const double distance
        (
          linal::Distance
          (
            sleff_tmpltaa.First(), sleff_tmpltaa.Second()->GetAtom( biol::GetAtomTypes().CB).GetCoordinates()
          )
        );

        BCL_MessageDbg( "current SLeffective->CB distance is " + util::Format()( distance));

        return distance;
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
    std::istream &ConeModel::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ConeModel::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

      //! @brief calculates the position of the unpaired electron in a spin label residue
      //!        the midpoint of the bond between the oxygen and nitrogen in the nitroxide moiety in the ring
      //! @param RESIDUE the residue where the position of the spin label will be calculated
      //! @return vector 3d which has the coordinates of the spin label (i.e. unpaired electron)
      linal::Vector3D ConeModel::GetUnpairedElectronCoordinates( const biol::AABase &RESIDUE)
      {
        // get the coordinates for the oxygen and nitrogen in the nitroxide moiety
        const util::SiPtrVector< const linal::Vector3D> o_coords
        (
          RESIDUE.GetAtomCoordinates( storage::Set< biol::AtomType>( biol::GetAtomTypes().O1))
        );
        BCL_Assert
        (
          o_coords.GetSize() == 1, util::Format()( o_coords.GetSize()) + " O1 coords found in " +
          RESIDUE.GetIdentification() + "\n" + util::Format()( RESIDUE) + "\nBe sure you are using AAComplete"
        );
        const util::SiPtrVector< const linal::Vector3D> n_coords
        (
          RESIDUE.GetAtomCoordinates( storage::Set< biol::AtomType>( biol::GetAtomTypes().N1))
        );

        BCL_Assert
        (
          n_coords.GetSize() == 1, util::Format()( n_coords.GetSize()) + " N1 coords found in " +
          RESIDUE.GetIdentification() + "\n" + util::Format()( RESIDUE) + "\nBe sure you are using AAComplete"
        );

        // calculate the position of the unpaired electron
        const linal::Vector3D ave_pos( ( *o_coords.FirstElement() + *n_coords.FirstElement()) / 2.0);

        // return SL coordinates
        return ave_pos;
      }

      //! @brief gives the single spin label residue that is within a protein model - asserts if there is more than one
      //! @param MODEL the model which has one spin label residue which will be given
      //! @return SiPtr to the spin label residue that is within the protein model
      util::SiPtr< const biol::AABase> ConeModel::GetSpinLabelResidue
      (
        const assemble::ProteinModel &MODEL
      )
      {
        // will collect the spin label residues from the protein models
        const assemble::CollectorAAType sl_collector( storage::Set< biol::AAType>( biol::GetAATypes().R1A));

        // get the spin label residue from model
        const util::SiPtrList< const biol::AABase> sl_list_a( sl_collector.Collect( MODEL.GetAminoAcids()));

        // make sure exactly one spin label side residue was found
        BCL_Assert( sl_list_a.GetSize() == 1, util::Format()( sl_list_a.GetSize()) + " spin labels found");

        // get the one spin label residue out of the list
        const util::SiPtr< const biol::AABase> sl_a( sl_list_a.FirstElement());

        // return the coordinates of the unpaired electron
        return sl_a;
      }

      //! @brief gives the atom types that are used to superimpose spin label residues when the effective spin label
      //!        position is being determined
      //! @return set which has the atom types used for superimposition
      const storage::Set< biol::AtomType> &ConeModel::GetSuperimposeAtomTypes()
      {
        static const storage::Set< biol::AtomType> s_types
        (
          storage::Set< biol::AtomType>::Create
          (
            biol::GetAtomTypes().CA, biol::GetAtomTypes().CB, biol::GetAtomTypes().N, biol::GetAtomTypes().C,
            biol::GetAtomTypes().O
          )
        );

        return s_types;
      }

      //! @brief gets a spin label residue from a protein model and superimposes it onto the provided coordinates
      //!        and returns the resulting position of the unpaired electron in that spin label
      //! @param MODEL the model from which the spin label residue is going to be gotten
      //! @param COORDS_TO_SUPERIMPOSE_ON the coordinates to which the spin label residue will be superimposed
      //! @return Vector3D which has the coordinates of the unpaired electron in the superimposed SL residue
      linal::Vector3D ConeModel::GetSuperimposedUnpairedElectronCoordinates
      (
        const assemble::ProteinModel &MODEL, const util::SiPtrVector< const linal::Vector3D> &COORDS_TO_SUPERIMPOSE_ON
      )
      {
        // get the spin label residue from model
        const util::SiPtr< const biol::AABase> sl( GetSpinLabelResidue( MODEL));

        // get the atom coordinates of the spin label b side chain
        const util::SiPtrVector< const linal::Vector3D> sl_coords( sl->GetAtomCoordinates( GetSuperimposeAtomTypes()));

        // get the transformation matrix necessary to superimpose the CA, CB, N, and C atoms of the two sl residues
        // transformation matrix that will superimpose sl_b_coords onto sl_a_coords
        const math::TransformationMatrix3D transform( quality::RMSD::SuperimposeCoordinates( COORDS_TO_SUPERIMPOSE_ON, sl_coords));

        // copy spin label b and superimpose sl_b onto sl_a
        util::ShPtr< biol::AABase> sl_copy( sl->Clone());
        sl_copy->Transform( transform);

        // get the spin label coordinates from spin label b copy that has been superimposed on spin label a
        const linal::Vector3D unpaired_electron_coords( GetUnpairedElectronCoordinates( *sl_copy));

        return unpaired_electron_coords;
      }

      //! @brief gives the coordinate of the effective spin label (unpaired electron) position and the amino acid
      //!        used as the basis for superimposing all the spin labels into the same space
      //! @param ENSEMBLE the ensemble for which the effective spin label position will be calculated
      //! @return Pair which has the coordinate of the effective spin label and the spin label used for superimposition
      storage::Pair< linal::Vector3D, util::SiPtr< const biol::AABase> >
      ConeModel::GetSLEffectivePositionAndTemplateResidue( const assemble::ProteinEnsemble &ENSEMBLE)
      {
        linal::Vector3D sl_effective_pos( 0, 0, 0);

        util::SiPtr< const biol::AABase> sl_a( GetSpinLabelResidue( **ENSEMBLE.Begin()));

        // get the atom coordinates of the spin label side chain
        const util::SiPtrVector< const linal::Vector3D> sl_a_coords
        (
          sl_a->GetAtomCoordinates( GetSuperimposeAtomTypes())
        );

        sl_effective_pos += GetUnpairedElectronCoordinates( *sl_a);

        // iterate through the ensemble - starting with second protein, since first was already counted
        for
        (
          assemble::ProteinEnsemble::const_iterator model_itr( ++ENSEMBLE.Begin()), model_itr_end( ENSEMBLE.End());
          model_itr != model_itr_end; ++model_itr
        )
        {
          sl_effective_pos += GetSuperimposedUnpairedElectronCoordinates( **model_itr, sl_a_coords);
        }

        // get the average position of all the unpaired electrons, this is the effective SL position
        sl_effective_pos /= double( ENSEMBLE.GetSize());

        BCL_MessageDbg( "effective sl position is " + util::Format()( sl_effective_pos));

        return storage::Pair< linal::Vector3D, util::SiPtr< const biol::AABase> >( sl_effective_pos, sl_a);
      }

  } // namespace restraint

} // namespace bcl
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
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "fold/bcl_fold_score_weight_set.h"
#include "restraint/bcl_restraint_contact_data.h"
#include "score/bcl_score_restraint_atom_attraction.h"
#include "score/bcl_score_restraint_atom_distance.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    // initialize distance cutoff for contact
    const double ContactData::s_ContactDistance( 8.0);

    // initialize score
    fold::Score ContactData::e_ScoreContactRestraint( fold::GetScores().e_Undefined);

    //! contact restraint atom distance score
    const util::SiPtr< const score::RestraintAtomDistanceAssignment> ContactData::s_ScoreContactRestraintAssignment
    (
      util::Enumerated< score::RestraintAtomDistanceAssignment>::AddInstance
      (
        new score::RestraintAtomAttraction
        (
          score::RestraintAtomAttraction::GetDefaultDepthRange(),
          0,
          s_ContactDistance,
          true,
          score::RestraintAtomAttraction::GetDefaultScheme() + "_" + ContactData().GetAlias()
        )
      )
    );

    const util::SiPtr< const score::RestraintAtomDistance> ContactData::s_ScoreContactDistance
    (
      util::Enumerated< score::RestraintAtomDistance>::AddInstance
      (
        new score::RestraintAtomDistance
        (
          *s_ScoreContactRestraintAssignment,
          1.0,
          "contact_restraint"
        )
      )
    );

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ContactData::s_Instance
    (
      util::Enumerated< Interface>::AddInstance( new ContactData())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ContactData::ContactData() :
      m_Handler( ".RR", "-SASACCN", 0.0, 10.0, 8.0),
      m_Restraints(),
      m_Fraction( 1.0)
    {
    }

    //! @brief Clone function
    //! @return pointer to new ContactData
    ContactData *ContactData::Clone() const
    {
      return new ContactData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ContactData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ContactData::GetAlias() const
    {
      static const std::string s_name( "Contact");
      return s_name;
    }

    //! @brief get the default file extension
    //! @return the default file extension
    const std::string &ContactData::GetDefaultExtension() const
    {
      static const std::string s_extension( ".RR");
      return s_extension;
    }

    //! @brief gives reference to the specific type of data this restraint uses
    //! @return gives reference to the specific type of data this restraint uses
    const util::ShPtrVector< AtomDistance> &ContactData::GetAtomDistanceRestraints() const
    {
      return m_Restraints;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initialize the scores and add them to Scores enumerator
    void ContactData::InitializeScores()
    {
      if( !e_ScoreContactRestraint.IsDefined())
      {
        // read in restraints
        if( m_Restraints.IsEmpty())
        {
          m_Restraints = m_Handler.ReadRestraintsFromFile();
        }

        e_ScoreContactRestraint = fold::GetScores().AddScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::RestraintAtomDistance
            (
              *s_ScoreContactRestraintAssignment,
              m_Fraction,
              "contact_restraint",
              util::CloneToShPtr( m_Restraints)
            )
          )
        );
      }
    }

    //! @brief sets the weights of scores in a weight set
    //! @param SCORE_WEIGHT_SET the score weight set that will be modified
    void ContactData::ModifyScoreWeightSet( fold::ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
      SCORE_WEIGHT_SET.SetWeight( e_ScoreContactRestraint, 500);
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void ContactData::MergeAndModifyMutateTree( fold::MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< fold::MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads restraints formatted for this restraint type from an istream
    //! @return istream restraints formatted for this restraint type were read from
    std::istream &ContactData::ReadRestraints( std::istream &ISTREAM)
    {
      m_Restraints = m_Handler.ReadRestraints( ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ContactData::GetSerializer() const
    {
      io::Serializer serial( m_Handler.GetSerializer());
      serial.SetClassDescription( "Predicted or known AA or atom contacts (within 8A by default)");
      serial.AddInitializer
      (
        "fraction",
        "fraction of contacts that the protein model is expected to satisfy. "
        "The best scoring contacts (eq. to this fraction of the contacts in the file) will be used when computing the "
        "score. Useful when contacts are less than certain and the folding process is expected to satisfy only a "
        "fraction of them",
        io::Serialization::GetAgent( &m_Fraction),
        "1.0"
      );
      return serial;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ContactData::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Restraints, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ContactData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Restraints, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_contains_body_origin.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

    template class BinaryFunctionInterface< assemble::SSEGeometryInterface, assemble::SSE, bool>;

  } // namespace util

  namespace restraint
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ContainsBodyOrigin::s_Instance
    (
      GetObjectInstances().AddInstance( new ContainsBodyOrigin())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ContainsBodyOrigin::ContainsBodyOrigin()
    {
    }

    //! @brief virtual copy constructor
    ContainsBodyOrigin *ContainsBodyOrigin::Clone() const
    {
      return new ContainsBodyOrigin( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ContainsBodyOrigin::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator() determines if the BODY is occupied by any of the geometries of the SSE
    //! @details this is done by checking if any of the fragments in the body will contain the center of any of the
    //!          fragments in the sse
    //! @param BODY representing the density rod
    //! @param SSE represents the secondary structure element which contains geometries/fragments
    //! @return returns true if BODY contains any of the SSE's geometries
    bool ContainsBodyOrigin::operator()( const assemble::SSEGeometryInterface &BODY, const assemble::SSE &SSE) const
    {
      // create local copy of fragments of body (this is necessary because GetGeometries() returns only local copy of fragments)
      const util::SiPtrVector< const assemble::SSEGeometryInterface> body_geometries( BODY.GetSSEGeometries());

      // create local copy of fragments of sse
      const util::SiPtrVector< const assemble::SSEGeometryInterface> sse_geometries( SSE.GetSSEGeometries());

      // iterate over the fragments of BODY
      for
      (
        util::SiPtrVector< const assemble::SSEGeometryInterface>::const_iterator
          body_itr( body_geometries.Begin()), body_itr_end( body_geometries.End());
        body_itr != body_itr_end; ++body_itr
      )
      {
        // iterate over sse fragments
        for
        (
          util::SiPtrVector< const assemble::SSEGeometryInterface>::const_iterator
            sse_itr( sse_geometries.Begin()), sse_itr_end( sse_geometries.End());
          sse_itr != sse_itr_end; ++sse_itr
        )
        {
          // check whether origin of sse fragment is within the body_fragment
          if( ContainsBodyOriginCheck( **body_itr, ( *sse_itr)->GetCenter()))
          {
            return true;
          }
        }
      }

      // non of the sse fragments was found to be within one of the body fragments
      return false;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read distance from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ContainsBodyOrigin::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write distance to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &ContainsBodyOrigin::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief ContainsBodyOriginCheck determines if the coord::GeometryInterface envelops the origin
    //! @param BODY coord::GeometryInterface which is the body involved
    //! @param ORIGIN origin for which is checked whether it is within body
    //! @return returns true if Body envelops the origin
    bool ContainsBodyOrigin::ContainsBodyOriginCheck( const assemble::SSEGeometryInterface &BODY, const linal::Vector3D &ORIGIN) const
    {
      // calculate the distance of the origin from the z-axis of the body
      const double dist_origin_from_zaxis
      (
        linal::CalculateDistancePointFromLine( ORIGIN, BODY.GetCenter(), BODY.GetAxis( coord::GetAxes().e_Z))
      );

      // calculate the footpoint of the origin on the z-axis of the body
      const linal::Vector3D footpoint_origin_on_zaxis
      (
        linal::CalculateFootpoint( ORIGIN, BODY.GetCenter(), BODY.GetAxis( coord::GetAxes().e_Z))
      );

      // return true if distance of origin to footpoint is less than the sse radius and if distance of footpoint
      // to center of body is less than the z extent
      return
      (
        dist_origin_from_zaxis <= BODY.GetExtent( coord::GetAxes().e_X) &&
        linal::Distance( footpoint_origin_on_zaxis, BODY.GetCenter()) <= BODY.GetExtent( coord::GetAxes().e_Z)
      );
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "restraint/bcl_restraint_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

    //! @brief return command line flag for specifying desired restraint types to use
    //! @return command line flag for specifying desired restraint types to use
    const util::ShPtr< command::FlagInterface> &GetFlagRestraintsTypes()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagDynamic
        (
          "restraint_types",
          "\tone or more restraint types to be used",
          command::Parameter
          (
            "restraint_type",
            "any restraint type from the list",
            command::ParameterCheckSerializable( Type())
          ),
          0,
          util::Enumerated< Interface>::GetSize()
        )
      );

      return s_flag;
    }

    //! @brief return command line flag for specifying the prefix prepended to each restraint's postfix
    //! @return command line flag for specifying the prefix prepended to each restraint's postfix
    util::ShPtr< command::FlagInterface> &GetFlagRestraintsFilePrefix()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "restraint_prefix", "\tthe prefix prepended to each restraint's postfix",
          command::Parameter( "prefix", "\tthe prefix prepended to each restraint's postfix", "")
        )
      );
      // end
      return s_flag;
    }

    //! @brief identifier for the name space
    //! @return the name of the namespace
    const std::string &GetNamespaceIdentifier()
    {
      static const std::string *s_namespace_name( new std::string( util::ExtractNamespaceIdentifier( __PRETTY_FUNCTION__)));
      return *s_namespace_name;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_data_pairwise.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_atom.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "storage/bcl_storage_triplet.h"
#include "util/bcl_util_string_replacement.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> DataPairwise::s_Instance
    (
      GetObjectInstances().AddInstance( new DataPairwise())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    DataPairwise::DataPairwise() :
      m_First(),
      m_Second()
    {
    }

    //! @brief constructor taking members
    //! @param LOCATOR_A first locator
    //! @param LOCATOR_B second locator
    DataPairwise::DataPairwise
    (
      const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &LOCATOR_A,
      const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &LOCATOR_B
    ) :
      m_First(),
      m_Second()
    {
      Set( LOCATOR_A, LOCATOR_B);
    }

    //! @brief Clone function
    //! @return pointer to new DataPairwise
    DataPairwise *DataPairwise::Clone() const
    {
      return new DataPairwise( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DataPairwise::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives formatted string describing the data pair
    //! @return formatted string describing the data pair
    std::string DataPairwise::GetIdentification() const
    {
      // true if either of the locators is not defined
      if( !IsSet())
      {
        return std::string();
      }

      std::stringstream write;
      write << First()->GetIdentification();
      io::Serialize::Write( std::string( "<=>"), write, 1);
      write << Second()->GetIdentification();
      // return identification string
      return write.str();
    }

    //! @brief reads formatted string describing the locator
    //! @return formatted string describing the locator
    std::istream &DataPairwise::ReadIdentification( std::istream &ISTREAM)
    {
      m_First->ReadIdentification( ISTREAM);
      BCL_MessageDbg( "read in " + m_First->GetIdentification());
      std::string separator;
      io::Serialize::Read( separator, ISTREAM);
      BCL_MessageDbg( "read in " + separator);
      m_Second->ReadIdentification( ISTREAM);
      BCL_MessageDbg( "read in " + m_Second->GetIdentification());
      return ISTREAM;
    }

    //! @brief gives the first point in the data pair
    //! @return the first point in the data pair
    const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &DataPairwise::First() const
    {
      // true if this is not set
      if( !IsSet())
      {
        static const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> s_undefined;
        return s_undefined;
      }

      // return the first locator
      return m_First;
    }

    //! @brief gives the second point in the data pair
    //! @return the second point in the data pair
    const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &DataPairwise::Second() const
    {
      // true if this is not set
      if( !IsSet())
      {
        static util::ShPtr< assemble::LocatorAtomCoordinatesInterface> s_undefined;
        return s_undefined;
      }

      // return the second locator
      return m_Second;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief sets the two data points involved in this data pair
    //! @param ATOM_A the  first data point that will be set
    //! @param ATOM_B the second data point that will be set
    void DataPairwise::Set
    (
      const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &ATOM_A,
      const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &ATOM_B
    )
    {
      if( assemble::LocatorAtomCoordinatesInterface::PtrLessThan()( ATOM_A, ATOM_B))
      {
        m_First = ATOM_A;
        m_Second = ATOM_B;
      }
      else
      {
        m_First = ATOM_B;
        m_Second = ATOM_A;
      }
    }

    //! @brief indicates whether or not the two data points have been set
    //! @return boolean true if the data pair has been successfully set - false otherwise
    bool DataPairwise::IsSet() const
    {
      return m_First.IsDefined() && m_Second.IsDefined();
    }

    //! @brief calculates the euclidian distance indicated by a protein model
    //! @param MODEL the model from which the distance will be calculated
    //! @return double which is the distance between the two coordinates indicated by this data pair
    double DataPairwise::EuclidianDistance( const assemble::ProteinModel &MODEL) const
    {
      // locate the coordinates for the two locators
      const linal::Vector3D coords_a( m_First->Locate( MODEL));
      const linal::Vector3D coords_b( m_Second->Locate( MODEL));

      // true if either of the coordinates are not defined
      if( !coords_a.IsDefined() || !coords_b.IsDefined())
      {
        // return undefined double
        return util::GetUndefinedDouble();
      }

      // holds the distance between the two coordinate points
      const double distance( linal::Distance( coords_a, coords_b));

      // return the calculated distance
      return distance;
    }

    //! @brief calculates statistics for the distance indicated by this across an ensemble of models
    //! @param ENSEMBLE the ensemble over which distances and statistics will be calculated
    //! @return pair of RunningAverageSD< double> and RunningMinMax< double> indicated distance statistics over ENSEMBLE
    storage::Pair< math::RunningAverageSD< double>, math::RunningMinMax< double> >
    DataPairwise::EuclidianDistance( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // will keep track of the min max mean and stddev statistics
      storage::Pair< math::RunningAverageSD< double>, math::RunningMinMax< double> > stats;

      // iterate through the models of the ensemble
      for
      (
        assemble::ProteinEnsemble::const_iterator ensemble_itr( ENSEMBLE.Begin()), ensemble_itr_end( ENSEMBLE.End());
        ensemble_itr != ensemble_itr_end;
        ++ensemble_itr
      )
      {
        // locate the coordinates for the two locators
        linal::Vector3D coords_a( m_First->Locate( **ensemble_itr));
        linal::Vector3D coords_b( m_Second->Locate( **ensemble_itr));

        // holds the distance between the two coordinate points
        const double distance( linal::Distance( coords_a, coords_b));

        // true if either of the coordinates is not defined or the distance is not defined
        if( !coords_a.IsDefined() || !coords_b.IsDefined() || !util::IsDefined( distance))
        {
          // go to next model in ensemble
          continue;
        }

        // add the distance two the two statistic objects
        stats.First() += distance;
        stats.Second() += distance;
      }

      // return the statistics objects
      return stats;
    }

    //! @brief calculates the sequence separation between the two data points in this data pair
    //! @return size_t which is the sequence separation between the two data points in this data pair
    size_t DataPairwise::SequenceSeparation() const
    {
      return CalculateSequenceSeparation
      (
        First()->GetChainID(), First()->GetSeqID(), Second()->GetChainID(), Second()->GetSeqID()
      );
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
    std::istream &DataPairwise::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_First, ISTREAM);
      io::Serialize::Read( m_Second, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DataPairwise::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_First, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Second, OSTREAM, INDENT) << '\n';

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculates sequence separation between two residues indicated by chain and seq ids
    //! @param CHAIN_ID_A the chain id of the first residue
    //! @param SEQ_ID_A the seq id of the first residue
    //! @param CHAIN_ID_B the chain id of the second residue
    //! @param SEQ_ID_B the seq id of the second residue
    size_t DataPairwise::CalculateSequenceSeparation
    (
      const char CHAIN_ID_A, const int SEQ_ID_A, const char CHAIN_ID_B, const int SEQ_ID_B
    )
    {
      if( CHAIN_ID_A != CHAIN_ID_B)
      {
        return util::GetUndefinedSize_t();
      }

      return math::Absolute( SEQ_ID_A - SEQ_ID_B);
    }

    //! @brief less than operator for comparing two DataPairwise
    //! @param LHS the first DataPairwise which will be compared against the second DataPairwise
    //! @param RHS the second DataPairwise which will be compared against the first DataPairwise
    //! @return boolean true if LHS is less than RHS - false otherwise
    bool operator <( const DataPairwise &LHS, const DataPairwise &RHS)
    {
      // true if LHS is not set - cannot be less than RHS
      if( !LHS.IsSet())
      {
        BCL_MessageDbg( "LHS " + LHS.GetIdentification() + " RHS " + RHS.GetIdentification());
        return false;
      }

      // true if RHS is not set but LHS is set - LHS is less than RHS
      if( LHS.IsSet() && !RHS.IsSet())
      {
        BCL_MessageDbg( "LHS " + LHS.GetIdentification() + " RHS " + RHS.GetIdentification());
        return true;
      }

      // make vector nd 2 of locator atom from the LHS argument
      const storage::VectorND< 2, assemble::LocatorAtom> lhs
      (
        assemble::LocatorAtom( LHS.First()->GetChainID(), LHS.First()->GetSeqID(), LHS.First()->GetAtomType()),
        assemble::LocatorAtom( LHS.Second()->GetChainID(), LHS.Second()->GetSeqID(), LHS.Second()->GetAtomType())
      );

      // make vector nd 2 of locator atom from the RHS argument
      const storage::VectorND< 2, assemble::LocatorAtom> rhs
      (
        assemble::LocatorAtom( RHS.First()->GetChainID(), RHS.First()->GetSeqID(), RHS.First()->GetAtomType()),
        assemble::LocatorAtom( RHS.Second()->GetChainID(), RHS.Second()->GetSeqID(), RHS.Second()->GetAtomType())
      );

      BCL_MessageDbg( "LHS " + LHS.GetIdentification() + " RHS " + RHS.GetIdentification());
      // use the less than operator of the locator atom and the vector nd
      return lhs < rhs;
    }

    //! @brief writes pymol script formatted data to a stream that can show distances in pymol
    //! @param OSTREAM the stream which will write the data
    //! @param DISTANCES the data that will be written along with their distance and the desired color of their line
    //! @return ostream the stream that wrote the data distance information
    std::ostream &ShowDistancesInPymol
    (
      std::ostream &OSTREAM,
      const storage::List< storage::Triplet< DataPairwise, double, linal::Vector3D> > &DISTANCES
    )
    {
      // iterate through the distances
      for
      (
        storage::List< storage::Triplet< DataPairwise, double, linal::Vector3D> >::const_iterator
          distance_itr( DISTANCES.Begin()), distance_itr_end( DISTANCES.End());
        distance_itr != distance_itr_end;
        ++distance_itr
      )
      {
        // get critical data
        const std::string first_pymol_name( distance_itr->First().First()->GetPymolName());
        const std::string second_pymol_name( distance_itr->First().Second()->GetPymolName());
        const std::string dist_name( first_pymol_name + "_" + second_pymol_name);
        const std::string color_name( "color_" + dist_name);
        const char chain_a( distance_itr->First().First()->GetChainID());
        const int  resi_a(  distance_itr->First().First()->GetSeqID());
        const char chain_b( distance_itr->First().Second()->GetChainID());
        const int  resi_b(  distance_itr->First().Second()->GetSeqID());
        std::string atom_a( distance_itr->First().First()->GetAtomType());
        std::string atom_b( distance_itr->First().Second()->GetAtomType());

        // no hydrogens so replace with CA for visualization purposes
        const util::StringReplacement replacer
        (
          util::StringReplacement::e_Any, biol::GetAtomTypes().HA2, biol::GetAtomTypes().CA
        );
        replacer.ReplaceEachIn( atom_a);
        replacer.ReplaceEachIn( atom_b);

        // write the distance selection
        OSTREAM << "distance " << dist_name << ", (chain " << chain_a << " and resi " << resi_a << " and name " << atom_a
        << "),(chain " << chain_b << " and resi " << resi_b << " and name " << atom_b << ")" << '\n';

        // define the color for the current line
        const linal::Vector3D color( distance_itr->Third());
        OSTREAM << "set_color " << color_name << ", [" << color.X() << "," << color.Y() << "," << color.Z() << "]"
        << '\n';

        // set the color of the dashes
        OSTREAM << " set dash_color, " << color_name << ", " << dist_name << '\n';

        // set the gap between dashes
        static const double s_dash_gap( 0);
        OSTREAM << " set dash_gap, " << s_dash_gap << ", " << dist_name << '\n';

        // set the dash width
        static const double s_dash_width( 5);
        OSTREAM << " set dash_width, " << s_dash_width << ", " << dist_name << '\n';
      }

      return OSTREAM;
    }

  } // namespace restraint

} // namespace bcl
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
#include "restraint/bcl_restraint_data_set_pairwise.h"

// includes from bcl - sorted alphabetically
#include "restraint/bcl_restraint_locator_coordinates_first_side_chain_atom.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> DataSetPairwise::s_Instance
    (
      GetObjectInstances().AddInstance( new DataSetPairwise())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    DataSetPairwise::DataSetPairwise() :
      m_DataSet()
    {
    }

    //! @brief Clone function
    //! @return pointer to new Sequence
    DataSetPairwise *DataSetPairwise::Clone() const
    {
      return new DataSetPairwise( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DataSetPairwise::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return iterator on begin
    //! @return iterator pointing to the beginning of the container, i.e. the first element
    DataSetPairwise::iterator DataSetPairwise::Begin()
    {
      return m_DataSet.Begin();
    }

    //! @brief return const_iterator on begin
    //! @return const_iterator pointing to the beginning of the container, i.e. the first element
    DataSetPairwise::const_iterator DataSetPairwise::Begin() const
    {
      return m_DataSet.Begin();
    }

    //! @brief return iterator on end
    //! @return iterator pointing to the end of the container, i.e. behind the last element
    DataSetPairwise::iterator DataSetPairwise::End()
    {
      return m_DataSet.End();
    }

    //! @brief return const_iterator on end
    //! @return const_iterator pointing to the end of the container, i.e. behind the last element
    DataSetPairwise::const_iterator DataSetPairwise::End() const
    {
      return m_DataSet.End();
    }

    //! @brief return iterator to reverse begin
    //! @return reverse_iterator pointing to the beginning of the reversed container, i.e. the last element
    DataSetPairwise::reverse_iterator DataSetPairwise::ReverseBegin()
    {
      return m_DataSet.ReverseBegin();
    }

    //! @brief return const_iterator to reverse begin
    //! @return const_reverse_iterator pointing to the beginning of the reversed container
    DataSetPairwise::const_reverse_iterator DataSetPairwise::ReverseBegin() const
    {
      return m_DataSet.ReverseBegin();
    }

    //! @brief return iterator to reverse end
    //! @return reverse_iterator pointing to the end of the reversed container, i.e. behind the first element
    DataSetPairwise::reverse_iterator DataSetPairwise::ReverseEnd()
    {
      return m_DataSet.ReverseEnd();
    }

    //! @brief return const_iterator to reverse end
    //! @return const_reverse_iterator pointing to the end of the reversed container
    DataSetPairwise::const_reverse_iterator DataSetPairwise::ReverseEnd() const
    {
      return m_DataSet.ReverseEnd();
    }

    //! @brief indicates if there is nothing in the data set or not
    //! @return bool true if dataset is empty -false otherwise
    bool DataSetPairwise::IsEmpty() const
    {
      return m_DataSet.IsEmpty();
    }

    //! @brief gives the number of elements in the data set
    //! @return size_t indicating the number of elements in the data set
    size_t DataSetPairwise::GetSize() const
    {
      return m_DataSet.GetSize();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief AddPair function for adding a pair of atoms to the dataset
    //! @param LOCATOR_A first of the two atoms to be added to the dataset
    //! @param LOCATOR_B second of the atoms to be paired with the first
    //! @return pair of bool true if insertion was successful - false otherwise and iterator pointing to insertion place
    std::pair< DataSetPairwise::iterator, bool> DataSetPairwise::Insert
    (
      const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &LOCATOR_A,
      const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &LOCATOR_B
    )
    {
      // make an empty pairwise data
      DataPairwise data;

      // set the data
      data.Set( LOCATOR_A, LOCATOR_B);

      // try to insert the data and return the result
      return std::pair< iterator, bool>( m_DataSet.Insert( data));
    }

    //! @brief AddPair function for adding a pair of atoms to the dataset
    //! @param DATA data pair to be added
    //! @return pair of bool  true if insertion was successful - false otherwise and iterator pointing to place of insertion
    std::pair< DataSetPairwise::iterator, bool> DataSetPairwise::Insert
    (
      const DataPairwise &DATA
    )
    {
      return std::pair< iterator, bool>( m_DataSet.Insert( DATA));
    }

    //! @brief Erase function for removing a data pair from the data set
    //! @param DATA_PAIR the data points to be removed
    //! @return bool true if DATA_PAIR was removed - false otherwise
    bool DataSetPairwise::Erase( const DataPairwise &DATA_PAIR)
    {
      return m_DataSet.Erase( DATA_PAIR);
    }

    //! @brief gives all of the single points that are unique in the dataset
    //! @return set which has the individual data points that are unique from the dataset
    storage::Set
    <
      util::ShPtr< assemble::LocatorAtomCoordinatesInterface>,
      assemble::LocatorAtomCoordinatesInterface::PtrLessThan
    > DataSetPairwise::GetUniqueDataPoints() const
    {
      // get the list of all data points
      storage::Set
      <
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, assemble::LocatorAtomCoordinatesInterface::PtrLessThan
      > data_points;

      // iterate through the data set to fill up all the data points
      for
      (
        DataSetPairwise::const_iterator data_itr( Begin()), data_itr_end( End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // try to insert the data points, only unique points will be able to be inserted
        data_points.Insert( data_itr->First());
        data_points.Insert( data_itr->Second());
      }

      return data_points;
    }

    //! @brief gives all of the single points in the dataset - leaves in duplicates
    //! @return set which has the individual data points leaving in duplicates
    util::ShPtrList< assemble::LocatorAtomCoordinatesInterface> DataSetPairwise::GetDataPoints() const
    {
      util::ShPtrList< assemble::LocatorAtomCoordinatesInterface> data_points;

      // iterate through the data set to fill up all the data points
      for
      (
        DataSetPairwise::const_iterator data_itr( Begin()), data_itr_end( End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // try to insert the data points
        data_points.PushBack( data_itr->First());
        data_points.PushBack( data_itr->Second());
      }

      // sort the data points by sequence
      data_points.Sort( assemble::LocatorAtomCoordinatesInterface::PtrLessThan());

      return data_points;
    }

    //! @brief Find takes a key and returns an iterator to the element denoted by key, or an off the end iterator
    //! @param KEY the key for which an iterator is desired
    //! @return returns an iterator to the element
    DataSetPairwise::iterator DataSetPairwise::Find( const DataPairwise &KEY)
    {
      return m_DataSet.Find( KEY);
    }

    //! @brief Find takes a key and returns a const_iterator to the element denoted by key, or off the end iterator
    //! @param KEY the key for which a const_iterator is desired
    //! @return returns an iterator to the element
    DataSetPairwise::const_iterator DataSetPairwise::Find( const DataPairwise &KEY) const
    {
      return m_DataSet.Find( KEY);
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
    std::istream &DataSetPairwise::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_DataSet, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DataSetPairwise::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_DataSet, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief creates a complete pairwise dataset from vector of aabases
    //! @param SEQUENCE the sequence for which a complete data set will be created
    //! @return DataSetPairwise which has the complete dataset for SEQUENCE
    util::ShPtr< DataSetPairwise> DataSetPairwise::GetCompleteDataSet( const util::ShPtrVector< biol::AABase> &SEQUENCE)
    {
      util::ShPtr< DataSetPairwise> dataset( new DataSetPairwise());

      // iterate through the sequence
      for
      (
        biol::AASequence::const_iterator aa_itr_a( SEQUENCE.Begin()), aa_itr_end( SEQUENCE.End());
        aa_itr_a != aa_itr_end;
        ++aa_itr_a
      )
      {
        // make locator a
        const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_a
        (
          assemble::LocatorAtomCoordinatesInterface::CreateLocator< LocatorCoordinatesFirstSideChainAtom>
          (
            **aa_itr_a
          )
        );

        // iterate through the sequence again
        for
        (
          biol::AASequence::const_iterator aa_itr_b( ++biol::AASequence::const_iterator( aa_itr_a));
          aa_itr_b != aa_itr_end;
          ++aa_itr_b
        )
        {
          // make locator b
          const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_b
          (
            assemble::LocatorAtomCoordinatesInterface::CreateLocator< LocatorCoordinatesFirstSideChainAtom>
            (
              **aa_itr_b
            )
          );

          // try to insert locator a and b into the dataset
          const bool inserted( dataset->Insert( locator_a, locator_b).second);

          // make sure the pair could be inserted
          BCL_Assert
          (
            inserted, "could not insert " + locator_a->GetIdentification() + " and " + locator_b->GetIdentification()
          );
        }
      }

      BCL_MessageDbg( "complete dataset size is " + util::Format()( dataset->GetSize()));

      // return the data set
      return dataset;
    }

    //! @brief gets the intersection out of a list of data sets
    //! @param DATA_SETS the data sets whose intersection will be determined
    //! @return a data set that is the intersection of all that have been passed
    DataSetPairwise DataSetPairwise::GetIntersection( const storage::List< DataSetPairwise> &DATA_SETS)
    {
      BCL_Assert( !DATA_SETS.IsEmpty(), "DATA_SETS is empty");

      // get first data set
      DataSetPairwise start_data_set( DATA_SETS.FirstElement());

      // will hold the current intersection
      std::set< DataPairwise> intersection;

      // iterate through the data sets
      for
      (
        storage::List< DataSetPairwise>::const_iterator
          set_itr( ++DATA_SETS.Begin()), set_itr_end( DATA_SETS.End());
        set_itr != set_itr_end; ++set_itr
      )
      {
        // the current set
        const DataSetPairwise &current_set( *set_itr);

        // get the current intersection
        std::set_intersection
        (
          start_data_set.Begin(),      start_data_set.End(),
          current_set.Begin(),         current_set.End(),
          std::inserter( intersection, intersection.begin())
        );

        BCL_MessageDbg( "intersection size is " + util::Format()( intersection.size()));

        // set to the current intersection
        start_data_set = DataSetPairwise( intersection.begin(), intersection.end());

        // clear the current intersection
        intersection.clear();
      }

      // return the intersection of all data sets
      return start_data_set;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator * to multiply mutate result of dataset pairwise with a scalar
    //!        repeats the list of mutates nodes on the argument the number of times indicated by the scalar
    //! @param SCALAR scalar to multiply mutate result by
    //! @param MUTATE_RESULT result that will be multiplied
    //! @return new mutate result that has been multiplied by the scalar
    math::MutateResult< DataSetPairwise>
    operator *
    (
      const double &SCALAR,
      const math::MutateResult< DataSetPairwise> &MUTATE_RESULT
    )
    {
      // get the list of mutates
      const util::SiPtrList< const math::MutateInterface< DataSetPairwise> > &mutates( MUTATE_RESULT.GetNodes());

      util::SiPtrList< const math::MutateInterface< DataSetPairwise> > all_mutates;

      util::ShPtr< DataSetPairwise> mutated_dataset( MUTATE_RESULT.GetArgument());
      BCL_Assert( mutated_dataset.IsDefined(), "mutated_dataset pointer is not defined");

      // repeat the list of mutates SCALAR number of times
      for( double repeat( 0); repeat < SCALAR - 1; ++repeat)
      {
        // iterate through the mutates
        for
        (
          util::SiPtrList< const math::MutateInterface< DataSetPairwise> >::const_iterator
            mutate_itr( mutates.Begin()), mutate_itr_end( mutates.End());
          mutate_itr != mutate_itr_end;
          ++mutate_itr
        )
        {
          // apply the current mutate
          const math::MutateResult< DataSetPairwise> current_mutate_result( ( *mutate_itr)->operator()( *mutated_dataset));

          mutated_dataset = current_mutate_result.GetArgument();
          BCL_Assert( mutated_dataset.IsDefined(), "mutated_dataset pointer is not defined");

          // add mutates to all_mutates
          all_mutates.Prepend( current_mutate_result.GetNodes());
        }
      }

      // make the final mutate result with final mutated dataset and all the mutates that led to it
      const math::MutateResult< DataSetPairwise> mutate_result( mutated_dataset, all_mutates);
      BCL_Assert( mutated_dataset.IsDefined(), "mutated_dataset pointer is not defined");

      // add the
      return mutate_result;
    }

    //! @brief operator * to multiply mutate result of dataset pairwise with a scalar
    //!        repeats the list of mutates nodes on the argument the number of times indicated by the scalar
    //! @param SCALAR scalar to multiply mutate result by
    //! @param MUTATE_RESULT result that will be multiplied
    //! @return new mutate result that has been multiplied by the scalar
    math::MutateResult< DataSetPairwise>
    operator *
    (
      const math::MutateResult< DataSetPairwise> &MUTATE_RESULT,
      const double &SCALAR
    )
    {
      return operator *( SCALAR, MUTATE_RESULT);
    }

    //! @brief operator + to add two mutate results of dataset pairwise
    //!        takes the intersection of the two mutate results and combines the mutates
    //! @param LHS first argument
    //! @param RHS second argument
    //! @return new mutate result that is the sum of LHS and RHS
    math::MutateResult< DataSetPairwise>
    operator +
    (
      const math::MutateResult< DataSetPairwise> &LHS,
      const math::MutateResult< DataSetPairwise> &RHS
    )
    {
      BCL_MessageStd( "plus");
      util::ShPtr< DataSetPairwise> intersection;
      if( LHS.GetArgument().IsDefined() && RHS.GetArgument().IsDefined())
      {
        BCL_MessageStd( "DataSetPairwise::GetIntersection");
        storage::List< DataSetPairwise> results( 1, *LHS.GetArgument());
        results.PushBack( *RHS.GetArgument());
        intersection = util::ShPtr< DataSetPairwise>( DataSetPairwise::GetIntersection( results).Clone());
      }
      else if( LHS.GetArgument().IsDefined())
      {
        BCL_MessageStd( "LHS.GetArgument().IsDefined()");
        intersection = LHS.GetArgument();
      }
      else if( RHS.GetArgument().IsDefined())
      {
        BCL_MessageStd( "RHS.GetArgument().IsDefined()");
        intersection = RHS.GetArgument();
      }
      util::SiPtrList< const math::MutateInterface< DataSetPairwise> > nodes( LHS.GetNodes());
      nodes.Prepend( RHS.GetNodes());
      BCL_Assert( intersection.IsDefined(), "data set pointer is not defined");
      const math::MutateResult< DataSetPairwise> result( intersection, nodes);

      return result;
    }

  } // namespace restraint

} // namespace bcl
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
#include "restraint/bcl_restraint_distance.h"

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Distance::s_Instance
    (
      GetObjectInstances().AddInstance( new Distance())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Distance::Distance() :
      m_Distance( util::GetUndefined< double>()),
      m_UpperBound( util::GetUndefined< double>()),
      m_LowerBound( util::GetUndefined< double>())
    {
    }

    //! @brief construct from a distance and a margin of error
    //! @param DISTANCE double which is the distance
    //! @param UPPER_BOUND is the upper bound possible for the distance
    //! @param LOWER_BOUND is the lower bound possible for the distance
    Distance::Distance( const double &DISTANCE, const double &UPPER_BOUND, const double &LOWER_BOUND) :
      m_Distance( DISTANCE),
      m_UpperBound( UPPER_BOUND),
      m_LowerBound( LOWER_BOUND)
    {
      BCL_Assert
      (
        m_UpperBound >= m_Distance && m_LowerBound <= m_Distance,
        "The lower bound " + util::Format()( LOWER_BOUND) + " should be lower than the distance " +
        util::Format()( DISTANCE) + " and the upper bound " + util::Format()( UPPER_BOUND) + " should be higher than the distance"
      );
    }

    //! @brief virtual copy constructor
    Distance *Distance::Clone() const
    {
      return new Distance( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Distance::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives formatted string describing this
    //! @return formatted string describing this
    const std::string Distance::GetIdentification() const
    {
      return util::Format()( LowerBound()) + " " +
        util::Format()( GetDistance()) + " " + util::Format()( UpperBound());
    }

    //! @brief GetDistance return m_Distance
    //! @return returns double which is m_Distance
    double Distance::GetDistance() const
    {
      return m_Distance;
    }

    //! @brief UpperBound returns m_UpperBound
    //! @return returns double which is m_UpperBound
    double Distance::UpperBound() const
    {
      return m_UpperBound;
    }

    //! @brief LowerBound returns m_LowerBound
    //! @return returns double which is m_LowerBound
    double Distance::LowerBound() const
    {
      return m_LowerBound;
    }

    //! @brief determines if the distance object is defined or not
    //! @return bool true if distance, upper bound, and lower bound are all defined
    bool Distance::IsDefined() const
    {
      return util::IsDefined( m_Distance) && util::IsDefined( m_LowerBound) && util::IsDefined( m_UpperBound);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read distance from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Distance::Read( std::istream &ISTREAM)
    {
      // write members
      io::Serialize::Read( m_Distance, ISTREAM);
      io::Serialize::Read( m_UpperBound, ISTREAM);
      io::Serialize::Read( m_LowerBound, ISTREAM);

      return ISTREAM;
    }

    //! @brief write distance to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &Distance::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Distance, OSTREAM, INDENT);
      io::Serialize::Write( m_UpperBound, OSTREAM, INDENT);
      io::Serialize::Write( m_LowerBound, OSTREAM, INDENT);

      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculates the upper error
    //! @return the upper error
    double Distance::GetUpperError() const
    {
      return m_UpperBound - m_Distance;
    }

    //! @brief calculates the lower error
    //! @return the lower error
    double Distance::GetLowerError() const
    {
      return m_Distance - m_LowerBound;
    }

  } // namespace restraint
} // namespace bcl

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
#include "restraint/bcl_restraint_epr_accessibility_data.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_vector.h"
#include "command/bcl_command_parameter.h"
#include "fold/bcl_fold_scores.h"
#include "pdb/bcl_pdb_factory.h"
#include "score/bcl_score_accessibility.h"
#include "score/bcl_score_accessibility_hydrophobic_moment_magnitude.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    // initialize scores
    fold::Score EPRAccessibilityData::e_ScoreExposureMomentNiEDDA( fold::GetScores().e_Undefined);
    fold::Score EPRAccessibilityData::e_ScoreExposureMomentOxygen( fold::GetScores().e_Undefined);
    fold::Score EPRAccessibilityData::e_ScoreExposureMomentMagnitudeNiEDDA( fold::GetScores().e_Undefined);
    fold::Score EPRAccessibilityData::e_ScoreExposureMomentMagnitudeOxygen( fold::GetScores().e_Undefined);

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> EPRAccessibilityData::s_Instance
    (
      util::Enumerated< Interface>::AddInstance( new EPRAccessibilityData())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor taking member variable parameters
    //! @param HANDLER the handler that will be used to read and write the restraints
    EPRAccessibilityData::EPRAccessibilityData
    (
      const HandlerAccessibilityAA &HANDLER
    ) :
      m_Handler( HANDLER),
      m_HydrophobicMomentWindowSizes( size_t( 5), size_t( 5), size_t( 5)),
      m_Restraints()
    {
    }

    //! @brief Clone function
    //! @return pointer to new EPRAccessibilityData
    EPRAccessibilityData *EPRAccessibilityData::Clone() const
    {
      return new EPRAccessibilityData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &EPRAccessibilityData::GetAlias() const
    {
      static const std::string s_name( "AccessibilityEPR");
      return s_name;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &EPRAccessibilityData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the default file extension
    //! @return the default file extension
    const std::string &EPRAccessibilityData::GetDefaultExtension() const
    {
      static const std::string s_extension( ".access_bcl");
      return s_extension;
    }

    //! @brief get the default restraint handler
    //! @return the default restraint handler
    const HandlerAccessibilityAA &EPRAccessibilityData::GetDefaultHandler()
    {
      static const HandlerAccessibilityAA s_handler
      (
        util::ShPtr< assemble::AAExposureInterface>( new assemble::AANeighborVector())
      );
      return s_handler;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initialize the scores and add them to Scores enumerator
    void EPRAccessibilityData::InitializeScores()
    {
      // read restraints from the restraint file if necessary
      if( !m_Restraints.IsDefined() && m_Handler.Exists())
      {
        m_Restraints = util::CloneToShPtr( m_Handler.ReadRestraintsFromFile());
      }
      const size_t &helix_window( m_HydrophobicMomentWindowSizes( biol::GetSSTypes().HELIX));
      const size_t &strand_window( m_HydrophobicMomentWindowSizes( biol::GetSSTypes().STRAND));
      const size_t &coil_window( m_HydrophobicMomentWindowSizes( biol::GetSSTypes().COIL));

      if( !e_ScoreExposureMomentNiEDDA.IsDefined())
      {
        const AccessibilityAA::EnvironmentEnum environment( AccessibilityAA::e_NiEDDA);
        const util::ShPtr< score::AccessibilityHydrophobicMoment> score
        (
          new score::AccessibilityHydrophobicMoment
          (
            environment,
            pdb::Factory::GetSSETypeMinSizes( helix_window, strand_window, coil_window)
          )
        );
        e_ScoreExposureMomentNiEDDA = fold::GetScores().AddScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::Accessibility
            (
              score,
              m_Restraints,
              score::Accessibility::GetDefaultScheme() + environment.GetString()
            )
          )
        );
      }
      if( !e_ScoreExposureMomentOxygen.IsDefined())
      {
        const AccessibilityAA::EnvironmentEnum environment( AccessibilityAA::e_Oxygen);
        const util::ShPtr< score::AccessibilityHydrophobicMoment> score
        (
          new score::AccessibilityHydrophobicMoment
          (
            environment,
            pdb::Factory::GetSSETypeMinSizes( helix_window, strand_window, coil_window)
          )
        );
        e_ScoreExposureMomentOxygen = fold::GetScores().AddScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::Accessibility
            (
              score, m_Restraints, score::Accessibility::GetDefaultScheme() + environment.GetString()
            )
          )
        );
      }
      if( !e_ScoreExposureMomentMagnitudeNiEDDA.IsDefined())
      {
        const AccessibilityAA::EnvironmentEnum environment( AccessibilityAA::e_NiEDDA);
        const util::ShPtr< score::AccessibilityHydrophobicMoment> score
        (
          new score::AccessibilityHydrophobicMomentMagnitude
          (
            environment,
            pdb::Factory::GetSSETypeMinSizes( helix_window, strand_window, coil_window)
          )
        );
        e_ScoreExposureMomentMagnitudeNiEDDA = fold::GetScores().AddScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::Accessibility
            (
              score,
              m_Restraints,
              score::Accessibility::GetDefaultScheme() + "_magn_" + environment.GetString()
            )
          )
        );
      }
      if( !e_ScoreExposureMomentMagnitudeOxygen.IsDefined())
      {
        const AccessibilityAA::EnvironmentEnum environment( AccessibilityAA::e_Oxygen);
        const util::ShPtr< score::AccessibilityHydrophobicMoment> score
        (
          new score::AccessibilityHydrophobicMomentMagnitude
          (
            environment,
            pdb::Factory::GetSSETypeMinSizes( helix_window, strand_window, coil_window)
          )
        );
        e_ScoreExposureMomentMagnitudeOxygen = fold::GetScores().AddScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::Accessibility
            (
              score, m_Restraints, score::Accessibility::GetDefaultScheme() + "_magn_" + environment.GetString()
            )
          )
        );
      }
    }

    //! @brief sets the weights of scores in a weight set
    //! @param SCORE_WEIGHT_SET the score weight set that will be modified
    void EPRAccessibilityData::ModifyScoreWeightSet( fold::ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void EPRAccessibilityData::MergeAndModifyMutateTree( fold::MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< fold::MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer EPRAccessibilityData::GetSerializer() const
    {
      io::Serializer serial( m_Handler.GetSerializer());
      serial.SetClassDescription( "Allows use of electron paramagnetic resonance-based accessibility data measurements");

      serial.AddInitializer
      (
        "helix window size",
        "Number of residues considered when computing hydrophobic moment in helices",
        io::Serialization::GetAgent( &m_HydrophobicMomentWindowSizes( 0)),
        "5"
      );

      serial.AddInitializer
      (
        "strand window size",
        "Number of residues considered when computing hydrophobic moment in strands",
        io::Serialization::GetAgent( &m_HydrophobicMomentWindowSizes( 1)),
        "5"
      );
      serial.AddInitializer
      (
        "loop window size",
        "Number of residues considered when computing hydrophobic moment in loops",
        io::Serialization::GetAgent( &m_HydrophobicMomentWindowSizes( 2)),
        "5"
      );
      return serial;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &EPRAccessibilityData::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Handler, ISTREAM);
      io::Serialize::Read( m_Restraints, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &EPRAccessibilityData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Handler, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Restraints, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_epr_decay.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> EPRDecay::s_Instance( GetObjectInstances().AddInstance( new EPRDecay()));

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    EPRDecay::EPRDecay() :
      m_FirstResidueIDs(),
      m_SecondResidueIDs(),
      m_Measurements(),
      m_MeasurementsFilePath()
    {
    }

    //! @brief construct from file name
    //! @param MEASUREMENT_FILE_PATH path to the file containing the results of the EPR decay measurements
    EPRDecay::EPRDecay( const std::string &MEASUREMENT_FILE_PATH) :
      m_FirstResidueIDs(),
      m_SecondResidueIDs(),
      m_Measurements(),
      m_MeasurementsFilePath( MEASUREMENT_FILE_PATH)
    {
      std::ostringstream oss;
      ReadInitializerSuccessHook( util::ObjectDataLabel(), oss);
    }

    //! @brief construct from members
    //! @param CHAIN_1 chain ID of the first residue
    //! @param SEQ_1 sequence ID of the first residue
    //! @param CHAIN_2 chain ID of the second residue
    //! @param SEQ_2 sequence ID of the second residue
    EPRDecay::EPRDecay( const char CHAIN_1, int SEQ_1, char CHAIN_2, int SEQ_2) :
      m_FirstResidueIDs( CHAIN_1, SEQ_1),
      m_SecondResidueIDs( CHAIN_2, SEQ_2),
      m_Measurements(),
      m_MeasurementsFilePath()
    {
    }

    //! @brief copy constructor
    //! @return pointer to a new EPRDecay
    EPRDecay *EPRDecay::Clone() const
    {
      return new EPRDecay( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &EPRDecay::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &EPRDecay::GetAlias() const
    {
      static const std::string s_alias( "EPRDecay");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer EPRDecay::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Stores EPR decay measurements.");
      serializer.AddInitializer
      (
        "file path",
        "path to the file containing the results of the EPR decay measurements",
        io::Serialization::GetAgent( &m_MeasurementsFilePath)
      );

      return serializer;
    }

    //! @brief returns the spin-labeling sites these measurements are describing
    //! @return pair of spin-labeling sites described by chain ID and sequence ID
    EPRDecay::SLPair EPRDecay::GetSpinLabelingSites() const
    {
      return SLPair( m_FirstResidueIDs, m_SecondResidueIDs);
    }

    //! @brief adds a data point to the measurements
    //! @param TIME time of the measurement
    //! @DECAY observed decay
    void EPRDecay::AddMeasurement( double TIME, double DECAY)
    {
      m_Measurements.PushBack( Measurement( TIME, DECAY));
    }

    //! @brief returns the EPR decay measurements for this spin-labeling pair
    //! @return the EPR decay measurements for this spin-labeling pair
    const storage::Vector< EPRDecay::Measurement> &EPRDecay::GetMeasurements() const
    {
      return m_Measurements;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    //! @return return code indicating success or failure
    bool EPRDecay::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      return true;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_epr_decay_simulation.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically
#include <math.h>

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> EPRDecaySimulation::s_Instance
    (
      GetObjectInstances().AddInstance( new EPRDecaySimulation())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    EPRDecaySimulation::EPRDecaySimulation() :
      m_SpinLabelingPairs()
    {
    }

    //! @brief construct from list of spin-labeling pairs
    //! @param SL_PAIRS list of spin-labeling pairs
    EPRDecaySimulation::EPRDecaySimulation( const storage::Vector< SLPair> &SL_PAIRS) :
      m_SpinLabelingPairs( SL_PAIRS)
    {
    }

    //! @brief copy constructor
    //! @return pointer to a new EPRDecaySimulation
    EPRDecaySimulation *EPRDecaySimulation::Clone() const
    {
      return new EPRDecaySimulation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &EPRDecaySimulation::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &EPRDecaySimulation::GetAlias() const
    {
      static const std::string s_alias( "EPRDecaySimulation");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer EPRDecaySimulation::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Simulates EPR decay patterns.");
      // serializer.AddInitializer
      // (
      //   "metropolis",
      //   "metropolis criterion to decide which mutates are accepted",
      //   io::Serialization::GetAgent( &m_Metropolis)
      // );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate decay for the provided distance and time
    //! @param DISTANCE spin-spin distance in angstroms
    //! @param TIME
    //! @return decay
    double EPRDecaySimulation::CalculateDecay( double DISTANCE, double TIME)
    {
      // calculate factors that are independent of the integration
      const double epr_constant( 326.08);
      const double inverse_distance( pow( DISTANCE / 10.0, 3));

      // integrate over the integration points
      double decay( 0.0);
      const size_t number_bins( 201);
      for( size_t bin_index( 0); bin_index < number_bins; ++bin_index)
      {
        const double fy( ( 1.0 - 3.0 * pow( 0.005 * bin_index, 2)) * epr_constant * inverse_distance);
        decay += cos( fy * TIME);
      }

      return decay;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_epr_distance_data.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_static.h"
#include "fold/bcl_fold_score_weight_set.h"
#include "score/bcl_score_restraint_atom_attraction.h"
#include "score/bcl_score_restraint_atom_distance.h"
#include "score/bcl_score_restraint_distance_spin_label.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    // initialize scores
    fold::Score EPRDistanceData::e_ScoreEPRDistanceRestraint( fold::GetScores().e_Undefined);
    fold::Score EPRDistanceData::e_ScoreEPRDistanceUpperPenalty( fold::GetScores().e_Undefined);
    fold::Score EPRDistanceData::e_ScoreEPRDistanceLowerPenalty( fold::GetScores().e_Undefined);

    const util::SiPtr< const score::RestraintAtomDistance> EPRDistanceData::s_SpinLabelScore
    (
      util::Enumerated< score::RestraintAtomDistance>::AddInstance
      (
        new score::RestraintAtomDistance
        (
          score::RestraintDistanceSpinLabel(),
          1.0,
          score::RestraintDistanceSpinLabel::GetDefaultScheme()
        )
      )
    );

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> EPRDistanceData::s_Instance
    (
      util::Enumerated< Interface>::AddInstance( new EPRDistanceData())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor taking member variable parameters
    //! @param HANDLER the handler that will be used to read and write the restraints
    EPRDistanceData::EPRDistanceData() :
      m_Handler( GetDefaultHandler()),
      m_Restraints()
    {
    }

    //! @brief Clone function
    //! @return pointer to new EPRDistanceData
    EPRDistanceData *EPRDistanceData::Clone() const
    {
      return new EPRDistanceData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &EPRDistanceData::GetAlias() const
    {
      static const std::string s_name( "DistanceEPR");
      return s_name;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &EPRDistanceData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the default handler
    //! @return default handler for this class
    const HandlerAtomDistanceAssigned &EPRDistanceData::GetDefaultHandler()
    {
      static const HandlerAtomDistanceAssigned s_handler( EPRDistanceData::GetDefaultExtension());
      return s_handler;
    }

    //! @brief get the default file extension
    //! @return the default file extension
    const std::string &EPRDistanceData::GetDefaultExtension()
    {
      static const std::string s_extension( ".epr_cst_bcl");
      return s_extension;
    }

    //! @brief gives reference to the specific type of data this restraint uses
    //! @return gives reference to the specific type of data this restraint uses
    const util::ShPtrVector< AtomDistance> &EPRDistanceData::GetAtomDistanceRestraints() const
    {
      return *m_Restraints;
    }

    //! @brief gives the scoring object that is used to score this type of restraint
    //! @return the scoring object that is used to score this type of restraint
    const util::ShPtrVector< score::ProteinModel> &EPRDistanceData::GetScores() const
    {
      // initialize scores
      static util::ShPtrVector< score::ProteinModel> scores;

      // if the vector has not been filled
      if( scores.IsEmpty())
      {
        // initialize scores then add them to vector
        EPRDistanceData epr_data;
        epr_data.InitializeScores();
        scores.PushBack( *e_ScoreEPRDistanceRestraint);
        scores.PushBack( *e_ScoreEPRDistanceLowerPenalty);
        scores.PushBack( *e_ScoreEPRDistanceUpperPenalty);
      }

      // end
      return scores;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initialize the scores and add them to Scores enumerator
    void EPRDistanceData::InitializeScores()
    {
      if( !e_ScoreEPRDistanceRestraint.IsDefined())
      {
        // read in restraints
        if( !m_Restraints.IsDefined() && m_Handler.Exists())
        {
          m_Restraints = util::CloneToShPtr( m_Handler.ReadRestraintsFromFile());
        }
        util::ShPtr< score::RestraintAtomDistance> score( s_SpinLabelScore->Clone());
        score->SetRestraints( m_Restraints);
        e_ScoreEPRDistanceRestraint = fold::GetScores().AddScore( score);

        // create Histogram "histogram" which will be used to hold the histogram of SL-CB distances
        math::Histogram histogram;

        // create IFStream "read"
        io::IFStream read;

        // open "read" and bind it to the histogram file containing SL-CB distances
        io::File::MustOpenIFStream
        (
          read,
          score::Score::AddHistogramPath( score::RestraintDistanceSpinLabel::GetDefaultHistogramFilename())
        );

        // read in from "read" into "histogram"
        read >> histogram;

        // close and clear read stream
        io::File::CloseClearFStream( read);

        // score for the transition on the left (negative x-axis) side of the kb potential
        util::ShPtr< score::RestraintAtomDistance> lower_penalty
        (
          new score::RestraintAtomDistance
          (
            score::RestraintAtomAttraction
            (
              score::RestraintAtomAttraction::GetDefaultDepthRange(),
              score::RestraintAtomAttraction::GetDefaultLeftEndWell( histogram),
              score::RestraintAtomAttraction::GetDefaultTransitionWidth(),
              true
            ),
            1.0,
            "epr_lower_penalty",
            m_Restraints
          )
        );
        util::Enumerated< score::RestraintAtomDistance>::AddInstance( lower_penalty->Clone());
        e_ScoreEPRDistanceLowerPenalty = fold::GetScores().AddScore( lower_penalty);

        // score for the transition on the right (right x-axis) side of the kb potential
        util::ShPtr< score::RestraintAtomDistance> upper_penalty
        (
          new score::RestraintAtomDistance
          (
            score::RestraintAtomAttraction
            (
              score::RestraintAtomAttraction::GetDefaultDepthRange(),
              score::RestraintAtomAttraction::GetDefaultRightEndWell( histogram),
              score::RestraintAtomAttraction::GetDefaultTransitionWidth(),
              false
            ),
            1.0,
            "epr_upper_penalty",
            m_Restraints
          )
        );
        util::Enumerated< score::RestraintAtomDistance>::AddInstance( upper_penalty->Clone());
        e_ScoreEPRDistanceUpperPenalty = fold::GetScores().AddScore( upper_penalty);
      }
    }

    //! @brief sets the weights of scores in a weight set
    //! @param SCORE_WEIGHT_SET the score weight set that will be modified
    void EPRDistanceData::ModifyScoreWeightSet( fold::ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
      SCORE_WEIGHT_SET.SetWeight( e_ScoreEPRDistanceRestraint, 5);
      SCORE_WEIGHT_SET.SetWeight( e_ScoreEPRDistanceUpperPenalty, 5);
      SCORE_WEIGHT_SET.SetWeight( e_ScoreEPRDistanceLowerPenalty, 5);
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void EPRDistanceData::MergeAndModifyMutateTree( fold::MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< fold::MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }
 
    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer EPRDistanceData::GetSerializer() const
    {
      io::Serializer serial( m_Handler.GetSerializer());
      serial.SetClassDescription( "Allows use of electron paramagnetic resonance-based accessibility data measurements");
      return serial;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &EPRDistanceData::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Handler, ISTREAM);
      io::Serialize::Read( m_Restraints, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &EPRDistanceData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Handler, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Restraints, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_handler_accessibility_aa.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_default_flags.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> HandlerAccessibilityAA::s_Instance
    (
      GetObjectInstances().AddInstance( new HandlerAccessibilityAA())
    );

    // initialize "s_FileHeader"
    const std::string HandlerAccessibilityAA::s_FileHeader( "Accessibility AA");

    // initialize "s_LineFormat"
    const std::string HandlerAccessibilityAA::s_LineFormat
    (
      "CHAIN_ID AA_SEQ_ID ENVIRONMENT MEASUREMENT ENVIRONMENT MEASUREMENT ENVIRONMENT MEASUREMENT ..."
    );

    // initialize "s_MinimumEntriesPerLine"
    const size_t HandlerAccessibilityAA::s_MinimumEntriesPerLine( 4);

    // initialize "s_ChainIDColumn"
    const size_t HandlerAccessibilityAA::s_ChainIDColumn( 0);

    // initialize "s_AASeqIDColumn"
    const size_t HandlerAccessibilityAA::s_AASeqIDColumn( 1);

    // initialize "s_FirstEnvironmentColumn"
    const size_t HandlerAccessibilityAA::s_FirstEnvironmentColumn( 2);

    // initialize "s_FirstMeasurementColumn"
    const size_t HandlerAccessibilityAA::s_FirstMeasurementColumn( 3);

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    HandlerAccessibilityAA::HandlerAccessibilityAA() :
      HandlerBase< AccessibilityProfile>( ".access_bcl"),
      m_ExposureCalculator()
    {
    }

    //! @brief constructor taking member variables
    //! @param SEQUENCE_EXCLUSION sequence exclusion to usewhen calculating exposure from a structure
    //! @param THRESHOLD_LOW_HIGH the min and max distance threshold to use for calculating exposure from a structure
    HandlerAccessibilityAA::HandlerAccessibilityAA
    (
      const util::ShPtr< assemble::AAExposureInterface> &EXPOSURE_CALCULATOR
    ) :
      HandlerBase< AccessibilityProfile>( ".access_bcl"),
      m_ExposureCalculator( EXPOSURE_CALCULATOR)
    {
    }

    //! @brief virtual copy constructor
    HandlerAccessibilityAA *HandlerAccessibilityAA::Clone() const
    {
      return new HandlerAccessibilityAA( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &HandlerAccessibilityAA::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &HandlerAccessibilityAA::GetAlias() const
    {
      static const std::string s_name( "AAAccessibility");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief CreateRestraints is the function which creates the Restraints from an istream
    //! @param ISTREAM is the istream from which the restraints will be created
    //! @return returns a ShPtrVector of RestraintInterfaces
    AccessibilityProfile HandlerAccessibilityAA::ReadRestraints( std::istream &ISTREAM) const
    {
      // create ShPtrVector "restraints" to hold the restraints created from "ISTREAM"
      storage::List< AccessibilityAA> restraints;

      // create string "file_header" to hold the string at the top of "ISTREAM"
      std::string file_header;

      // get the first line of "STREAM"
      std::getline( ISTREAM, file_header);

      // assert that "file_header" contains the expected string
      BCL_Assert
      (
        file_header == s_FileHeader,
        "first line of file should be \"" + s_FileHeader + "\" but instead is \"" + file_header + "\""
      );

      // create storage::Vector of strings "lines" and initialize it with all the lines in "ISTREAM"
      storage::Vector< std::string> lines( util::StringLineListFromIStream( ISTREAM));

      // iterate through "lines" to create AccessibilityAA restraints
      for
      (
        storage::Vector< std::string>::const_iterator line_itr( lines.Begin()), line_itr_end( lines.End());
        line_itr != line_itr_end;
        ++line_itr
      )
      {
        // create Vector of strings "split_line" and initialize with the individual strings of the "current_line"
        const storage::Vector< std::string> split_line( util::SplitString( util::TrimString( *line_itr)));

        // create size_t "split_line_size" and initialize with the size of "split_line"
        const size_t split_line_size( split_line.GetSize());

        // make sure "split_line" has at least the minimum size necessary
        BCL_Assert
        (
          split_line_size >= s_MinimumEntriesPerLine,
          "File line should have at least " + util::Format()( s_MinimumEntriesPerLine) +
          " entries, but\n" + util::Format()( split_line) + "\ninstead has " +
          util::Format()( split_line_size) + " entries"
        );

        // check that "split_line" format possibly makes sense
        BCL_Assert
        (
              util::LengthOfIntegerType( split_line( s_AASeqIDColumn))          //< should be aa seq id
          && !util::LengthOfIntegerType( split_line( s_FirstEnvironmentColumn)) //< should be environment type
          &&  util::IsNumerical( split_line( s_FirstMeasurementColumn)),      //< should be the measurement
          "Line format should be +\"" + s_LineFormat + "\" but instead the line started with \""
          + util::Format()( split_line( 0)) + " " + util::Format()( split_line( 1)) + " "
          + util::Format()( split_line( 2)) + " " + util::Format()( split_line( 3)) + "\""
        );

        // create char "chain" and initialize with the chain the AccessibilityAA corresponds to
        BCL_Assert( split_line( s_ChainIDColumn).size() == 1, "chain id has not the right number of chars");
        const char chain( split_line( s_ChainIDColumn)[ 0]);

        // create int "aa_seq_id" initialize with amino acid (identified by the sequence id) of the AccessibilityAA
        const int aa_seq_id( util::ConvertStringToNumericalValue< int>( split_line( s_AASeqIDColumn)));

        // create ShPtr to find::LocatorInterface "aa_locator" and initialize with assemble::LocatorAA
        const util::ShPtr< assemble::LocatorAA> aa_locator
        (
          new assemble::LocatorAA( chain, aa_seq_id, fold::DefaultFlags::GetFlagPDBIDNumbering()->GetFlag())
        );

        storage::Map< AccessibilityAA::EnvironmentEnum, double> data;

        // iterate through "split_line" to fill "accessibility_data"
        for( size_t column( s_FirstEnvironmentColumn); column < split_line_size; ++column)
        {
          // create std::string "environment_string" and initialize with the string in "column" of "split_line"
          const std::string &environment_string( split_line( column));

          // create AccessibilityAA::EnvironmentType "current_environment_type" to hold the environment the current
          // accessibility was measured in and initialize with enum found with "environment_string"
          const AccessibilityAA::EnvironmentEnum current_environment_type( environment_string);

          // increment "column" to the measurement which immediately follows the corresponding environment
          ++column;

          // create double "accessibility_measurement" to hold the actual accessibility measurement
          const double accessibility_measurement
          (
            util::ConvertStringToNumericalValue< double>( split_line( column))
          );

          // insert "current_environment_type" and "current_accessibility" into "accessibility_data"
          BCL_Assert( data.Insert( std::make_pair( current_environment_type, accessibility_measurement)).second, "insertion failed");
          // check to make sure that the insertion was successful.
        }

        // add to "restraints" an AccessibilityAA constructed from "accessibility_data" and "aa_locator"
        restraints.PushBack( AccessibilityAA( data, aa_locator, m_ExposureCalculator));
      }

      const AccessibilityProfile profile( restraints);

      // return "restraints"
      return profile;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read restraint from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &HandlerAccessibilityAA::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ExposureCalculator, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write restraint to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &HandlerAccessibilityAA::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ExposureCalculator, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    std::ostream &HandlerAccessibilityAA::WriteRestraints( std::ostream &OSTREAM, const AccessibilityProfile &PROFILE)
    {
      // write the header
      OSTREAM << s_FileHeader << '\n';

      // iterate through the accessibilities
      for
      (
        storage::List< AccessibilityAA>::const_iterator
          accessibility_itr( PROFILE.GetAccessibilities().Begin()),
          accessibility_itr_end( PROFILE.GetAccessibilities().End());
        accessibility_itr != accessibility_itr_end;
        ++accessibility_itr
      )
      {
        const util::ShPtr< assemble::LocatorAA> &aa_locator( accessibility_itr->GetAA());
        BCL_Assert
        (
          aa_locator.IsDefined(),
          "could not cast " + util::Format()( accessibility_itr->GetAA()) + "\nto util::ShPtr< assemble::LocatorAA>"
        );

        const char chain( aa_locator->GetLocatorChain().GetChainID());
        const int aa_id( aa_locator->GetAAID());

        // write chain and aa id
        io::Serialize::Write( chain, OSTREAM, 0);
        OSTREAM << " ";
        io::Serialize::Write( aa_id, OSTREAM, 0);

        // write out the environment information
        for
        (
          storage::Map< AccessibilityAA::EnvironmentEnum, double>::const_iterator
            data_itr( accessibility_itr->GetAccessibilityAAs().Begin()),
            data_itr_end( accessibility_itr->GetAccessibilityAAs().End());
          data_itr != data_itr_end;
          ++data_itr
        )
        {
          OSTREAM << " ";
          const double accessibility_value( data_itr->second);
          OSTREAM << data_itr->first.GetString() << " ";
          io::Serialize::Write( accessibility_value, OSTREAM, 0);
        }
        OSTREAM << '\n';
      }

      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_handler_atom_distance_assigned.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_atom.h"
#include "fold/bcl_fold_default_flags.h"
#include "math/bcl_math_running_average_sd.h"
#include "restraint/bcl_restraint_atom_distance.h"
#include "restraint/bcl_restraint_data_set_pairwise.h"
#include "restraint/bcl_restraint_locator_coordinates_first_side_chain_atom.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> HandlerAtomDistanceAssigned::s_Instance
    (
      util::Enumerated< HandlerBase< util::ShPtrVector< AtomDistance> > >::AddInstance
      (
        new HandlerAtomDistanceAssigned()
      )
    );

    //! @brief gives the identifying string at the top of the restraint file
    //! @return string that identifies the file as an atom distance restraint file
    const std::string HandlerAtomDistanceAssigned::GetFileHeader()
    {
      static const std::string s_file_header( "Atom Distance Assigned");
      return s_file_header;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new HandlerAtomDistanceAssigned
    HandlerAtomDistanceAssigned *HandlerAtomDistanceAssigned::Clone() const
    {
      return new HandlerAtomDistanceAssigned( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a short name for this class
    //! @return a short name for this class
    HandlerAtomDistanceAssigned::HandlerAtomDistanceAssigned
    (
      const std::string &DEFAULT_EXTENSION,
      const std::string &DEFAULT_FORMAT,
      const double &LOWER_BOUND,
      const double &UPPER_BOUND,
      const double &DISTANCE
    ) :
      HandlerBase< util::ShPtrVector< AtomDistance> >( DEFAULT_EXTENSION),
      m_DefaultFormat( DEFAULT_FORMAT),
      m_DefaultLowerBound( LOWER_BOUND),
      m_DefaultUpperBound( UPPER_BOUND),
      m_DefaultDistance( DISTANCE)
    {
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &HandlerAtomDistanceAssigned::GetAlias() const
    {
      static const std::string s_name( "AtomDistance");
      return s_name;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &HandlerAtomDistanceAssigned::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief WriteRestraints writes restraint information to an ostream
    //! TODO use list of atom distance object
    //! @param OSTREAM the stream to which the restraint information will be written
    //! @param RESTRAINT_LIST the list of restraint information that will be written to OSTREAM
    //!        the two triplets have the two chains, seqids, and atoms needed to specify the objects of the restraint
    //!        the three double in the vectornd<3> has the distance, upper bound, and lower bound, respectively
    //! @return ostream
    std::ostream &HandlerAtomDistanceAssigned::WriteRestraints
    (
      std::ostream &OSTREAM,
      const storage::Vector
      <
        storage::Pair
        <
          storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
        >
      > &RESTRAINT_LIST
    ) const
    {
      // write the file header to the ostream
      OSTREAM << GetFileHeader() << '\n';

      // iterate through the restraint information and write it out to "OSTREAM"
      for
      (
        storage::Vector
        <
          storage::Pair
          <
            storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
          >
        >::const_iterator info_itr( RESTRAINT_LIST.Begin()), info_itr_end( RESTRAINT_LIST.End());
        info_itr != info_itr_end;
        ++info_itr
      )
      {
        // write chain id of first object in restraint
        io::Serialize::Write( info_itr->First().First().First(), OSTREAM) << ' ';

        // write the seq id of the first object in the restraint
        OSTREAM << info_itr->First().First().Second() << ' ';

        // write the atom type of the first object in the restraint
        OSTREAM << info_itr->First().First().Third().GetType().GetName() << ' ';

        // write chain id of second object in restraint
        io::Serialize::Write( info_itr->First().Second().First(), OSTREAM) << ' ';

        // write the seq id of the second object in the restraint
        OSTREAM << info_itr->First().Second().Second() << ' ';

        // write the atom type of the second object in the restraint
        OSTREAM << info_itr->First().Second().Third().GetType().GetName() << ' ';

        io::Serialize::Write( info_itr->Second().First(), OSTREAM) << ' ';
        io::Serialize::Write( info_itr->Second().Second(), OSTREAM) << ' ';
        io::Serialize::Write( info_itr->Second().Third(), OSTREAM) << '\n';
      }

      return OSTREAM;
    }

    //! @brief WriteRestraints writes restraint information to an ostream
    //! TODO use list of atom distance object
    //! @param OSTREAM the stream to which the restraint information will be written
    //! @param RESTRAINT_LIST the list of restraints information that will be written to OSTREAM
    //! @return ostream
    std::ostream &HandlerAtomDistanceAssigned::WriteRestraints
    (
      std::ostream &OSTREAM,
      const util::ShPtrVector< AtomDistance> &RESTRAINT_LIST,
      const bool &INCLUDE_ATOM_TYPE,
      const bool &INCLUDE_AA_TYPE
    )
    {
      // write the file header to the ostream
      OSTREAM << GetFileHeader() << '\n';

      // iterate through the restraints and write them to the stream
      for
      (
        util::ShPtrVector< AtomDistance>::const_iterator
          restraint_itr( RESTRAINT_LIST.Begin()), restraint_itr_end( RESTRAINT_LIST.End());
        restraint_itr != restraint_itr_end;
        ++restraint_itr
      )
      {
        // get the two data points associated with this restraint
        const AtomDistance &atom_distance( **restraint_itr);
        const assemble::LocatorAtomCoordinatesInterface &data_a( *atom_distance.GetData().First());
        const assemble::LocatorAtomCoordinatesInterface &data_b( *atom_distance.GetData().Second());

        // write chain id of first object in restraint
        io::Serialize::Write( data_a.GetChainID(), OSTREAM) << ' ';

        // write the seq id of the first object in the restraint
        OSTREAM << data_a.GetSeqID() << ' ';

        if( INCLUDE_ATOM_TYPE)
        {
          OSTREAM << data_a.GetAtomType().GetName() << ' ';
        }
        if( INCLUDE_AA_TYPE)
        {
          OSTREAM << ( data_a.GetAAType().IsDefined() ? data_a.GetAAType()->GetOneLetterCode() : 'X') << ' ';
        }
        // write chain id of second object in restraint
        io::Serialize::Write( data_b.GetChainID(), OSTREAM) << ' ';

        // write the seq id of the second object in the restraint
        OSTREAM << data_b.GetSeqID() << ' ';

        // write the atom type of the second object in the restraint
        if( INCLUDE_ATOM_TYPE)
        {
          OSTREAM << data_b.GetAtomType().GetName() << ' ';
        }
        if( INCLUDE_AA_TYPE)
        {
          OSTREAM << ( data_b.GetAAType().IsDefined() ? data_b.GetAAType()->GetOneLetterCode() : 'X') << ' ';
        }

        // get the distance and upper and lower bounds
        const double distance( atom_distance.GetDistance()->GetDistance());
        io::Serialize::Write( distance, OSTREAM) << ' ';
        io::Serialize::Write( std::max( atom_distance.GetUpperBound(), distance), OSTREAM) << ' ';
        io::Serialize::Write( std::min( atom_distance.GetLowerBound(), distance), OSTREAM) << '\n';
      }

      // return the stream
      return OSTREAM;
    }

    //! @brief creates atom distance restraints given a set of data with distances calculated from a model ensemble
    //! @param ENSEMBLE the protein model ensemble from which distances for the data pairs will be calculated
    //! @param DATA_PAIRS the list of restraints whose distances will be calculated from the model
    //! @return list of atom distance restraints that were calculated from the data pairs and the model ensemble
    util::ShPtrVector< AtomDistance> HandlerAtomDistanceAssigned::CreateRestraints
    (
      const assemble::ProteinEnsemble &ENSEMBLE,
      const DataSetPairwise &DATA_PAIRS
    )
    {
      // will hold the restraints
      util::ShPtrVector< AtomDistance> restraints;

      // iterate through the data pairs to calculate distances from the model and create restraints
      for
      (
        DataSetPairwise::const_iterator data_itr( DATA_PAIRS.Begin()), data_itr_end( DATA_PAIRS.End());
        data_itr != data_itr_end; ++data_itr
      )
      {
        // calculate the distance information for the current pair
        const storage::Pair< math::RunningAverageSD< double>, math::RunningMinMax< double> > distance_info
        (
          data_itr->EuclidianDistance( ENSEMBLE)
        );

        // reference on mean sd object
        const math::RunningAverageSD< double> &mean_sd( distance_info.First());

        // reference on min max object
        const math::RunningMinMax< double> &min_max( distance_info.Second());

        // true if distance statistics for ensemble were successful i.e. data pair could be found in the ensemble
        if( mean_sd.GetWeight() && min_max.GetMax() >= min_max.GetMin())
        {
          // create distance object
          util::ShPtr< Distance> distance( new Distance( mean_sd.GetAverage(), min_max.GetMax(), min_max.GetMin()));

          // create the restraint and add it to the list of restraints
          restraints.PushBack( util::ShPtr< AtomDistance>( new AtomDistance( *data_itr, distance)));
        }
      }

      return restraints;
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief writes a list of restraint information in rosetta format
    //! @param RESTRAINTS the list of restraint information
    //! @param OSTREAM the stream to which the restraint will be written
    //! @return std::ostream which was passed as parameter
    std::ostream &HandlerAtomDistanceAssigned::WriteDistanceRestraintsRosettaFormat
    (
      std::ostream &OSTREAM,
      const util::ShPtrVector< AtomDistance> &RESTRAINT_LIST
    )
    {
      // iterate through the vector of restraint information in order to write it to "OSTREAM"
      for
      (
        util::ShPtrVector< AtomDistance>::const_iterator restraint_itr( RESTRAINT_LIST.Begin()),
          restraint_itr_end( RESTRAINT_LIST.End());
        restraint_itr != restraint_itr_end;
        ++restraint_itr
      )
      {
        OSTREAM << util::Format().W( 10).L()( "AtomPair")
        << util::Format().W( 4).R()( ( *restraint_itr)->GetData().First()->GetAtomType().GetName()) //< first atom type
        << util::Format().W( 4).R()( ( *restraint_itr)->GetData().First()->GetSeqID()) //< first residue number
        << util::Format().W( 4).R()( ( *restraint_itr)->GetData().Second()->GetAtomType().GetName()) //< second atom type
        << util::Format().W( 4).R()( ( *restraint_itr)->GetData().Second()->GetSeqID()) //< second residue number
        << "  SPLINE  " << "EPR_DISTANCE  "
        << util::Format().W( 10).R()( ( *restraint_itr)->GetDistance()->GetDistance()) //< actual distance
        << " 1.0 " //< restraint weight
        << " 0.5"  //< bin size of epr histogram
        << '\n';
      }

      return OSTREAM;
    }

    namespace
    {
      //! @brief helper function to detect whether a set of characters contains only characters that are valid AA 1-letter codes
      bool ContainsNonAATypes( const std::set< char> &SET)
      {
        static const std::string s_non_aatypes( "BJOZ"); // only four capital letters that cannot be interpreted as an aa type
        for( std::set< char>::const_iterator itr( SET.begin()), itr_end( SET.end()); itr != itr_end; ++itr)
        {
          if( !std::isupper( *itr) || s_non_aatypes.find( *itr) != std::string::npos)
          {
            return true;
          }
        }
        return false;
      }

      //! @brief function to test if all strings in a vector have the same size, and if so, return the common size
      //! @param STRINGS the strings to test for constant field size
      //! @return the field size if all strings have the same size, undefined size_t otherwise
      size_t GetCommonFieldSize( const storage::Vector< std::string> &STRINGS)
      {
        if( STRINGS.IsEmpty())
        {
          return 0;
        }
        const size_t common_size( STRINGS( 0).size());
        for
        (
          storage::Vector< std::string>::const_iterator itr( STRINGS.Begin() + 1), itr_end( STRINGS.End());
          itr != itr_end;
          ++itr
        )
        {
          if( itr->size() != common_size)
          {
            return util::GetUndefined< size_t>();
          }
        }
        return common_size;
      }

      //! @brief function to test if all strings in a vector are integers
      //! @param STRINGS the strings to test for integralness
      //! @return true if all strings in the vector are integers
      bool AreIntegral( const storage::Vector< std::string> &STRINGS)
      {
        for
        (
          storage::Vector< std::string>::const_iterator itr( STRINGS.Begin()), itr_end( STRINGS.End());
          itr != itr_end;
          ++itr
        )
        {
          if( util::LengthOfIntegerType( *itr) != itr->size())
          {
            return false;
          }
        }
        return true;
      }

      //! @brief function to test if all strings in a vector are numerical (floating point or integers)
      //! @param STRINGS the strings to test for floating point or integers
      //! @return true if all strings in the vector represent floating point or integers
      bool AreNumerical( const storage::Vector< std::string> &STRINGS)
      {
        for
        (
          storage::Vector< std::string>::const_iterator itr( STRINGS.Begin()), itr_end( STRINGS.End());
          itr != itr_end;
          ++itr
        )
        {
          if( !util::IsNumerical( *itr))
          {
            return false;
          }
        }
        return true;
      }

      //! @return 0 if very likely chain id, 1 if very likely aa one letter code, -1 if no idea
      int IsMoreLikelyAAOneLetterCodeThanChainId
      (
        const storage::Vector< std::string> &LETTERS,
        const size_t &N_RECORDS,
        const bool &WARN_ON_AMBIGUOUS
      )
      {
        // create sets for each of the fields
        std::set< char> vals;
        for( storage::Vector< std::string>::const_iterator itr( LETTERS.Begin()), itr_end( LETTERS.End()); itr != itr_end; ++itr)
        {
          vals.insert( ( *itr)[ 0]);
        }

        if( ContainsNonAATypes( vals))
        {
          // there are letters that are invalid for amino acid types; so the letters must be chain ids
          return 0;
        }
        else if( vals.size() == size_t( 1))
        {
          if( vals.count( 'X'))
          {
            if( WARN_ON_AMBIGUOUS)
            {
              // only one letter and it is X; probably means that the format writer didn't know the aa types,
              // and so considers them unknown
              BCL_MessageCrt
              (
                "Warning: Interpreting column in contact file containing only X as indicating the amino acid type. "
                "If these are actually the chain id, then the format must be given"
              );
            }
            return 1;
          }
          else
          {
            return 0;
          }
        }
        else if( vals.size() > size_t( 10))
        {
          // more than 10 unique letters, almost certainly aa types since there should very rarely be that many chains
          return 1;
        }
        else if( N_RECORDS >= size_t( 20))
        {
          // >= 20 contacts very low chance that there weren't at least 10 unique aa types
          return 0;
        }
        else if( N_RECORDS > size_t( 2) && vals.size() >= N_RECORDS)
        {
          // fewer records than unique values. most likely aa types
          return 1;
        }
        // < 10 contacts; yet 2-4 distinct values. Could be cross-linking data or sparse predicted contacts.
        // Assert out and force the user to set the format
        if( WARN_ON_AMBIGUOUS)
        {
          BCL_MessageCrt( "Could not determine whether field is a chain ID or AA type, so format must be provided");
        }
        return -1;
      }

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //! @class HandlerAtomDistanceSerializer
      //! @brief reads contacts from most one-contact-per-line based file formats; infers format if it was not given
      //! @author mendenjl
      //! @date Nov 18, 2014
      //! @detail class is rather specific for reading contact type files, which may be at atomic or AA resolution
      //!         likewise, class is not exposed to external users since they should use HandlerAtomDistanceAssigned
      //!         instead
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      class HandlerAtomDistanceSerializer
      {
      private:

        //! Format string, given by user (and used to set the indexes below)
        std::string m_Format;

        //! Format string, if it was inferred
        bool m_FormatWasInferred;

        double m_DefaultLowerBound; //!< Default value used for lower bound if it is not contact-specific
        double m_DefaultUpperBound; //!< Default value used for upper bound if it is not contact-specific
        double m_DefaultDistance;   //!< Default value used for distance if it is not contact-specific

        enum Index
        {
          e_SeqIdA,         //!< Field Index for Sequence id, first AA
          e_SeqIdB,         //!< Field Index for Sequence id, second AA
          e_ChainIdA,       //!< Field Index for Chain ID, first AA
          e_ChainIdB,       //!< Field Index for Chain ID, second AA
          e_AATypeA,        //!< Field Index for AA Type, first AA
          e_AATypeB,        //!< Field Index for AA Type, second AA
          e_AtomTypeA,      //!< Field Index for Atom Type, first AA
          e_AtomTypeB,      //!< Field Index for Atom Type, second AA
          e_Distance,       //!< Field Index for Distance of the contact
          e_LowerLimit,     //!< Field Index for Lower limit of the contact
          e_UpperLimit,     //!< Field Index for Upper limit of the contact
          e_Confidence,     //!< Field Index for Confidence
          s_NumberIndices,
          s_NumberAASpecificIndices = e_Distance //!< Number of indexes that are specific to aa type
        };

        //! Format section; inferred from the given file or the given section
        size_t m_FormatIndices[ s_NumberIndices];

        //! whether the format has each of the fields
        bool   m_HasField[ s_NumberIndices];

        //! minimum number of fields that must be given
        size_t m_MinNumberFields;

        //! @brief get all letters in the format alphabet
        static const std::string &GetFormatAlphabet()
        {
          static const std::string s_format( "SSCCAATTDLUN-");
          return s_format;
        }

        //! Whether the given field applies to only one of the two amino acids
        static bool GetIsAASpecific( const Index &IDX)
        {
          return int( IDX) < int( s_NumberAASpecificIndices);
        }

        //! @brief get all letters in the format alphabet
        static char GetFormatLetter( const Index &IDX)
        {
          return GetFormatAlphabet()[ int( IDX)];
        }

      public:

        //! @brief get all letters in the format alphabet
        static const std::string &GetAlphabetDescription()
        {
          static const std::string s_desc
          (
            "S -> sequence id\n"
            "A -> amino acid 1 or 3 letter code\n"
            "C -> chain id; defaults to A\n"
            "T -> PDB atom type; defaults to 1st sidechain atom if not given\n"
            "D -> Distance -> the most likely distance of the contact\n"
            "L -> Lower limit of the contact (Angstrom); defaults to Distance - 2\n"
            "U -> Upper limit of the contact (Angstrom); defaults to Distance + 2\n"
            "N -> Confidence in range (0,1]; defaults to 1\n"
            "- -> Arbitrary word; not parsed or used\n"
            "AA-specific types (represented by letters SACT) can appear twice in a given format string, e.g. "
            "-CATSCATSLUN is a valid format. - can appear any number of times in the format string"
          );
          return s_desc;
        }

        //! @brief constructor from default lower, upper bounds and distance
        HandlerAtomDistanceSerializer
        (
          const double &LOWER_BOUND,
          const double &UPPER_BOUND,
          const double &DISTANCE,
          const std::string &FORMAT = ""
        ) :
          m_FormatWasInferred( true),
          m_DefaultLowerBound( LOWER_BOUND),
          m_DefaultUpperBound( UPPER_BOUND),
          m_DefaultDistance( DISTANCE),
          m_MinNumberFields( 0)
        {
          if( !FORMAT.empty())
          {
            SetFormat( FORMAT);
          }
        }

      /////////////////
      // data access //
      /////////////////

        //! @brief get the format string (inferred string or the format string given by the user)
        const std::string &GetFormat() const
        {
          return m_Format;
        }

        //! @brief tell whether the format string was inferred from the inputs (true) or provided to this class (false)
        bool GetWasFormatInferred() const
        {
          return m_FormatWasInferred;
        }

        //! @brief set format
        //! @param FORMAT the format string containing field codes given in GetFormatAlphabet(), and described in
        //!
        void SetFormat( const std::string &FORMAT);

        //! @brief infer the format from a given set of tokens
        //! @param TOKENS tokens from a file
        //! @return the inferred format string
        static std::string InferFormat( const storage::Vector< storage::Vector< std::string> > &TOKENS);

      ////////////////
      // operations //
      ////////////////

        //! @brief creates contact restraints
        //! @param TOKENS tokens from a file
        util::ShPtrVector< AtomDistance> CreateContacts
        (
          const storage::Vector< storage::Vector< std::string> > &TOKENS
        );

      };

      //! @brief set format
      //! @param FORMAT the format string; should only contain letters in
      void HandlerAtomDistanceSerializer::SetFormat( const std::string &FORMAT)
      {
        for( size_t i( 0); i < size_t( s_NumberIndices); ++i)
        {
          m_FormatIndices[ i] = util::GetUndefined< size_t>();
          m_HasField[ i] = false;
        }
        m_MinNumberFields = FORMAT.size();
        for( size_t i( 0), sz( FORMAT.size()); i < sz; ++i)
        {
          size_t pos( GetFormatAlphabet().find( FORMAT[ i]));
          BCL_Assert
          (
            pos != std::string::npos,
            "Unknown format character; allowed characters must be one of\n" + GetAlphabetDescription()
          );
          if( pos >= s_NumberIndices) // letter -
          {
            continue;
          }
          if( m_HasField[ pos])
          {
            BCL_Assert
            (
              pos < int( s_NumberAASpecificIndices),
              "Duplicate field for " + util::Format()( FORMAT[ i]) + "; should only be one of these values per contact!"
            );

            ++pos;
            BCL_Assert
            (
              !m_HasField[ pos],
              "AA-specific field: " + util::Format()( FORMAT[ i])
              + " cannot be given more than twice in given format string: " + FORMAT
            );
          }
          m_HasField[ pos] = true;
          m_FormatIndices[ pos] = i;
        }
        BCL_Assert
        (
          m_HasField[ int( e_SeqIdA)] && m_HasField[ int( e_SeqIdB)],
          "At least two sequence ids must be present in the given format"
        );
        BCL_Assert
        (
          m_HasField[ int( e_AATypeA)] == m_HasField[ int( e_AATypeB)],
          "Either 0 or 2 amino acid types must be given, not 1"
        );
        BCL_Assert
        (
          m_HasField[ int( e_AtomTypeA)] == m_HasField[ int( e_AtomTypeB)],
          "Either 0 or 2 atom types are required in the format string, not 1"
        );
        // allow a single chain id; since intra-chain contacts are the most common
        if( m_HasField[ int( e_ChainIdA)] != m_HasField[ int( e_ChainIdB)])
        {
          m_HasField[ int( e_ChainIdB)] = m_HasField[ int( e_ChainIdA)];
          m_FormatIndices[ int( e_ChainIdB)] = m_FormatIndices[ int( e_ChainIdA)];
        }
        m_FormatWasInferred = false;
      }

      // possibilities for each field
      enum BasicPossibilities
      {
        e_SingleLetter   , // a letter A-Za-z
        e_ThreeLetterCode, // three letter AA code
        e_SingleNumber   , // single digit number
        e_MultipleNumber , // number with multiple digits
        e_FloatingPoint  , // floating point number
        e_AtomType       , // biol::AtomType
        e_Word           , // anything else
        e_Unknown          // initial assignment
      };

      //! @brief get the three letter codes for all amino acids as a set
      storage::Set< std::string> GetThreeLetterCodeSet()
      {
        storage::Set< std::string> three_letter_codes;
        for
        (
          biol::AATypes::const_iterator itr( biol::GetAATypes().Begin()), itr_end( biol::GetAATypes().End());
          itr != itr_end;
          ++itr
        )
        {
          three_letter_codes.Insert( ( *itr)->GetThreeLetterCode());
        }
        return three_letter_codes;
      }

      //! @brief infer the format from a given set of tokens
      //! @param TOKENS tokens from a file
      //! @return the inferred format string
      std::string HandlerAtomDistanceSerializer::InferFormat
      (
        const storage::Vector< storage::Vector< std::string> > &TOKENS
      )
      {
        BCL_Assert( TOKENS.GetSize() && TOKENS( 0).GetSize(), "Cannot infer format when no tokens were given!");

        // number of fields (TOKENS should be ordered as [field][contact])
        const size_t n_fields( TOKENS.GetSize());
        const size_t n_contacts( TOKENS( 0).GetSize());
        std::vector< BasicPossibilities> basic_field_types( n_fields, e_Unknown);

        // determine field types using all contacts for that field
        for( size_t field_n( 0); field_n < n_fields; ++field_n)
        {
          // get all tokens for this field
          const storage::Vector< std::string> &field_tokens( TOKENS( field_n));

          // initial assigment for the first contact
          const std::string &first_contact( field_tokens( 0));

          // determine whether this field has constant size
          const size_t common_size( GetCommonFieldSize( field_tokens));

          BasicPossibilities &type( basic_field_types[ field_n]);
          if( common_size == size_t( 1))
          {
            type = AreIntegral( field_tokens) ? e_SingleNumber : e_SingleLetter;
          }
          else if( AreIntegral( field_tokens))
          {
            type = e_MultipleNumber;
          }
          else if( AreNumerical( field_tokens))
          {
            type = e_FloatingPoint;
          }
          else
          {
            // appears to be a word type
            type = e_Word;

            // test for special words : aa types and atom types
            const bool could_be_aa_types
            (
              common_size == size_t( 3)
              && biol::GetAATypes().AATypeFromThreeLetterCode( first_contact).IsDefined()
            );
            const bool could_be_atom_types( biol::GetAtomTypes().HaveEnumWithName( first_contact));
            if( could_be_aa_types || could_be_atom_types)
            {
              // plausible aa type or atom type. Check that all values are like that
              // get the set of all strings in the vector
              storage::Set< std::string> unique_field_vals( field_tokens.Begin(), field_tokens.End());
              if( could_be_aa_types)
              {
                static const storage::Set< std::string> aa_type_three_letter_codes( GetThreeLetterCodeSet());
                // test whether each string is an aa type
                if( unique_field_vals.IsSubsetOf( aa_type_three_letter_codes))
                {
                  type = e_ThreeLetterCode;
                }
              }
              else if( could_be_atom_types)
              {
                // test whether each string is an atom type
                static const storage::Set< std::string> atom_types( biol::GetAtomTypes().Begin(), biol::GetAtomTypes().End());
                if( unique_field_vals.IsSubsetOf( atom_types))
                {
                  type = e_AtomType;
                }
              }
            }
          }
        }

        // record the index of each field based on its basic possibility type
        storage::Vector< size_t> single_letter_fields, three_letter_fields, single_digit_fields, multi_digit_fields;
        storage::Vector< size_t> flt_fields, atom_type_fields;
        for( size_t i( 0); i < n_fields; ++i)
        {
          switch( basic_field_types[ i])
          {
            case e_SingleLetter:    single_letter_fields.PushBack( i); break;
            case e_ThreeLetterCode: three_letter_fields.PushBack( i);  break;
            case e_SingleNumber:    single_digit_fields.PushBack( i);  break;
            case e_MultipleNumber:  multi_digit_fields.PushBack( i);   break;
            case e_FloatingPoint:   flt_fields.PushBack( i);           break;
            case e_AtomType:        atom_type_fields.PushBack( i);     break;
            default: break;
          };
        }

        enum FieldType
        {
          e_FTOneLetterAACode,
          e_FTChainId,
          e_FTThreeLetterAACode,
          e_FTSeqId,
          e_FTBiolAtomType,
          e_FTDistance,
          e_FTLowerLimit,
          e_FTUpperLimit,
          e_FTConfidence,
          s_FTNumberFieldTypes
        };

        storage::Vector< storage::Vector< size_t> > field_type_indices
        (
          static_cast< size_t>( s_FTNumberFieldTypes),
          storage::Vector< size_t>()
        );

        // handle the case where the sequence ids cannot reasonably be inferred from the field types
        if( multi_digit_fields.GetSize() < 2)
        {
          if( single_digit_fields.GetSize() + multi_digit_fields.GetSize() >= 2)
          {
            multi_digit_fields.PushBack( single_digit_fields( 0));
            single_digit_fields.RemoveElements( 0, 1);
            if( multi_digit_fields.GetSize() == 1)
            {
              multi_digit_fields.PushBack( single_digit_fields( 1));
              single_digit_fields.RemoveElements( 0, 1);
            }
          }
          else
          {
            BCL_Exit( "There were not at least two columns that could be sequence ids! Specify format", -1);
          }
        }
        // remaining single digit fields can be treated as floats since we're not considering numeric chain ids
        flt_fields.Append( single_digit_fields);
        single_digit_fields.Reset();
        while( multi_digit_fields.GetSize() > 2 && flt_fields.GetSize() < 3)
        {
          flt_fields.PushBack( multi_digit_fields.LastElement());
          multi_digit_fields.PopBack();
        }
        BCL_Assert( multi_digit_fields.GetSize() == size_t( 2), "Exactly two sequence ids are required");
        field_type_indices( size_t( e_FTSeqId)) = multi_digit_fields;

        BCL_Assert( flt_fields.GetSize() < size_t( 5), "Cannot use more than 4 floating point fields");
        BCL_Assert( atom_type_fields.GetSize() < size_t( 3), "Cannot use more than 2 atom types");
        BCL_Assert( three_letter_fields.GetSize() < size_t( 3), "Cannot use more than 2 three letter codes");

        field_type_indices( size_t( e_FTBiolAtomType)) = atom_type_fields;
        field_type_indices( size_t( e_FTThreeLetterAACode)) = three_letter_fields;

        // next, try to determine which, if any, fields match the chain id and aa type. For inference, assume that chain
        // ids cannot be numeric, since single digits could be too many other things
        if( single_letter_fields.GetSize())
        {
          BCL_Assert
          (
            single_letter_fields.GetSize() < 5
            && single_letter_fields.GetSize() != size_t( 3),
            "Must be 1, 2, or 4 single letter fields (aa types and chain ids)"
          );

          // Determine which fields match the aa types vs. the chain id.
          // if three letter codes are given and only two one letter fields are present, assume they are chain ids
          if( three_letter_fields.GetSize() && single_letter_fields.GetSize() == size_t( 2))
          {
            field_type_indices( size_t( e_FTChainId)) = single_letter_fields;
          }
          else if( single_letter_fields.GetSize() == size_t( 1))
          {
            // only one one letter code; almost certainly chain id
            field_type_indices( size_t( e_FTChainId)) = single_letter_fields;
          }
          else if( single_letter_fields.GetSize() == size_t( 2))
          {
            storage::Vector< std::string> all_one_letter_codes( TOKENS( single_letter_fields( 0)));
            all_one_letter_codes.Append( TOKENS( single_letter_fields( 1)));
            const int response( IsMoreLikelyAAOneLetterCodeThanChainId( all_one_letter_codes, n_contacts, true));
            if( response == 1 || response == -1)
            {
              field_type_indices( size_t( e_FTOneLetterAACode)) = single_letter_fields;
            }
            else if( response == 0 || response == -2)
            {
              field_type_indices( size_t( e_FTChainId)) = single_letter_fields;
            }
            else
            {
              BCL_Exit( "Type was too ambiguous; format must be specified", -1);
            }
          }
          else // four fields; two must be chain ids, two must be aa types
          {
            std::vector< int> responses;
            responses.push_back( IsMoreLikelyAAOneLetterCodeThanChainId( TOKENS( single_letter_fields( 0)), n_contacts, false));
            responses.push_back( IsMoreLikelyAAOneLetterCodeThanChainId( TOKENS( single_letter_fields( 1)), n_contacts, false));
            responses.push_back( IsMoreLikelyAAOneLetterCodeThanChainId( TOKENS( single_letter_fields( 2)), n_contacts, false));
            responses.push_back( IsMoreLikelyAAOneLetterCodeThanChainId( TOKENS( single_letter_fields( 3)), n_contacts, false));
            for( int i( 0); i < 4; ++i)
            {
              if( responses[ i] < 0)
              {
                responses[ i] += 2;
              }
            }
            std::set< int> responses_set( responses.begin(), responses.end());
            if( responses_set.size() == size_t( 1))
            {
              BCL_Exit( "Type was too ambiguous; format must be specified", -1);
            }
            else
            {
              if( responses_set.size() == size_t( 3))
              {
                if( std::count( responses.begin(), responses.end(), -1) != 1)
                {
                  BCL_Exit( "Type was too ambiguous; format must be specified", -1);
                }
                if( std::count( responses.begin(), responses.end(), 0) == 1)
                {
                  std::replace( responses.begin(), responses.end(), -1, 0);
                }
                else
                {
                  std::replace( responses.begin(), responses.end(), -1, 1);
                }
              }
              const int val_1( responses[ 0]);
              if( std::count( responses.begin(), responses.end(), responses[ 0]) != 2)
              {
                BCL_Exit( "Type was too ambiguous; format must be specified", -1);
              }
              const int val_2( ( responses[ 1] + responses[ 2] + responses[ 3] - val_1) / 2);
              int val_max( std::max( val_1, val_2));
              int val_min( std::min( val_1, val_2));
              int aa_type_id_col_id;
              if( val_max == 1)
              {
                aa_type_id_col_id = val_max;
              }
              else // if( val_max == 0)
              {
                aa_type_id_col_id = val_min;
              }
              for( int i( 0); i < 4; ++i)
              {
                field_type_indices
                (
                  size_t( responses[ i] == aa_type_id_col_id ? e_FTOneLetterAACode : e_FTChainId)
                ).PushBack( single_letter_fields( i));
              }
            }
          }
        }

        // at this point, numeric fields e_SingleNumber must be one of the following:
        // A. Lower limit for the contact restraint (e.g. 0)
        // B. Confidence for the contact restraint (e.g. 1)
        // C. Upper limit for the contact restraint (e.g. 8)
        // D. Chain ID (e.g. 2) In practice this is rare, so it'll be ignored.
        // e_FloatingPoint can be any of the first three

        storage::Vector< linal::Vector< double> > floating_point_vectors
        (
          flt_fields.GetSize(),
          linal::Vector< double>( n_contacts)
        );
        size_t vecn( 0);
        linal::Vector< double> max_flts( flt_fields.GetSize()), min_flts( flt_fields.GetSize());
        for
        (
          storage::Vector< size_t>::const_iterator itr_field( flt_fields.Begin()), itr_field_end( flt_fields.End());
          itr_field != itr_field_end;
          ++itr_field, ++vecn
        )
        {
          linal::Vector< double> &vec_ref( floating_point_vectors( vecn));
          const storage::Vector< std::string> &vec_strings( TOKENS( *itr_field));
          for( size_t contact_n( 0); contact_n < n_contacts; ++contact_n)
          {
            vec_ref( contact_n) = util::ConvertStringToNumericalValue< double>( vec_strings( contact_n));
          }
          max_flts( vecn) = vec_ref.Max();
          min_flts( vecn) = vec_ref.Min();
        }

        storage::Vector< size_t> confidence_candidates, limit_candidates;
        storage::Vector< double> limit_maxes;
        for( size_t i( 0); i < flt_fields.GetSize(); ++i)
        {
          // determine whether the field could indicate the confidence for the contact, which should be in the range (0,1]
          if
          (
            min_flts( i) > double( 0.0)
            && max_flts( i) <= double( 1.0)
            && ( min_flts( i) != max_flts( i) || max_flts( i) == double( 1.0))
          )
          {
            confidence_candidates.PushBack( flt_fields( i));
          }
          else
          {
            // must be a distance (upper/lower limit or expected distance)
            limit_candidates.PushBack( flt_fields( i));
            limit_maxes.PushBack( max_flts( i));
          }
        }
        BCL_Assert
        (
          confidence_candidates.GetSize() <= size_t( 1),
          "More than one column appears to be a confidence value; format will need to be specified"
        );
        BCL_Assert
        (
          limit_candidates.GetSize() <= size_t( 3),
          "More than three columns appears to be limits/distances; format will need to be specified"
        );
        if( confidence_candidates.GetSize() == size_t( 1))
        {
          field_type_indices( e_FTConfidence).PushBack( confidence_candidates( 0));
        }
        // if there is only one field, it will either be considered the distance or the confidence
        if( limit_candidates.GetSize() == size_t( 1))
        {
          field_type_indices( e_FTDistance).PushBack( limit_candidates( 0));
        }
        // if there are two fields, then the numbers could represent and upper and lower limit or a distance
        // and upper or lower limit. Since there is a logical lower limit (0) but not a logical upper limit, the safest
        // assumption seems to be that one of them represents a distance and the other the upper limit
        else if( limit_candidates.GetSize() == size_t( 2))
        {
          if( limit_maxes( 0) >= limit_maxes( 1))
          {
            field_type_indices( e_FTUpperLimit).PushBack( limit_candidates( 0));
            field_type_indices( e_FTDistance).PushBack( limit_candidates( 1));
          }
          else
          {
            field_type_indices( e_FTUpperLimit).PushBack( limit_candidates( 1));
            field_type_indices( e_FTDistance).PushBack( limit_candidates( 0));
          }
        }
        else if( limit_candidates.GetSize() == size_t( 3))
        {
          std::vector< std::pair< double, size_t> > val_to_index( 3);
          val_to_index[ 0] = std::make_pair( limit_maxes( 0), limit_candidates( 0));
          val_to_index[ 1] = std::make_pair( limit_maxes( 1), limit_candidates( 1));
          val_to_index[ 2] = std::make_pair( limit_maxes( 2), limit_candidates( 2));
          std::sort( val_to_index.begin(), val_to_index.end());
          field_type_indices( e_FTLowerLimit).PushBack( val_to_index[ 0].second);
          field_type_indices( e_FTDistance).PushBack( val_to_index[ 1].second);
          field_type_indices( e_FTUpperLimit).PushBack( val_to_index[ 2].second);
        }
        std::string format_str( n_fields, '-');
        const std::string field_alphabet( "ACASTDLUN");
        for( size_t field_id( 0); field_id < size_t( s_FTNumberFieldTypes); ++field_id)
        {
          // get the character for this field type
          const char c( field_alphabet[ field_id]);
          for
          (
            storage::Vector< size_t>::const_iterator
              itr_field_id( field_type_indices( field_id).Begin()),
              itr_field_id_end( field_type_indices( field_id).End());
            itr_field_id != itr_field_id_end;
            ++itr_field_id
          )
          {
            format_str[ *itr_field_id] = c;
          }
        }
        BCL_MessageVrb( "Inferred format: " + format_str);
        return format_str;
      }

      //! @brief creates contact restraints
      //! @param TOKENS tokens from a file
      util::ShPtrVector< AtomDistance> HandlerAtomDistanceSerializer::CreateContacts
      (
        const storage::Vector< storage::Vector< std::string> > &TOKENS
      )
      {
        if( !m_MinNumberFields)
        {
          SetFormat( InferFormat( TOKENS));
        }
        BCL_Assert
        (
          TOKENS.GetSize() >= m_MinNumberFields,
          "At least " + util::Format()( m_MinNumberFields) + " of fields should have been present in the file"
        );
        const size_t n_contacts( TOKENS( 0).GetSize());
        util::ShPtrVector< AtomDistance> restraints;
        restraints.AllocateMemory( n_contacts);
        storage::Vector< std::string> empty_vec;
        const storage::Vector< std::string> &chain_ids_a
        (
          m_HasField[ e_ChainIdA]
          ? TOKENS( m_FormatIndices[ size_t( e_ChainIdA)])
          : empty_vec
        );
        const storage::Vector< std::string> &chain_ids_b
        (
          m_HasField[ e_ChainIdB]
          ? TOKENS( m_FormatIndices[ size_t( e_ChainIdB)])
          : empty_vec
        );
        const storage::Vector< std::string> &seq_ids_a( TOKENS( m_FormatIndices[ e_SeqIdA]));
        const storage::Vector< std::string> &seq_ids_b( TOKENS( m_FormatIndices[ e_SeqIdB]));
        const storage::Vector< std::string> &aa_types_a
        (
          m_HasField[ e_AATypeA]
          ? TOKENS( m_FormatIndices[ e_AATypeA])
          : empty_vec
        );
        const storage::Vector< std::string> &aa_types_b
        (
          m_HasField[ e_AATypeB]
          ? TOKENS( m_FormatIndices[ e_AATypeB])
          : empty_vec
        );
        const storage::Vector< std::string> &atom_types_a
        (
          m_HasField[ e_AtomTypeA]
          ? TOKENS( m_FormatIndices[ e_AtomTypeA])
          : empty_vec
        );
        const storage::Vector< std::string> &atom_types_b
        (
          m_HasField[ e_AtomTypeB]
          ? TOKENS( m_FormatIndices[ e_AtomTypeB])
          : empty_vec
        );
        const storage::Vector< std::string> &distances
        (
          m_HasField[ e_Distance]
          ? TOKENS( m_FormatIndices[ e_Distance])
          : empty_vec
        );
        const storage::Vector< std::string> &lower_limits
        (
          m_HasField[ e_LowerLimit]
          ? TOKENS( m_FormatIndices[ e_LowerLimit])
          : empty_vec
        );
        const storage::Vector< std::string> &upper_limits
        (
          m_HasField[ e_UpperLimit]
          ? TOKENS( m_FormatIndices[ e_UpperLimit])
          : empty_vec
        );
        const storage::Vector< std::string> &confidences
        (
          m_HasField[ e_Confidence]
          ? TOKENS( m_FormatIndices[ e_Confidence])
          : empty_vec
        );

        if( !m_FormatWasInferred)
        {
          // format given by user, validate string sizes where appropriate
          if( !aa_types_a.IsEmpty())
          {
            const size_t field_size_a( GetCommonFieldSize( aa_types_a));
            const size_t field_size_b( GetCommonFieldSize( aa_types_b));
            BCL_Assert
            (
              ( field_size_a == size_t( 1) || field_size_a == 3) && field_size_a == field_size_b,
              "AA types should have the same size and be either 1 or 3 letters in length"
            );
          }
          if( !chain_ids_a.IsEmpty())
          {
            BCL_Assert
            (
              GetCommonFieldSize( chain_ids_a) == size_t( 1) && GetCommonFieldSize( chain_ids_b) == size_t( 1),
              "Chain ids should both be 1 letter in length. Contact format specified incorrectly"
            );
          }
        }

        const bool is_pdb_id( fold::DefaultFlags::GetFlagPDBIDNumbering()->GetFlag());
        for( size_t contact_n( 0); contact_n < n_contacts; ++contact_n)
        {
          const char chain_a( chain_ids_a.IsEmpty() ? 'A' : chain_ids_a( contact_n)[ 0]);
          const char chain_b( chain_ids_b.IsEmpty() ? 'A' : chain_ids_b( contact_n)[ 0]);
          const int seq_id_a( util::ConvertStringToNumericalValue< int>( seq_ids_a( contact_n)));
          const int seq_id_b( util::ConvertStringToNumericalValue< int>( seq_ids_b( contact_n)));
          biol::AtomType atom_type_a, atom_type_b;
          if( !atom_types_a.IsEmpty())
          {
            atom_type_a = biol::AtomType( atom_types_a( contact_n));
            atom_type_b = biol::AtomType( atom_types_b( contact_n));
          }

          // determine AA types
          biol::AAType aa_type_a, aa_type_b;
          if( !aa_types_a.IsEmpty())
          {
            aa_type_a = aa_types_a( contact_n).size() == size_t( 1)
                        ? biol::GetAATypes().AATypeFromOneLetterCode( aa_types_a( contact_n)[ 0])
                        : biol::GetAATypes().AATypeFromThreeLetterCode( aa_types_a( contact_n));
            aa_type_b = aa_types_b( contact_n).size() == size_t( 1)
                        ? biol::GetAATypes().AATypeFromOneLetterCode( aa_types_b( contact_n)[ 0])
                        : biol::GetAATypes().AATypeFromThreeLetterCode( aa_types_b( contact_n));
          }

          // serialize distance-related strings
          const double distance
          (
            distances.IsEmpty()
            ? m_DefaultDistance
            : util::ConvertStringToNumericalValue< double>( distances( contact_n))
          );
          const double lower_bound
          (
            lower_limits.IsEmpty()
            ? m_DefaultLowerBound
            : util::ConvertStringToNumericalValue< double>( lower_limits( contact_n))
          );
          const double upper_bound
          (
            upper_limits.IsEmpty()
            ? m_DefaultUpperBound
            : util::ConvertStringToNumericalValue< double>( upper_limits( contact_n))
          );
          util::ShPtr< Distance> distance_obj( new Distance( distance, upper_bound, lower_bound));

          // create atom or just aa locators, depending on what was given
          util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_a, locator_b;
          if( !atom_types_a.IsEmpty())
          {
            locator_a =
              util::ShPtr< assemble::LocatorAtomCoordinatesInterface>
              (
                new assemble::LocatorAtom( chain_a, seq_id_a, atom_type_a, aa_type_a, is_pdb_id)
              );
            locator_b =
              util::ShPtr< assemble::LocatorAtomCoordinatesInterface>
              (
                new assemble::LocatorAtom( chain_b, seq_id_b, atom_type_b, aa_type_b, is_pdb_id)
              );
          }
          else
          {
            locator_a =
              util::ShPtr< assemble::LocatorAtomCoordinatesInterface>
              (
                new LocatorCoordinatesFirstSideChainAtom( chain_a, seq_id_a, aa_type_a, is_pdb_id)
              );
            locator_b =
              util::ShPtr< assemble::LocatorAtomCoordinatesInterface>
              (
                new LocatorCoordinatesFirstSideChainAtom( chain_b, seq_id_b, aa_type_b, is_pdb_id)
              );
          }
          const double confidence
          (
            confidences.IsEmpty()
            ? 1.0
            : util::ConvertStringToNumericalValue< double>( confidences( contact_n))
          );

          // PushBack a ShPtr to a new DistanceAssigned restraint into "restraint"
          restraints.PushBack( util::ShPtr< AtomDistance>( new AtomDistance( locator_a, locator_b, distance_obj, confidence)));
        }
        return restraints;
      }
    }

    //! @brief reads atom distance restraints from an input stream
    //! @brief input stream to read the restraints from
    //! @return the read in restraints
    util::ShPtrVector< AtomDistance> HandlerAtomDistanceAssigned::ReadRestraints( std::istream &ISTREAM) const
    {
      // create all lines
      storage::Vector< std::string> lines( util::StringLineListFromIStream( ISTREAM));

      // walk through each line.  If a number is followed by a letter, separate them with a space.
      // This handles formats that directly append the residue type to the residue number
      // convert all punctuation other than _ and . to a space. Erase lines that start with #,/, or REMARK

      storage::Vector< std::string> valid_lines;
      valid_lines.AllocateMemory( lines.GetSize());

      // record valid fields
      storage::Vector< storage::Vector< std::string> > valid_fields;

      size_t last_line_nfields( 0);
      for
      (
        storage::Vector< std::string>::const_iterator itr( lines.Begin()), itr_end( lines.End());
        itr != itr_end;
        ++itr
      )
      {
        if( itr->empty())
        {
          continue;
        }
        const std::string s( *itr);

        // ignore comments
        if( s[ 0] == '!' || s[ 0] == '#')
        {
          continue;
        }
        if( util::StartsWith( s, "REMARK"))
        {
          // skip remark lines
          continue;
        }

        // If a number is followed by a letter, separate them with a space.
        // This handles formats that directly append the residue type to the residue number

        std::string s_new;
        s_new.reserve( s.size());
        // in fp indicates whether a . was seen prior to the last separator
        bool last_was_digit( false), in_fp( false);
        size_t current_digit_field_size( 0);
        size_t pos( 0);
        size_t n_numbers( 0);
        for( std::string::const_iterator itr_s( s.begin()), itr_send( s.end()); itr_s != itr_send; ++itr_s, ++pos)
        {
          bool isa( std::isalpha( *itr_s));
          bool isd( !isa && std::isdigit( *itr_s));
          if( isa && last_was_digit)
          {
            s_new += ' ';
          }

          if( !isd)
          {
            current_digit_field_size = 0;
          }
          else if( !in_fp)
          {
            if( !current_digit_field_size)
            {
              ++n_numbers;
            }
            current_digit_field_size += 1;
            if( current_digit_field_size == size_t( 6))
            {
              // issue a message, since this is a somewhat shady operation
              BCL_MessageCrt
              (
                "Separating line " + s + " at index " + util::Format()( pos) + "; assuming that these are sequence ids "
                "of residues with seq id > 9999"
              );
              // maximum PDB ID length is 5; so if something that is 5 digits long is found, automatically insert a space
              s_new += ' ';
              current_digit_field_size = 0;
            }
          }

          // convert punctuation to space
          if( std::ispunct( *itr_s) && *itr_s != '_' && *itr_s != '.')
          {
            s_new += ' ';
          }
          else
          {
            s_new += *itr_s;
          }
          last_was_digit = isd;
          if( *itr_s == '.')
          {
            in_fp = true;
          }
          else if( !isd)
          {
            in_fp = false;
          }
        }

        if( n_numbers < 2)
        {
          // skip lines without at least two integers on them
          continue;
        }

        // split the line
        storage::Vector< std::string> split_line( util::SplitString( s_new, " \t"));

        if( last_line_nfields && split_line.GetSize() != last_line_nfields)
        {
          BCL_Exit
          (
            "Could not read format for restraint file; lines had different #s of fields even after normalization. "
            "Last line read was: " + s,
            -1
          );
        }
        last_line_nfields = split_line.GetSize();
        if( valid_fields.IsEmpty())
        {
          valid_fields.Resize( last_line_nfields);
        }
        for( size_t i( 0); i < last_line_nfields; ++i)
        {
          valid_fields( i).PushBack( split_line( i));
        }
      }

      HandlerAtomDistanceSerializer contacts_serializer
      (
        m_DefaultLowerBound,
        m_DefaultUpperBound,
        m_DefaultDistance,
        m_Format
      );
      return contacts_serializer.CreateContacts( valid_fields);
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer HandlerAtomDistanceAssigned::GetSerializer() const
    {
      io::Serializer serializer( HandlerInterface::GetSerializer());
      serializer.SetClassDescription
      (
        "Reads most contact formats that have one contact per line. Format type can be inferred from the given files in "
        "the vast majority of cases. Format can be provided or adjusted by user to handle unusual cases, such as sparse "
        "contacts between multiple chains. Lines that lack at least two integers, as well as lines that begin with # or !, "
        "are automatically ignored. For formats that lack an explicitly given upper and lower bounds, or distance, a default "
        "may be provided"
      );
      serializer.AddInitializer
      (
        "format",
        "format for the file; allows overriding the inferred format type if it is incorrect or for files for which the "
        "types of each field cannot be inferred. Valid format strings use letters to indicate the various fields: "
        + HandlerAtomDistanceSerializer::GetAlphabetDescription() + "\n"
        "It is not necessary or allowed to indicate the delimiters between each field; they will be determined automatically. "
        + std::string( m_DefaultFormat.empty() ? "" : " The most common format for files of this type is: " + m_DefaultFormat),
        io::Serialization::GetAgent( &m_Format),
        ""
      );
      serializer.AddInitializer
      (
        "lower bound",
        "Lower bound for all contact restraints, if not given in the file",
        io::Serialization::GetAgent( &m_DefaultLowerBound),
        "0"
      );
      serializer.AddInitializer
      (
        "upper bound",
        "Upper bound for all contact restraints, if not given in the file",
        io::Serialization::GetAgent( &m_DefaultUpperBound),
        "10"
      );
      serializer.AddInitializer
      (
        "distance",
        "Distance for all contact restraints, if not given in the file",
        io::Serialization::GetAgent( &m_DefaultDistance),
        "8"
      );
      return serializer;
    }

  } // namespace restraint

} // namespace bcl
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
#include "restraint/bcl_restraint_handler_body.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "restraint/bcl_restraint_body.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> HandlerBody::s_Instance
    (
      GetObjectInstances().AddInstance( new HandlerBody())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    HandlerBody::HandlerBody() :
      m_DetermineOccupancy()
    {
    }

    //! @brief constructor from a ShPtr to a FunctionInterface defining how occupancy is determined
    //! @param OCCUPANCY ShPtr to a FunctionInterface defining how occupancy is determined
    HandlerBody::HandlerBody
    (
      const util::ShPtr< util::BinaryFunctionInterface< assemble::SSEGeometryInterface, assemble::SSE, bool> > &OCCUPANCY
    ) :
      m_DetermineOccupancy( OCCUPANCY)
    {
      BCL_Assert
      (
        m_DetermineOccupancy.IsDefined(),
        "HandlerBody::HandlerBody occupancy constructor : m_DetermineOccupancy is not defined"
      );
    }

    //! @brief virtual copy constructor
    HandlerBody *HandlerBody::Clone() const
    {
      return new HandlerBody( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &HandlerBody::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief CreateRestraintsBody is the function which creates the Restraints from an istream
    //! @param ISTREAM is the istream from which the restraints will be created
    //! @return returns a ShPtrVector of RestraintInterfaces
    util::ShPtrVector< Body>
    HandlerBody::CreateRestraintsBody( std::istream &ISTREAM) const
    {
      // create the density rods
      const util::ShPtr< util::ShPtrVector< assemble::SSEGeometryInterface> > restraint_bodies( GenerateDensityRods( ISTREAM));

      // put density rods into body restraint and return body restraint
      return util::ShPtrVector< Body>( 1, util::ShPtr< Body>( new Body( restraint_bodies, m_DetermineOccupancy)));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read restraint from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &HandlerBody::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_DetermineOccupancy, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write restraint to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &HandlerBody::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_DetermineOccupancy, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief CreateRestraintsBody is the function which creates the Restraints from an istream
    //! @param ISTREAM is the istream from which the restraints will be created
    //! @return returns a ShPtrVector of RestraintInterfaces
    util::ShPtr< util::ShPtrVector< assemble::SSEGeometryInterface> > HandlerBody::GenerateDensityRods( std::istream &ISTREAM) const
    {
      // read pdb file
      pdb::Handler pdb( ISTREAM);
      ISTREAM.clear();

      // create "factory" to create protein model with amino acids of type AABackBone
      pdb::Factory factory( biol::GetAAClasses().e_AABackBone);

      // create ProteinModel "protein_model" from "pdb"
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 9;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel protein_model( factory.ProteinModelFromPDB( pdb, ssetype_min_size));

      // get the secondary structure elements of "protein_model" which will be the body restraints
      util::SiPtrVector< const assemble::SSE> sses( protein_model.GetSSEs());

      BCL_MessageStd( "number of SSEs to create bodies from: " + util::Format()( sses.GetSize()));
      BCL_MessageStd( "number of chains: " + util::Format()( protein_model.GetNumberOfChains()));
      BCL_MessageStd( "number of sses: " + util::Format()( protein_model.GetNumberSSEs()));
      //BCL_MessageStd( "number of sses in chain a: " + util::Format()( protein_model.GetChain('A')->GetNumberSSEs()));
      //BCL_MessageStd( "number of sses in chain b: " + util::Format()( protein_model.GetChain('B')->GetNumberSSEs()));
      BCL_MessageStd( "number of aas in protein model: " + util::Format()( protein_model.GetNumberAAs()));

      // create a ShPtrVector of bodies which will be created from "bodies" and will be used to create the body
      // restraint
      util::ShPtr< util::ShPtrVector< assemble::SSEGeometryInterface> > restraint_bodies
      (
        new util::ShPtrVector< assemble::SSEGeometryInterface>()
      );

      // create coord::Bodies out of the SSEs of "sses" and insert the coord::Bodies into "bodies"
      for
      (
        util::SiPtrVector< const assemble::SSE>::iterator sse_itr( sses.Begin()), sse_itr_end( sses.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        util::ShPtr< assemble::SSE> current_body( ( *sse_itr)->Clone());

        BCL_MessageDbg
        (
          "sse aa size : " + util::Format()( current_body->GetSize()) + " and length : "
          + util::Format()( 2 * current_body->GetExtent( coord::GetAxes().e_Z))
        );

        // if the inserted body corresponded to a helix, then shorten its fragments
        if( current_body->GetType() == biol::GetSSTypes().HELIX)
        {
          BCL_MessageDbg
          (
            "restraint helix fragment size before change: " +
              util::Format()( current_body->GetGeometries().FirstElement()->GetExtent( coord::GetAxes().e_Z))
          );

          // reduce the size of the helix fragment so that end points fall into density rod
          current_body->SetFragmentGeometries( 5);

          BCL_MessageDbg
          (
            "restraint helix fragment size after change: " +
              util::Format()( current_body->GetGeometries().FirstElement()->GetExtent( coord::GetAxes().e_Z))
          );
        }

        // if the inserted body corresponded to a strand, then change its x and y extents and shorten its fragments
        else if( current_body->GetType() == biol::GetSSTypes().STRAND)
        {
          BCL_MessageDbg
          (
            "restraint strand fragment size before change: " +
              util::Format()( current_body->GetGeometries().FirstElement()->GetExtent( coord::GetAxes().e_Z))
          );

          // reduce the size of the strand fragment so that end points fall into density rod
          current_body->SetFragmentGeometries( 2);

          BCL_MessageDbg
          (
            "restraint strand fragment size after change: " +
            util::Format()( current_body->GetGeometries().FirstElement()->GetExtent( coord::GetAxes().e_Z))
          );

//             // set x and y to 1.6, leave z extent unchanged
//             // (this is necessary for the correct occupancy check, otherwise more than one body restraint could
//             // be occupied by the center of one sse, due to idealization))
//             bodies.LastElement().SetExtents
//             (
//               linal::Vector3D( 1.6, 1.6, bodies.LastElement().GetExtent( coord::GetAxes().e_Z))
//             );
        }
        restraint_bodies->PushBack( current_body);
      }

      BCL_MessageCrt
      (
        "number of body restraints created : " + util::Format()( restraint_bodies->GetSize())
      );

      return restraint_bodies;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_handler_data_set_pairwise_identifiers.h"
#include "restraint/bcl_restraint_locator_coordinates_first_side_chain_atom.h"

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
    const util::SiPtr< const util::ObjectInterface> HandlerDataSetPairwiseIdentifiers::s_Instance
    (
      GetObjectInstances().AddInstance( new HandlerDataSetPairwiseIdentifiers())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    HandlerDataSetPairwiseIdentifiers::HandlerDataSetPairwiseIdentifiers() :
      m_Score( util::GetUndefinedDouble()),
      m_DataSetPairwise()
    {
    }

    //! @brief constructor taking member variable parameters
    //! @param SCORE the score of the data set for reading or writing
    HandlerDataSetPairwiseIdentifiers::HandlerDataSetPairwiseIdentifiers( const double SCORE) :
      m_Score( SCORE),
      m_DataSetPairwise()
    {
    }

    //! @brief Clone function
    //! @return pointer to new HandlerDataSetPairwiseIdentifiers
    HandlerDataSetPairwiseIdentifiers *HandlerDataSetPairwiseIdentifiers::Clone() const
    {
      return new HandlerDataSetPairwiseIdentifiers( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &HandlerDataSetPairwiseIdentifiers::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the score
    //! @return double which is the score of the dataset
    const double HandlerDataSetPairwiseIdentifiers::GetScore() const
    {
      return m_Score;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief reads data set from a stream
    //! @param ISTREAM the stream the data set will be read from
    //! @return stream that the data set was read from
    std::istream &HandlerDataSetPairwiseIdentifiers::ReadDataSetPairwise( std::istream &ISTREAM)
    {
      std::string score_id;
      ISTREAM >> score_id >> m_Score;

      // read in the data set
      while( !ISTREAM.eof() && ISTREAM.peek() != std::istream::traits_type::eof())
      {
        DataPairwise data_pair
        (
          util::ShPtr< assemble::LocatorAtomCoordinatesInterface>( new LocatorCoordinatesFirstSideChainAtom()),
          util::ShPtr< assemble::LocatorAtomCoordinatesInterface>( new LocatorCoordinatesFirstSideChainAtom())
        );

        // read in current data pair
        data_pair.ReadIdentification( ISTREAM);

        // insert the data pair into the data set
        BCL_MessageDbg( "trying to insert " + data_pair.GetIdentification() + " into " + util::Format()( m_DataSetPairwise));
        BCL_Assert( m_DataSetPairwise.Insert( data_pair).second, "could not insert data pair " + data_pair.GetIdentification());
        BCL_MessageDbg( "successfully inserted " + data_pair.GetIdentification());
        // get next char
        char tmp;
        ISTREAM.get( tmp);
      }

      return ISTREAM;
    }

    //! @brief provides the data set that the handler created
    //! @return const reference to an data set object
    const DataSetPairwise &HandlerDataSetPairwiseIdentifiers::GetDataSetPairwise() const
    {
      return m_DataSetPairwise;
    }

    //! @brief writes restraints to a stream
    //! @param OSTREAM the stream the data set will be written to
    //! @param DATA_SET the data set that will be written to the stream
    //! @return stream the data set were written to
    std::ostream &HandlerDataSetPairwiseIdentifiers::WriteDataSetPairwise
    (
      std::ostream &OSTREAM, const DataSetPairwise &DATA_SET
    ) const
    {
      // print score
      OSTREAM << "score: " << m_Score << '\n';

      // iterate through the sorted data
      for
      (
        DataSetPairwise::const_iterator
        itr( DATA_SET.Begin()), itr_end( DATA_SET.End());
        itr != itr_end; ++itr
      )
      {
        OSTREAM << itr->GetIdentification() << '\n';
      }

      return OSTREAM;
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
    std::istream &HandlerDataSetPairwiseIdentifiers::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Score, ISTREAM);
      io::Serialize::Read( m_DataSetPairwise, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &HandlerDataSetPairwiseIdentifiers::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Score, OSTREAM, INDENT);
      io::Serialize::Write( m_DataSetPairwise, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_handler_data_set_pairwise_interface.h"
#include "util/bcl_util_class_descriptor.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &HandlerDataSetPairwiseInterface::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_handler_epr_decay.h"

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_hash_map.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> HandlerEPRDecay::s_Instance
    (
      util::Enumerated< HandlerBase< storage::Vector< EPRDecay> > >::AddInstance( new HandlerEPRDecay())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from defaults
    HandlerEPRDecay::HandlerEPRDecay()
    {
    }

    //! @brief Clone function
    //! @return pointer to new HandlerEPRDecay
    HandlerEPRDecay *HandlerEPRDecay::Clone() const
    {
      return new HandlerEPRDecay( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &HandlerEPRDecay::GetAlias() const
    {
      static const std::string s_name( "EPRDecay");
      return s_name;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer HandlerEPRDecay::GetSerializer() const
    {
      io::Serializer serializer( HandlerInterface::GetSerializer());
      serializer.SetClassDescription
      (
        "Reads files listing the results of EPR decay measurements. Lines that begin with # or ! are automatically ignored."
      );
      return serializer;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &HandlerEPRDecay::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief reads atom distance restraints from an input stream
    //! @brief input stream to read the restraints from
    //! @return the read in restraints
    storage::Vector< EPRDecay> HandlerEPRDecay::ReadRestraints( std::istream &ISTREAM) const
    {
      // create all lines
      storage::Vector< std::string> lines( util::StringLineListFromIStream( ISTREAM));

      // create a mapping from the spin-labeling site to the EPR decay data
      storage::HashMap< std::string, EPRDecay> epr_decay_map;

      // iterate over all lines and add each valid restraint to the corresponding EPRDecay instance
      for( auto line_itr( lines.Begin()), line_itr_end( lines.End()); line_itr != line_itr_end; ++line_itr)
      {
        // skip empty lines or comment lines indicated by a leading # or !
        if( line_itr->empty() || ( *line_itr)[ 0] == '!' || ( *line_itr)[ 0] == '#')
        {
          continue;
        }

        // split the line and convert the values
        storage::Vector< std::string> split_line( util::SplitString( *line_itr, " \t"));
        char chain_id_first, chain_id_second;
        int seq_id_first, seq_id_second;
        double time, decay;
        util::TryConvertFromString( chain_id_first, *split_line[ 0], util::GetLogger());
        util::TryConvertFromString( seq_id_first, *split_line[ 1], util::GetLogger());
        util::TryConvertFromString( chain_id_second, *split_line[ 2], util::GetLogger());
        util::TryConvertFromString( seq_id_second, *split_line[ 3], util::GetLogger());
        util::TryConvertFromString( time, *split_line[ 4], util::GetLogger());
        util::TryConvertFromString( decay, *split_line[ 5], util::GetLogger());

        // compute the unique key for this spin-labeling pair
        const std::string key( *split_line[ 0] + *split_line[ 1] + *split_line[ 2] + *split_line[ 3]);

        // if EPR decay data has already been added for this spin-labeling pair, add it to the same set
        if( epr_decay_map.Find( key) == epr_decay_map.End())
        {
          EPRDecay epr_decay_set( chain_id_first, seq_id_first, chain_id_second, seq_id_second);
          epr_decay_set.AddMeasurement( time, decay);
          epr_decay_map.Insert( storage::Pair< std::string, EPRDecay>( key, epr_decay_set));
        }
        else
        {
          epr_decay_map[ key].AddMeasurement( time, decay);
        }
      }

      // create the final list containing all EPR decay data sets
      storage::Vector< EPRDecay> decay_data;
      for( auto map_it( epr_decay_map.Begin()), map_it_end( epr_decay_map.End()); map_it != map_it_end; ++map_it)
      {
        decay_data.PushBack( map_it->second);
      }

      return decay_data;
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_handler_interface.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_dynamic.h"
#include "io/bcl_io_directory_entry.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    //! @brief get the restraints filename
    //! @return filename for this restraints; uses -restraint_prefix flag and GetFilenameExtension()
    std::string HandlerInterface::GetFilename() const
    {
      return GetFlagRestraintsFilePrefix()->GetFirstParameter()->GetValue() + GetFilenameExtension();
    }

    //! @brief test for the existence of the restraints file
    bool HandlerInterface::Exists() const
    {
      return io::DirectoryEntry( GetFilename()).DoesExist();
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer HandlerInterface::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        "Allows use of " + GetAlias() + " restraints. The restraints are read in from X."
        + m_DefaultExtension + ", where X is the path/prefix given to -restraint_prefix "
      );
      if( m_DefaultExtension.empty())
      {
        serializer.AddInitializer( "extension", "File extension for this restraint type. The restraints are read in from {path/prefix given to -restraint_prefix}.{extension}", io::Serialization::GetAgent( &m_Extension));
      }
      else
      {
        serializer.AddInitializer
        (
          "extension",
          "File extension, defaults to " + m_DefaultExtension
          + ". The restraints are read in from {path/prefix given to -restraint_prefix}.{extension}",
          io::Serialization::GetAgent( &m_Extension),
          ""
        );
      }
      return serializer;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_locator_coordinates_first_side_chain_atom.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> LocatorCoordinatesFirstSideChainAtom::s_Instance
    (
      util::Enumerated< find::LocatorCoordinatesInterface< assemble::ProteinModel> >::AddInstance
      (
        util::Enumerated< assemble::LocatorAtomCoordinatesInterface>::AddInstance
        (
          new LocatorCoordinatesFirstSideChainAtom()
        )
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LocatorCoordinatesFirstSideChainAtom::LocatorCoordinatesFirstSideChainAtom() :
      m_LocatorAA(),
      m_AtomType()
    {
    }

    //! @brief constructor from chain ID, amino acid ID
    //! @param CHAIN_ID char which indicates the chain
    //! @param SEQ_ID int which indicates the SeqID of the amino acid
    //! @param USE_PDB_ID true if the pdd id should be used to locate the residue instead of the seqid
    LocatorCoordinatesFirstSideChainAtom::LocatorCoordinatesFirstSideChainAtom( const char CHAIN_ID, const int SEQ_ID, bool USE_PDB_ID) :
      m_LocatorAA( CHAIN_ID, SEQ_ID, USE_PDB_ID),
      m_AtomType( biol::GetAtomTypes().GetUndefined())
    {
    }

    //! @brief constructor from chain ID, amino acid ID, and residue type
    //! @param CHAIN_ID char which indicates the chain
    //! @param SEQ_ID int which indicates the SeqID of the amino acid
    //! @param AA_TYPE AtomName which indicates the residue type
    //! @param USE_PDB_ID true if the pdd id should be used to locate the residue instead of the seqid
    LocatorCoordinatesFirstSideChainAtom::LocatorCoordinatesFirstSideChainAtom
    (
      const char CHAIN_ID, const int SEQ_ID, const biol::AAType &AA_TYPE, bool USE_PDB_ID
    ) :
      m_LocatorAA( CHAIN_ID, SEQ_ID, AA_TYPE, USE_PDB_ID),
      m_AtomType( biol::GetAtomTypes().GetUndefined())
    {
    }

    //! @brief constructor from chain ID, amino acid ID, and atom type
    //! @param CHAIN_ID char which indicates the chain
    //! @param SEQ_ID int which indicates the SeqID of the amino acid
    //! @param ATOM_TYPE indicates the atom type
    //! @param USE_PDB_ID true if the pdd id should be used to locate the residue instead of the seqid
    LocatorCoordinatesFirstSideChainAtom::LocatorCoordinatesFirstSideChainAtom
    (
      const char CHAIN_ID, const int SEQ_ID, const biol::AtomType &ATOM_TYPE, bool USE_PDB_ID
    ) :
      m_LocatorAA( CHAIN_ID, SEQ_ID, USE_PDB_ID),
      m_AtomType( ATOM_TYPE)
    {
    }

    //! @brief constructor from chain ID, amino acid ID, and atom type and residue type
    //! @param CHAIN_ID char which indicates the chain
    //! @param SEQ_ID int which indicates the SeqID of the amino acid
    //! @param ATOM_TYPE AtomName which indicates the atom type
    //! @param AA_TYPE AtomName which indicates the residue type
    //! @param USE_PDB_ID true if the pdd id should be used to locate the residue instead of the seqid
    LocatorCoordinatesFirstSideChainAtom::LocatorCoordinatesFirstSideChainAtom
    (
      const char CHAIN_ID, const int SEQ_ID, const biol::AtomType &ATOM_TYPE, const biol::AAType &AA_TYPE, bool USE_PDB_ID
    ) :
      m_LocatorAA( CHAIN_ID, SEQ_ID, AA_TYPE, USE_PDB_ID),
      m_AtomType( ATOM_TYPE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new LocatorCoordinatesFirstSideChainAtom
    LocatorCoordinatesFirstSideChainAtom *LocatorCoordinatesFirstSideChainAtom::Clone() const
    {
      return new LocatorCoordinatesFirstSideChainAtom( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LocatorCoordinatesFirstSideChainAtom::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &LocatorCoordinatesFirstSideChainAtom::GetAlias() const
    {
      static std::string s_alias( "FirstSideChainAtom");
      return s_alias;
    }

    //! @brief gives formatted string describing the locator
    //! @return formatted string describing the locator
    std::string LocatorCoordinatesFirstSideChainAtom::GetIdentification() const
    {
      return m_LocatorAA.GetIdentification();
    }

    //! @brief reads formatted string describing the locator
    //! @return formatted string describing the locator
    std::istream &LocatorCoordinatesFirstSideChainAtom::ReadIdentification( std::istream &ISTREAM)
    {
      m_LocatorAA.ReadIdentification( ISTREAM);
      return ISTREAM;
    }

    //! @brief gives the chain id the locator corresponds to
    //! @return the chain id the locator corresponds to
    char LocatorCoordinatesFirstSideChainAtom::GetChainID() const
    {
      return m_LocatorAA.GetLocatorChain().GetChainID();
    }

    //! @brief gives the seq id the locator corresponds to
    //! @return the seq id the locator corresponds to
    int LocatorCoordinatesFirstSideChainAtom::GetSeqID() const
    {
      return m_LocatorAA.GetAAID();
    }

    //! @brief gives the atom type the locator corresponds to
    //! @return the atom type the locator corresponds to
    const biol::AtomType &LocatorCoordinatesFirstSideChainAtom::GetAtomType() const
    {
      return m_AtomType;
    }

    //! @brief gives the aa type the locator corresponds to
    //! @return the aa type the locator corresponds to
    const biol::AAType &LocatorCoordinatesFirstSideChainAtom::GetAAType() const
    {
      return m_LocatorAA.GetAAType();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief locates an atom from a protein model
    //! @param PROTEIN_MODEL model from which the atom will be located
    //! @return siptr to atom that has been located from PROTEIN_MODEL
    util::SiPtr< const biol::Atom>
    LocatorCoordinatesFirstSideChainAtom::LocateAtom( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // locate the aa containing the atom
      const util::SiPtr< const biol::AABase> aa( m_LocatorAA.Locate( PROTEIN_MODEL));

      // true if the aa could not be located
      if( !aa.IsDefined())
      {
        // return undefined vector3d
        return util::SiPtr< const biol::Atom>();
      }

      // get the first side chain atom of the aa
      const util::SiPtr< const biol::Atom> atom( aa->GetFirstSidechainAtom());

      // return the coordinates
      return atom;
    }

    //! @brief locates an atom from a protein model
    //! @param PROTEIN_MODEL model from which the atom will be located
    //! @return copy of atom that has been located from PROTEIN_MODEL
    biol::Atom LocatorCoordinatesFirstSideChainAtom::LocateAtomCopy( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // try to locate atom
      const util::SiPtr< const biol::Atom> atom( LocateAtom( PROTEIN_MODEL));

      // true if atom could not be found
      if( !atom.IsDefined())
      {
        // return empty atom
        return biol::Atom();
      }

      // atom is defined so return atom behind pointer
      return *atom;
    }

    //! @brief locates the desired coordinates from a protein model
    //! @param PROTEIN_MODEL ProteinModel which the LocatorAtom refers to
    //! @return the coordinates of the atom denoted by the LocatorAtom
    linal::Vector3D LocatorCoordinatesFirstSideChainAtom::Locate( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // locate the aa containing the atom
      const util::SiPtr< const biol::AABase> aa( m_LocatorAA.Locate( PROTEIN_MODEL));

      // true if the aa could not be located
      if( !aa.IsDefined())
      {
        // return undefined vector3d
        return linal::Vector3D( util::GetUndefinedDouble());
      }

      // get the coordinates from the first side chain atom of the aa
      const linal::Vector3D coords( aa->GetFirstSidechainAtom().GetCoordinates());

      // return the coordinates
      return coords;
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
    std::istream &LocatorCoordinatesFirstSideChainAtom::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_LocatorAA, ISTREAM);
      io::Serialize::Read( m_AtomType,  ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &LocatorCoordinatesFirstSideChainAtom::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_LocatorAA, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AtomType,  OSTREAM, INDENT) << '\n';

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer LocatorCoordinatesFirstSideChainAtom::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "locates the coordinates of the first side chain atom of the indicated residue");
      parameters.AddInitializer
      (
        "locator",
        "locator for the residue of interest",
        io::Serialization::GetAgent( &m_LocatorAA)
      );
      parameters.AddOptionalInitializer
      (
        "atom_type",
        "atom type if the first side chain atom atom-type is known",
        io::Serialization::GetAgent( &m_AtomType)
      );
      return parameters;
    }

  } // namespace restraint

} // namespace bcl
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
#include "restraint/bcl_restraint_locator_coordinates_hydrogen.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_atom.h"
#include "biol/bcl_biol_aa_back_bone_completer.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically
#include <algorithm>
#include <cctype>
#include <functional>

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> LocatorCoordinatesHydrogen::s_Instance
    (
      util::Enumerated< find::LocatorCoordinatesInterface< assemble::ProteinModel> >::AddInstance
      (
        util::Enumerated< assemble::LocatorAtomCoordinatesInterface>::AddInstance
        (
          new LocatorCoordinatesHydrogen()
        )
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LocatorCoordinatesHydrogen::LocatorCoordinatesHydrogen() :
      m_ChainID( 'A'),
      m_SeqID( util::GetUndefined< int>()),
      m_AtomTypeString( biol::GetAtomTypes().e_Undefined),
      m_AtomType( biol::GetAtomTypes().e_Undefined)
    {
    }

    //! @brief construct from chain id, sequence id, and atom type
    //! @param CHAIN_ID chain id
    //! @param SEQ_ID sequence id
    //! @param ATOM_TYPE atom type string
    LocatorCoordinatesHydrogen::LocatorCoordinatesHydrogen
    (
      const char CHAIN_ID,
      const int SEQ_ID,
      const std::string &ATOM_TYPE
    ) :
      m_ChainID( CHAIN_ID),
      m_SeqID( SEQ_ID),
      m_AtomTypeString( ATOM_TYPE),
      m_AtomType( biol::GetAtomTypes().e_Undefined)
    {
      // make string upper case
      std::transform
      (
        m_AtomTypeString.begin(),
        m_AtomTypeString.end(),
        m_AtomTypeString.begin(),
        []( unsigned char c) { return std::toupper( c);}
      );

      // try to get the real atom type
      m_AtomType = GetAtomTypeFromString();
    }

    //! @brief Clone function
    //! @return pointer to new LocatorCoordinatesHydrogen
    LocatorCoordinatesHydrogen *LocatorCoordinatesHydrogen::Clone() const
    {
      return new LocatorCoordinatesHydrogen( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LocatorCoordinatesHydrogen::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &LocatorCoordinatesHydrogen::GetAlias() const
    {
      static std::string s_alias( "Hydrogen");
      return s_alias;
    }

    //! @brief gives formatted string describing the locator
    //! @return formatted string describing the locator
    std::string LocatorCoordinatesHydrogen::GetIdentification() const
    {
      return util::Format()( m_ChainID) + " " + util::Format()( m_SeqID) + " " + m_AtomTypeString;
    }

    //! @brief reads formatted string describing the locator
    //! @return formatted string describing the locator
    std::istream &LocatorCoordinatesHydrogen::ReadIdentification( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_ChainID, ISTREAM);
      io::Serialize::Read( m_SeqID, ISTREAM);
      io::Serialize::Read( m_AtomTypeString, ISTREAM);

      return ISTREAM;
    }

    //! @brief returns reference to undefined AA TYPE, this function is required by the interface but not used for
    //!        this class
    //! @return reference to undefined AA TYPE
    const biol::AAType &LocatorCoordinatesHydrogen::GetAAType() const
    {
      return biol::GetAATypes().e_Undefined;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief locates the desired coordinates from a protein model
    //! @param PROTEIN_MODEL ProteinModel which the LocatorAtom refers to
    //! @return the coordinates of the atom denoted by the LocatorAtom
    linal::Vector3D LocatorCoordinatesHydrogen::Locate( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      return LocateAtomCopy( PROTEIN_MODEL).GetCoordinates();
    }

    //! @brief locates an atom from a protein model, creating H or HA atoms as needed - other hydrogen atoms will have
    //!        m_AtomType as an atom type but with coords of the first side chain atom
    //! @param PROTEIN_MODEL model from which the atom will be located
    //! @return copy of atom that has been located from PROTEIN_MODEL
    biol::Atom LocatorCoordinatesHydrogen::LocateAtomCopy( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // initialize static undefined atom
      static const biol::Atom s_undefined_atom;

      // try to locate the residue that contains the atom
      const util::SiPtr< const biol::AABase> sp_residue( LocateAA( PROTEIN_MODEL));

      // if the residue was not located
      if( !sp_residue.IsDefined())
      {
        // return an undefined atom
        return s_undefined_atom;
      }

      // convert the string to an atom type
      biol::AtomType atom_type( m_AtomType.IsDefined() ? m_AtomType : GetAtomTypeFromString( sp_residue->GetType()));

      // if the atom type is undefined
      if( atom_type == biol::GetAtomTypes().e_Undefined)
      {
        // return undefined atom
        return s_undefined_atom;
      }

      // try to locate the atom
      const util::SiPtr< const biol::Atom> sp_atom
      (
        assemble::LocatorAtom::LocateAtomFromModel( PROTEIN_MODEL, m_ChainID, m_SeqID, atom_type)
      );

      // if the atom was located
      if( sp_atom.IsDefined() && sp_atom->GetCoordinates().IsDefined())
      {
        // return the atom
        return *sp_atom;
      }

      // if the atom to be located is not hydrogen
      if
      (
        atom_type->GetElementType() != chemistry::GetElementTypes().e_Hydrogen &&
        atom_type != biol::GetAtomTypes().O1)
      {
        // return a default atom
        return s_undefined_atom;
      }

      // if the atom type is H
      if( atom_type == biol::GetAtomTypes().H)
      {
        // locate the previous residue
        const assemble::LocatorAA aa_locator( m_ChainID, m_SeqID - 1);
        const util::SiPtr< const biol::AABase> sp_prev_aa( aa_locator.Locate( PROTEIN_MODEL));

        // if the previous AA was not located
        if( !sp_prev_aa.IsDefined())
        {
          // return an undefined atom
          return s_undefined_atom;
        }

        // generate an H atom
        return biol::AABackBoneCompleter::GenerateHydrogen
        (
          *sp_residue, sp_prev_aa->GetAtom( biol::GetAtomTypes().C).GetCoordinates()
        );
      }
      // if the atom type is HA
      else if( atom_type == biol::GetAtomTypes().HA || atom_type == biol::GetAtomTypes().HA3)
      {
        // generate an HA atom
        return biol::AABackBoneCompleter::GenerateHA( *sp_residue);
      }

      // return an atom with the first side chain atom coords but with the requested atom type in ATOM_LOCATOR
      return biol::Atom( sp_residue->GetFirstSidechainAtom().GetCoordinates(), atom_type);
    }

    //! @brief takes a string of the Atom type as well as the amino acid type and converts them to the
    //!        correct IUPAC type
    //! @param AA_TYPE the AAType of the atom from the restraint file
    //! @return the IUPAC formatted atom type
    biol::AtomType LocatorCoordinatesHydrogen::GetAtomTypeFromString( const biol::AAType &AA_TYPE) const
    {
      // initialize return type
      biol::AtomType atom_type
      (
        AA_TYPE.IsDefined() ?
          AA_TYPE->GetAtomTypeFromAtomName( m_AtomTypeString) :
          biol::GetAtomTypes().TypeFromPDBAtomName( m_AtomTypeString)
      );

      // if the atom type is defined
      if( atom_type.IsDefined())
      {
        // return it
        return atom_type;
      }

      // find the AA in the map
      const storage::Map< std::string, biol::AtomType> &atom_strings_map( GetPseudoAtomMap().GetValue( AA_TYPE));

      // check to see if the atom_string is in the map
      if( atom_strings_map.Has( m_AtomTypeString))
      {
        // return the type
        return atom_strings_map.GetValue( m_AtomTypeString);
      }

      // if atom is not in the map but has the *, # or % symbol after the atom type
      if
      (
        !m_AtomTypeString.empty() &&
        (
          m_AtomTypeString[ m_AtomTypeString.length() - 1] == '*' ||
          m_AtomTypeString[ m_AtomTypeString.length() - 1] == '#' ||
          m_AtomTypeString[ m_AtomTypeString.length() - 1] == '%'
        )
      )
      {
        // try replacing the last char
        atom_type = ReplacePseudoAtomChar( m_AtomTypeString, AA_TYPE);

        // return the type if it is defined
        if( atom_type.IsDefined())
        {
          return atom_type;
        }

        // iterate from 1 to 3
        for( char i( '1'); i != '4'; ++i)
        {
          // set the new string HG* to HG11
          std::string atom_string( m_AtomTypeString.substr( 0, m_AtomTypeString.length() - 1) + i + i);

          // try replacing the last char again
          atom_type = ReplacePseudoAtomChar( atom_string, AA_TYPE);

          // return the type if it is defined
          if( atom_type.IsDefined())
          {
            return atom_type;
          }
        }
      }

      // end
      return atom_type;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LocatorCoordinatesHydrogen::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ChainID, ISTREAM);
      io::Serialize::Read( m_SeqID, ISTREAM);
      io::Serialize::Read( m_AtomTypeString, ISTREAM);
      io::Serialize::Read( m_AtomType, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &LocatorCoordinatesHydrogen::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ChainID, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SeqID, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AtomTypeString, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AtomType, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief tries replacing the last character in an atom name string with 1, 2, or 3 and getting the atom type
    //! @param ATOM_NAME atom name string
    //! @param AA_TYPE the AAType of the atom from the restraint file
    //! @return atom type
    biol::AtomType LocatorCoordinatesHydrogen::ReplacePseudoAtomChar
    (
      const std::string &ATOM_NAME,
      const biol::AAType &AA_TYPE
    )
    {
      // set atom type to undefined
      biol::AtomType atom_type( biol::GetAtomTypes().e_Undefined);

      // copy the string
      std::string atom_string( ATOM_NAME);

      // if the AA type is not defined
      if( AA_TYPE == biol::GetAATypes().e_Undefined)
      {
        // change the last char to 1 and return
        atom_string[ atom_string.length() - 1] = '1';
        return biol::GetAtomTypes().TypeFromPDBAtomName( atom_string);
      }

      // iterate from 1 to 3
      for( char i( '1'); i != '4'; ++i)
      {
        // change the last character
        atom_string[ atom_string.length() - 1] = i;

        // try to get the type again
        atom_type = AA_TYPE->GetAtomTypeFromAtomName( atom_string);

        // if the atom_type is now defined
        if( atom_type.IsDefined())
        {
          // return it
          return atom_type;
        }
      }

      // return undefined
      return atom_type;
    }

    //! @brief gets a map of PsuedoAtoms that can be mapped back to an IUPAC type
    //! @return the map of these PsuedoAtoms
    const storage::Map
    <
      biol::AAType, storage::Map< std::string, biol::AtomType>
    > &LocatorCoordinatesHydrogen::GetPseudoAtomMap()
    {
      // create a map to store all of the data
      static storage::Map< biol::AAType, storage::Map< std::string, biol::AtomType> > s_map;
      if( s_map.IsEmpty())
      {
        // undefined (backbone or not residue-specific)
        storage::Map< std::string, biol::AtomType> undefined_data;
        undefined_data[ "HN"] = biol::GetAtomTypes().H;
        undefined_data[ "MTS"] = biol::GetAtomTypes().O1;
        undefined_data[ "MTSL"] = biol::GetAtomTypes().O1;
        s_map[ biol::GetAATypes().e_Undefined] = undefined_data;
        // GLY
        storage::Map< std::string, biol::AtomType> gly_data;
        gly_data[ "HN"] = biol::GetAtomTypes().H;
        gly_data[ "QPA"] = biol::GetAtomTypes().HA2;
        gly_data[ "QA"] = biol::GetAtomTypes().HA2;
        s_map[ biol::GetAATypes().GLY] = gly_data;
        // ALA
        storage::Map< std::string, biol::AtomType> ala_data;
        ala_data[ "HN"] = biol::GetAtomTypes().H;
        ala_data[ "MB"] = biol::GetAtomTypes().HB1;
        ala_data[ "QB"] = biol::GetAtomTypes().HB1;
        s_map[ biol::GetAATypes().ALA] = ala_data;
        // VAL
        storage::Map< std::string, biol::AtomType> val_data;
        val_data[ "HN"] = biol::GetAtomTypes().H;
        val_data[ "MG1"] = biol::GetAtomTypes().HG11;
        val_data[ "QG1"] = biol::GetAtomTypes().HG11;
        val_data[ "MG2"] = biol::GetAtomTypes().HG11;
        val_data[ "QG2"] = biol::GetAtomTypes().HG11;
        val_data[ "QQG"] = biol::GetAtomTypes().HG11;
        val_data[ "HG#"] = biol::GetAtomTypes().HG11;
        s_map[ biol::GetAATypes().VAL] = val_data;
        // ILE
        storage::Map< std::string, biol::AtomType> ile_data;
        ile_data[ "HN"] = biol::GetAtomTypes().H;
        ile_data[ "HG2"] = biol::GetAtomTypes().HG12;
        ile_data[ "HG3"] = biol::GetAtomTypes().HG13;
        ile_data[ "QPG"] = biol::GetAtomTypes().HG12;
        ile_data[ "QG1"] = biol::GetAtomTypes().HG12;
        ile_data[ "MG2"] = biol::GetAtomTypes().HG21;
        ile_data[ "QG2"] = biol::GetAtomTypes().HG21;
        ile_data[ "MG"] = biol::GetAtomTypes().HG21;
        ile_data[ "QG"] = biol::GetAtomTypes().HG21;
        ile_data[ "MD"] = biol::GetAtomTypes().HG11;
        ile_data[ "QD1"] = biol::GetAtomTypes().HG11;
        ile_data[ "HD#"] = biol::GetAtomTypes().HD11;
        s_map[ biol::GetAATypes().ILE] = ile_data;
        // LEU
        storage::Map< std::string, biol::AtomType> leu_data;
        leu_data[ "HN"] = biol::GetAtomTypes().H;
        leu_data[ "QPB"] = biol::GetAtomTypes().HB2;
        leu_data[ "QB"] = biol::GetAtomTypes().HB2;
        leu_data[ "MD1"] = biol::GetAtomTypes().HD11;
        leu_data[ "QD1"] = biol::GetAtomTypes().HD11;
        leu_data[ "MD2"] = biol::GetAtomTypes().HD21;
        leu_data[ "QD2"] = biol::GetAtomTypes().HD21;
        leu_data[ "CD1"] = biol::GetAtomTypes().CD1;
        leu_data[ "QQD"] = biol::GetAtomTypes().HD11;
        leu_data[ "MDX"] = biol::GetAtomTypes().HD11;
        leu_data[ "HD#"] = biol::GetAtomTypes().HD11;
        s_map[ biol::GetAATypes().LEU] = leu_data;
        // PHE
        storage::Map< std::string, biol::AtomType> phe_data;
        phe_data[ "HN"] = biol::GetAtomTypes().H;
        phe_data[ "QPB"] = biol::GetAtomTypes().HB2;
        phe_data[ "QB"] = biol::GetAtomTypes().HB2;
        phe_data[ "HD"] = biol::GetAtomTypes().HD1;
        phe_data[ "CG"] = biol::GetAtomTypes().HD1;
        phe_data[ "HE"] = biol::GetAtomTypes().HE1;
        phe_data[ "CZ"] = biol::GetAtomTypes().HE1;
        phe_data[ "HZ"] = biol::GetAtomTypes().HZ;
        phe_data[ "QD"] = biol::GetAtomTypes().HD1;
        phe_data[ "QE"] = biol::GetAtomTypes().HE1;
        s_map[ biol::GetAATypes().PHE] = phe_data;
        // PRO
        storage::Map< std::string, biol::AtomType> pro_data;
        pro_data[ "QPB"] = biol::GetAtomTypes().HB2;
        pro_data[ "QB"] = biol::GetAtomTypes().HB2;
        pro_data[ "QPG"] = biol::GetAtomTypes().HG2;
        pro_data[ "QG"] = biol::GetAtomTypes().HG2;
        pro_data[ "QPD"] = biol::GetAtomTypes().HD2;
        pro_data[ "QD"] = biol::GetAtomTypes().HD2;
        s_map[ biol::GetAATypes().PRO] = pro_data;
        // MET
        storage::Map< std::string, biol::AtomType> met_data;
        met_data[ "HN"] = biol::GetAtomTypes().H;
        met_data[ "QPB"] = biol::GetAtomTypes().HB2;
        met_data[ "QB"] = biol::GetAtomTypes().HB2;
        met_data[ "QPG"] = biol::GetAtomTypes().HG2;
        met_data[ "QG"] = biol::GetAtomTypes().HG2;
        met_data[ "ME"] = biol::GetAtomTypes().HE1;
        met_data[ "QE"] = biol::GetAtomTypes().HE1;
        s_map[ biol::GetAATypes().MET] = met_data;
        // TRP
        storage::Map< std::string, biol::AtomType> trp_data;
        trp_data[ "HN"] = biol::GetAtomTypes().H;
        trp_data[ "QPB"] = biol::GetAtomTypes().HB2;
        trp_data[ "QB"] = biol::GetAtomTypes().HB2;
        trp_data[ "HNE"] = biol::GetAtomTypes().HE1;
        trp_data[ "HE1"] = biol::GetAtomTypes().HE1;
        trp_data[ "HD"] = biol::GetAtomTypes().HD1;
        trp_data[ "HE"] = biol::GetAtomTypes().HE3;
        trp_data[ "HH"] = biol::GetAtomTypes().HH2;
        trp_data[ "N1"] = biol::GetAtomTypes().HE1;
        trp_data[ "N1H"] = biol::GetAtomTypes().HE1;
        trp_data[ "H2"] = biol::GetAtomTypes().HD1;
        trp_data[ "C2"] = biol::GetAtomTypes().HD1;
        trp_data[ "C2H"] = biol::GetAtomTypes().HD1;
        trp_data[ "H4"] = biol::GetAtomTypes().HE3;
        trp_data[ "C4"] = biol::GetAtomTypes().HE3;
        trp_data[ "C4H"] = biol::GetAtomTypes().HE3;
        trp_data[ "H5"] = biol::GetAtomTypes().HZ3;
        trp_data[ "C5"] = biol::GetAtomTypes().HZ3;
        trp_data[ "C5H"] = biol::GetAtomTypes().HZ3;
        trp_data[ "H6"] = biol::GetAtomTypes().HH2;
        trp_data[ "C6"] = biol::GetAtomTypes().HH2;
        trp_data[ "C6H"] = biol::GetAtomTypes().HH2;
        trp_data[ "H7"] = biol::GetAtomTypes().HZ2;
        trp_data[ "C7"] = biol::GetAtomTypes().HZ2;
        trp_data[ "C7H"] = biol::GetAtomTypes().HZ2;
        s_map[ biol::GetAATypes().TRP] = trp_data;
        // SER
        storage::Map< std::string, biol::AtomType> ser_data;
        ser_data[ "HN"] = biol::GetAtomTypes().H;
        ser_data[ "QPB"] = biol::GetAtomTypes().HB2;
        ser_data[ "QB"] = biol::GetAtomTypes().HB2;
        ser_data[ "OG"] = biol::GetAtomTypes().HG;
        ser_data[ "HOG"] = biol::GetAtomTypes().HG;
        s_map[ biol::GetAATypes().SER] = ser_data;
        // THR
        storage::Map< std::string, biol::AtomType> thr_data;
        thr_data[ "HN"] = biol::GetAtomTypes().H;
        thr_data[ "HOG"] = biol::GetAtomTypes().HG1;
        thr_data[ "MG"] = biol::GetAtomTypes().HG21;
        thr_data[ "QG"] = biol::GetAtomTypes().HG21;
        thr_data[ "MG2"] = biol::GetAtomTypes().HG21;
        thr_data[ "QG2"] = biol::GetAtomTypes().HG21;
        s_map[ biol::GetAATypes().THR] = thr_data;
        // ASN
        storage::Map< std::string, biol::AtomType> asn_data;
        asn_data[ "HN"] = biol::GetAtomTypes().H;
        asn_data[ "QPB"] = biol::GetAtomTypes().HB2;
        asn_data[ "QB"] = biol::GetAtomTypes().HB2;
        asn_data[ "HND1"] = biol::GetAtomTypes().HD21;
        asn_data[ "HND2"] = biol::GetAtomTypes().HD22;
        asn_data[ "ND2"] = biol::GetAtomTypes().HD21;
        asn_data[ "HND"] = biol::GetAtomTypes().HD21;
        asn_data[ "QD2"] = biol::GetAtomTypes().HD21;
        s_map[ biol::GetAATypes().ASN] = asn_data;
        // GLN
        storage::Map< std::string, biol::AtomType> gln_data;
        gln_data[ "HN"] = biol::GetAtomTypes().H;
        gln_data[ "QPB"] = biol::GetAtomTypes().HB2;
        gln_data[ "QB"] = biol::GetAtomTypes().HB2;
        gln_data[ "QPG"] = biol::GetAtomTypes().HG2;
        gln_data[ "QG"] = biol::GetAtomTypes().HG2;
        gln_data[ "HNE1"] = biol::GetAtomTypes().HE21;
        gln_data[ "HNE2"] = biol::GetAtomTypes().HE22;
        gln_data[ "HE21"] = biol::GetAtomTypes().HE21;
        gln_data[ "HE22"] = biol::GetAtomTypes().HE22;
        gln_data[ "NE2"] = biol::GetAtomTypes().HE21;
        gln_data[ "HNE"] = biol::GetAtomTypes().HE21;
        gln_data[ "QE2"] = biol::GetAtomTypes().HE21;
        s_map[ biol::GetAATypes().GLN] = gln_data;
        // TYR
        storage::Map< std::string, biol::AtomType> tyr_data;
        tyr_data[ "HN"] = biol::GetAtomTypes().H;
        tyr_data[ "QPB"] = biol::GetAtomTypes().HB2;
        tyr_data[ "QB"] = biol::GetAtomTypes().HB2;
        tyr_data[ "HD"] = biol::GetAtomTypes().HD1;
        tyr_data[ "CG"] = biol::GetAtomTypes().HD1;
        tyr_data[ "HE"] = biol::GetAtomTypes().HE1;
        tyr_data[ "CZ"] = biol::GetAtomTypes().HE1;
        tyr_data[ "HOH"] = biol::GetAtomTypes().HH;
        tyr_data[ "QD"] = biol::GetAtomTypes().HD1;
        tyr_data[ "QE"] = biol::GetAtomTypes().HE1;
        s_map[ biol::GetAATypes().TYR] = tyr_data;
        // HIS
        storage::Map< std::string, biol::AtomType> his_data;
        his_data[ "HN"] = biol::GetAtomTypes().H;
        his_data[ "QPB"] = biol::GetAtomTypes().HB2;
        his_data[ "QB"] = biol::GetAtomTypes().HB2;
        his_data[ "HND"] = biol::GetAtomTypes().HD1;
        his_data[ "HE"] = biol::GetAtomTypes().HE1;
        his_data[ "HNE"] = biol::GetAtomTypes().HE2;
        his_data[ "HD"] = biol::GetAtomTypes().HD2;
        s_map[ biol::GetAATypes().HIS] = his_data;
        // ASP
        storage::Map< std::string, biol::AtomType> asp_data;
        asp_data[ "HN"] = biol::GetAtomTypes().H;
        asp_data[ "QPB"] = biol::GetAtomTypes().HB2;
        asp_data[ "QB"] = biol::GetAtomTypes().HB2;
        s_map[ biol::GetAATypes().ASP] = asp_data;
        // GLU
        storage::Map< std::string, biol::AtomType> glu_data;
        glu_data[ "HN"] = biol::GetAtomTypes().H;
        glu_data[ "QPB"] = biol::GetAtomTypes().HB2;
        glu_data[ "QB"] = biol::GetAtomTypes().HB2;
        glu_data[ "QPG"] = biol::GetAtomTypes().HG2;
        glu_data[ "QG"] = biol::GetAtomTypes().HG2;
        s_map[ biol::GetAATypes().GLU] = glu_data;
        // LYS
        storage::Map< std::string, biol::AtomType> lys_data;
        lys_data[ "HN"] = biol::GetAtomTypes().H;
        lys_data[ "QPB"] = biol::GetAtomTypes().HB2;
        lys_data[ "QB"] = biol::GetAtomTypes().HB2;
        lys_data[ "QPG"] = biol::GetAtomTypes().HG2;
        lys_data[ "QG"] = biol::GetAtomTypes().HG2;
        lys_data[ "QPD"] = biol::GetAtomTypes().HD2;
        lys_data[ "QD"] = biol::GetAtomTypes().HD2;
        lys_data[ "QPE"] = biol::GetAtomTypes().HE2;
        lys_data[ "QE"] = biol::GetAtomTypes().HE2;
        lys_data[ "NZ"] = biol::GetAtomTypes().HZ1;
        lys_data[ "HNZ"] = biol::GetAtomTypes().HZ1;
        lys_data[ "MNZ"] = biol::GetAtomTypes().HZ1;
        lys_data[ "QNZ"] = biol::GetAtomTypes().HZ1;
        s_map[ biol::GetAATypes().LYS] = lys_data;
        // ARG
        storage::Map< std::string, biol::AtomType> arg_data;
        arg_data[ "HN"] = biol::GetAtomTypes().H;
        arg_data[ "QPB"] = biol::GetAtomTypes().HB2;
        arg_data[ "QB"] = biol::GetAtomTypes().HB2;
        arg_data[ "QPG"] = biol::GetAtomTypes().HG2;
        arg_data[ "QG"] = biol::GetAtomTypes().HG2;
        arg_data[ "QPD"] = biol::GetAtomTypes().HD2;
        arg_data[ "QD"] = biol::GetAtomTypes().HD2;
        arg_data[ "CD"] = biol::GetAtomTypes().CD;
        arg_data[ "NE"] = biol::GetAtomTypes().HE;
        arg_data[ "HNE"] = biol::GetAtomTypes().HE;
        arg_data[ "HE"] = biol::GetAtomTypes().HE;
        arg_data[ "NH1"] = biol::GetAtomTypes().HH11;
        arg_data[ "NH2"] = biol::GetAtomTypes().HH21;
        s_map[ biol::GetAATypes().ARG] = arg_data;
        // CYS
        storage::Map< std::string, biol::AtomType> cys_data;
        cys_data[ "HN"] = biol::GetAtomTypes().H;
        cys_data[ "QPB"] = biol::GetAtomTypes().HB2;
        cys_data[ "QB"] = biol::GetAtomTypes().HB2;
        cys_data[ "HSG"] = biol::GetAtomTypes().HG;
        s_map[ biol::GetAATypes().CYS] = cys_data;
      }
      return s_map;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer LocatorCoordinatesHydrogen::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Locates atoms/coordinates, including hydrogens, from a protein model\n"
        "Handles locating hydrogen atoms since they are not typically included in BCL protein models but are\n"
        "often used in experimental restraints, most notably, NMR.  Classes that use this class may elect to locate\n"
        "the hydrogens (using LocateAtomCopy), otherwise, the CB position will be returned"
      );
      parameters.AddInitializer
      (
        "chain",
        "the desired chain id",
        io::Serialization::GetAgentWithRange( &m_ChainID, 'A', 'Z')
      );
      parameters.AddInitializer
      (
        "sequence",
        "the sequence id of the desired residue",
        io::Serialization::GetAgentWithMin( &m_SeqID, 0)
      );
      parameters.AddInitializer
      (
        "atom type",
        "atom type of interest, if known",
        io::Serialization::GetAgent( &m_AtomTypeString)
      );
      parameters.AddDataMember
      (
        "actual_atom_type",
        io::Serialization::GetAgent( &m_AtomType)
      );
      return parameters;
    }

    //! @brief hook that derived classes can override, called after TryRead successfully reads a serializer containing only initializer info (no data variables)
    //! @param SERIALIZER the serializer object with initialization information
    //! @param ERR_STREAM stream to write out errors to
    //! @return true, unless there were new errors
    bool LocatorCoordinatesHydrogen::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &SERIALIZER,
      std::ostream &ERR_STREAM
    )
    {
      // make string upper case
      std::transform
      (
        m_AtomTypeString.begin(),
        m_AtomTypeString.end(),
        m_AtomTypeString.begin(),
        []( unsigned char c) { return std::toupper( c);}
      );

      // try to get the real atom type
      m_AtomType = GetAtomTypeFromString();
      return true;
    }

  } // namespace restraint

} // namespace bcl
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
#include "restraint/bcl_restraint_mutate_data_set_pairwise_add.h"

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
    const util::SiPtr< const util::ObjectInterface> MutateDataSetPairwiseAdd::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateDataSetPairwiseAdd())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &MutateDataSetPairwiseAdd::GetDefaultScheme()
    {
      static const std::string s_scheme( "add");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from optional scheme
    //! @param SCHEME optional scheme
    MutateDataSetPairwiseAdd::MutateDataSetPairwiseAdd( const std::string &SCHEME) :
      m_CompleteDataSet(),
      m_SizeRange(),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking sequence
    //! @param POOL_DATA_SET pool of possible data points to add
    //! @param SCHEME optional scheme
    MutateDataSetPairwiseAdd::MutateDataSetPairwiseAdd
    (
      const util::ShPtr< DataSetPairwise> &POOL_DATA_SET, const std::string &SCHEME
    ) :
      m_CompleteDataSet( POOL_DATA_SET),
      m_SizeRange( 1, 1),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking sequence
    //! @param POOL_DATA_SET pool of possible data points to add
    //! @param MIN min possible number that will be added
    //! @param MAX max possible number that will be added
    //! @param SCHEME optional scheme
    MutateDataSetPairwiseAdd::MutateDataSetPairwiseAdd
    (
      const util::ShPtr< DataSetPairwise> &POOL_DATA_SET, const size_t MIN, const size_t MAX, const std::string &SCHEME
    ) :
      m_CompleteDataSet( POOL_DATA_SET),
      m_SizeRange( MIN, MAX),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateDataSetPairwiseAdd
    MutateDataSetPairwiseAdd *MutateDataSetPairwiseAdd::Clone() const
    {
      return new MutateDataSetPairwiseAdd( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateDataSetPairwiseAdd::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this mutate
    //! @return the scheme for this mutate
    const std::string &MutateDataSetPairwiseAdd::GetScheme() const
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
    MutateDataSetPairwiseAdd::operator()( const DataSetPairwise &DATA_SET) const
    {
//      // will hold the data points that are available to add to to the current dataset
//      std::vector< DataPairwise> difference;
//
//      BCL_Assert( m_CompleteDataSet.IsDefined(), "complete data set pointer is not defined");
//
//      std::set_difference
//      (
//        m_CompleteDataSet->Begin(), m_CompleteDataSet->End(),
//        DATA_SET.Begin(), DATA_SET.End(),
//        std::inserter( difference, difference.begin())
//      );

      // static empty data set
      static util::ShPtr< DataSetPairwise> s_empty_data_set;

      // true if no available data pairs to add
      if( m_CompleteDataSet->IsEmpty()) //difference.empty())
      {
        BCL_MessageDbg( "difference is empty");
        // return skipped move
        return math::MutateResult< DataSetPairwise>( s_empty_data_set, *this);
      }

      util::ShPtr< DataSetPairwise> new_data_set( DATA_SET.Clone());

//      // randomly shuffle the difference data pairs
//      std::random_shuffle( difference.begin(), difference.end());

      // get random number of elements to add within range but no more than are in difference
      const size_t num_to_add( std::min( random::GetGlobalRandom().Random( m_SizeRange), m_CompleteDataSet->GetSize())); // difference.size()));

      // add as many elements are necessary to the new data set
      for( size_t num_added( 0); num_added < num_to_add; ++num_added)
      {
//        // add the as many elements as needed
//        new_data_set->InsertElements( difference.begin(), difference.begin() + num_to_add);

        // get random iterator to one of the data pairs in the pool of possibles
        DataSetPairwise::const_iterator random_itr
        (
          random::GetGlobalRandom().Iterator
          (
            m_CompleteDataSet->Begin(), m_CompleteDataSet->End(), m_CompleteDataSet->GetSize()
          )
        );
        BCL_Assert( random_itr != m_CompleteDataSet->End(), "random iterator at end of total dataset");

        // insert random data pair from the complete data set into the new data set
        new_data_set->Insert( *random_itr);
      }
      BCL_MessageDbg
      (
        "MutateDataSetPairwiseAdd data set size is " + util::Format()( new_data_set->GetSize())
      );
      BCL_Assert( new_data_set.IsDefined(), "new_data_set pointer is not defined");
      return math::MutateResult< DataSetPairwise>( new_data_set, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateDataSetPairwiseAdd::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_CompleteDataSet, ISTREAM);
      io::Serialize::Read( m_SizeRange, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateDataSetPairwiseAdd::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_CompleteDataSet, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SizeRange, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint
  
} // namespace bcl
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
#include "restraint/bcl_restraint_mutate_data_set_pairwise_filter_aa_type.h"

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
    const util::SiPtr< const util::ObjectInterface> MutateDataSetPairwiseFilterAAType::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateDataSetPairwiseFilterAAType())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &MutateDataSetPairwiseFilterAAType::GetDefaultScheme()
    {
      static const std::string s_scheme( "filter_aa_type_excl");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateDataSetPairwiseFilterAAType::MutateDataSetPairwiseFilterAAType( const std::string &SCHEME) :
      m_AATypesToExclude(),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking scheme
    //! @param AA_TYPES the types of amino acids which are undesirable
    //! @param SCHEME the scheme for this scoring function
    MutateDataSetPairwiseFilterAAType::MutateDataSetPairwiseFilterAAType
    (
      const storage::Set< biol::AAType> AA_TYPES, const std::string &SCHEME
    ) :
      m_AATypesToExclude( AA_TYPES),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateDataSetPairwiseFilterAAType
    MutateDataSetPairwiseFilterAAType *MutateDataSetPairwiseFilterAAType::Clone() const
    {
      return new MutateDataSetPairwiseFilterAAType( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateDataSetPairwiseFilterAAType::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &MutateDataSetPairwiseFilterAAType::GetScheme() const
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
    //! @param SEQUENCE DataSetPairwise of interest that will be mutated
    //! @return MutateResult that results from mutating to the argument DATA_SET
    math::MutateResult< DataSetPairwise>
    MutateDataSetPairwiseFilterAAType::operator()( const DataSetPairwise &SEQUENCE) const
    {
      static util::ShPtr< DataSetPairwise> s_empty;

      if( SEQUENCE.IsEmpty() || m_AATypesToExclude.IsEmpty())
      {
        BCL_MessageStd
        (
          "data set or AATypes to exclude is empty. data set size : " +
          util::Format()( SEQUENCE.GetSize()) + " AATypes to exclude size is : " +
          util::Format()( m_AATypesToExclude.GetSize())
        );
        return math::MutateResult< DataSetPairwise>( s_empty, *this);
      }

      // the dataset filtered
      util::ShPtr< DataSetPairwise> filtered_data( new DataSetPairwise());

      // iterate through the data set to filter out data pairs that don't meet exposure criteria
      for
      (
        DataSetPairwise::const_iterator data_itr( SEQUENCE.Begin()), data_itr_end( SEQUENCE.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // get the aa type of the locators
        const biol::AAType &aa_type_a( data_itr->First()->GetAAType());
        const biol::AAType &aa_type_b( data_itr->Second()->GetAAType());

        // true if neither of the aa types are in m_AATypesToExclude
        if( !m_AATypesToExclude.Contains( aa_type_a) && !m_AATypesToExclude.Contains( aa_type_b))
        {
          // add data pair to filtered_data
          filtered_data->Insert( *data_itr);
        }
      }

      BCL_MessageStd
      (
        "MutateDataSetPairwiseFilterAAType num unique locators " +
        util::Format()( filtered_data->GetUniqueDataPoints().GetSize())
      );

      BCL_MessageStd( "end MutateDataSetPairwiseFilterAAType::operator()");
      return math::MutateResult< DataSetPairwise>( filtered_data, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateDataSetPairwiseFilterAAType::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_AATypesToExclude, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateDataSetPairwiseFilterAAType::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_AATypesToExclude, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint
  
} // namespace bcl
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
#include "math/bcl_math_mutate_result.h"
#include "restraint/bcl_restraint_data_set_pairwise.h"
#include "restraint/bcl_restraint_mutate_data_set_pairwise_filter_coordinate_exclusion.h"
#include "score/bcl_score_data_set_pairwise_coordinate_exclusion.h"

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
    const util::SiPtr< const util::ObjectInterface> MutateDataSetPairwiseFilterCoordinateExclusion::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateDataSetPairwiseFilterCoordinateExclusion())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &MutateDataSetPairwiseFilterCoordinateExclusion::GetDefaultScheme()
    {
      static const std::string s_scheme( "filter_coord_exclusion");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateDataSetPairwiseFilterCoordinateExclusion::MutateDataSetPairwiseFilterCoordinateExclusion( const std::string &SCHEME) :
        m_Radius(),
        m_DistanceMap(),
        m_Scheme( SCHEME)
      {
      }

    //! @brief constructor taking members
    //! @param EXCLUSION_RADIUS data with atoms inside this radius of any exclusion coordinates will be penalized
    //! @param X_COORD_COLUMN column in istream that has the x-coordinate - columns start at 0
    //! @param Y_COORD_COLUMN column in istream that has the y-coordinate - columns start at 0
    //! @param Z_COORD_COLUMN column in istream that has the z-coordinate - columns start at 0
    //! @param ALL_POSSIBLE_DATA_POINTS the set of data points that should be subjected to this filter
    //! @param ENSEMBLE ensemble for which coordinates will be checked to make sure they aren't near exclusion coords
    //! @param SCHEME the scheme for this scoring function
    MutateDataSetPairwiseFilterCoordinateExclusion::MutateDataSetPairwiseFilterCoordinateExclusion
    (
      const double &EXCLUSION_RADIUS,
      std::istream &READ,
      const size_t X_COORD_COLUMN,
      const size_t Y_COORD_COLUMN,
      const size_t Z_COORD_COLUMN,
      const storage::Set
      <
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>,
        assemble::LocatorAtomCoordinatesInterface::PtrLessThan
      > &ALL_POSSIBLE_DATA_POINTS,
      const assemble::ProteinEnsemble &ENSEMBLE,
      const std::string &SCHEME
    ) :
      m_Radius( EXCLUSION_RADIUS),
      m_DistanceMap(),
      m_Scheme( SCHEME)
    {
      m_DistanceMap = score::DataSetPairwiseCoordinateExclusion::FillDistanceMap
        ( READ, X_COORD_COLUMN, Y_COORD_COLUMN, Z_COORD_COLUMN, ALL_POSSIBLE_DATA_POINTS, ENSEMBLE);
    }

    //! @brief Clone function
    //! @return pointer to new MutateDataSetPairwiseFilterCoordinateExclusion
    MutateDataSetPairwiseFilterCoordinateExclusion *MutateDataSetPairwiseFilterCoordinateExclusion::Clone() const
    {
      return new MutateDataSetPairwiseFilterCoordinateExclusion( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateDataSetPairwiseFilterCoordinateExclusion::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &MutateDataSetPairwiseFilterCoordinateExclusion::GetScheme() const
    {
      return m_Scheme;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an DATA and returning a mutated object of DATA
    //! @param DATA DATA of interest that will be mutated
    //! @return MutateResult that results from mutating to the argument DATA
    math::MutateResult< DataSetPairwise>
    MutateDataSetPairwiseFilterCoordinateExclusion::operator()( const DataSetPairwise &DATA) const
    {
      static util::ShPtr< DataSetPairwise> s_empty;

      if( DATA.IsEmpty())
      {
        BCL_MessageStd( "MutateDataSetPairwiseFilterCoordinateExclusion data set is empty");
        return math::MutateResult< DataSetPairwise>( s_empty, *this);
      }

      // the dataset filtered
      util::ShPtr< DataSetPairwise> filtered_data( new DataSetPairwise());

      // iterate through the sorted dataset
      for
      (
        DataSetPairwise::const_iterator data_itr( DATA.Begin()), data_itr_end( DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        const double score_first
        (
          score::DataSetPairwiseCoordinateExclusion::CalculateExclusionScore( data_itr->First(), m_DistanceMap, m_Radius)
        );
        const double score_second
        (
          score::DataSetPairwiseCoordinateExclusion::CalculateExclusionScore( data_itr->Second(), m_DistanceMap, m_Radius)
        );

        // true if both the scores are defined and zero - meaning they meet the radius cutoff around exclusion coords
        if
        (
          util::IsDefined( score_first) && util::IsDefined( score_second) &&
          score_first == double( 0) && score_first == score_second
        )
        {
          // add data pair to filtered_data
          filtered_data->Insert( *data_itr);
        }
      }

      return math::MutateResult< DataSetPairwise>( filtered_data, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateDataSetPairwiseFilterCoordinateExclusion::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Radius, ISTREAM);
      io::Serialize::Read( m_DistanceMap, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateDataSetPairwiseFilterCoordinateExclusion::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Radius, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DistanceMap, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint
  
} // namespace bcl
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
#include "restraint/bcl_restraint_mutate_data_set_pairwise_filter_euclidian_distance.h"

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
    const util::SiPtr< const util::ObjectInterface> MutateDataSetPairwiseFilterEuclidianDistance::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateDataSetPairwiseFilterEuclidianDistance())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &MutateDataSetPairwiseFilterEuclidianDistance::GetDefaultScheme()
    {
      static const std::string s_scheme( "distance_range");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor taking optional scheme
    //! @param SCHEME the scheme for this scoring function
    MutateDataSetPairwiseFilterEuclidianDistance::MutateDataSetPairwiseFilterEuclidianDistance
    (
      const std::string &SCHEME
    ) :
      m_DistanceRange(),
      m_Ensemble(),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking optional scheme
    //! @param RANGE desired range of euclidian separation between data points
    //! @param ENSEMBLE ensemble distances will be calculated from
    //! @param SCHEME the scheme for this scoring function
    MutateDataSetPairwiseFilterEuclidianDistance::MutateDataSetPairwiseFilterEuclidianDistance
    (
      const math::Range< double> &RANGE,
      const util::ShPtr< assemble::ProteinEnsemble> &ENSEMBLE,
      const std::string &SCHEME
    ) :
      m_DistanceRange( RANGE),
      m_Ensemble( ENSEMBLE),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateDataSetPairwiseFilterEuclidianDistance
    MutateDataSetPairwiseFilterEuclidianDistance *MutateDataSetPairwiseFilterEuclidianDistance::Clone() const
    {
      return new MutateDataSetPairwiseFilterEuclidianDistance( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateDataSetPairwiseFilterEuclidianDistance::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &MutateDataSetPairwiseFilterEuclidianDistance::GetScheme() const
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
    //! @param DATA DataSetPairwise of interest that will be mutated
    //! @return MutateResult that results from mutating to the argument DATA_SET
    math::MutateResult< DataSetPairwise>
    MutateDataSetPairwiseFilterEuclidianDistance::operator()( const DataSetPairwise &DATA) const
    {

      static util::ShPtr< DataSetPairwise> s_empty;

      if( DATA.IsEmpty())
      {
        BCL_MessageDbg( "MutateDataSetPairwiseFilterEuclidianDistance data set is empty");
        return math::MutateResult< DataSetPairwise>( s_empty, *this);
      }

      // the data set filtered
      util::ShPtr< DataSetPairwise> filtered_data( new DataSetPairwise());

      // iterate through the data set to filter out data pairs that don't meet exposure criteria
      for
      (
        DataSetPairwise::const_iterator data_itr( DATA.Begin()), data_itr_end( DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // get the average distance of the current data pair in m_Ensemble
        const math::RunningAverageSD< double> mean_sd( data_itr->EuclidianDistance( *m_Ensemble).First());

        // if there are counts get the mean, otherwise mean is undefined
        const double mean_distance( mean_sd.GetWeight() ? mean_sd.GetAverage() : util::GetUndefinedDouble());

        // true if size is within the desired range and defined
        if( m_DistanceRange.IsWithin( mean_distance) && util::IsDefined( mean_distance))
        {
          // add data pair to filtered_data
          filtered_data->Insert( *data_itr);
        }

      }

      BCL_MessageDbg
      (
        "MutateDataSetPairwiseFilterEuclidianDistance data set size " +
        util::Format()( filtered_data->GetSize())
      );

      return math::MutateResult< DataSetPairwise>( filtered_data, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateDataSetPairwiseFilterEuclidianDistance::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_DistanceRange, ISTREAM);
      io::Serialize::Read( m_Ensemble, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateDataSetPairwiseFilterEuclidianDistance::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_DistanceRange, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Ensemble, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint
  
} // namespace bcl
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
#include "restraint/bcl_restraint_mutate_data_set_pairwise_filter_exposure.h"

// includes from bcl - sorted alphabetically
#include "restraint/bcl_restraint_data_set_pairwise.h"
#include "score/bcl_score_data_set_pairwise_structural_exposure.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateDataSetPairwiseFilterExposure::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateDataSetPairwiseFilterExposure())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &MutateDataSetPairwiseFilterExposure::GetDefaultScheme()
    {
      static const std::string s_scheme( "filter_exposure");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor taking optional scheme
    //! @param SCHEME the scheme for this scoring function
    MutateDataSetPairwiseFilterExposure::MutateDataSetPairwiseFilterExposure( const std::string &SCHEME) :
      m_ExposureCutoff(),
      m_ExposureMap(),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking members
    //! @param EXPOSURE_CUTOFF data involving residues outside of this exposure amount will be penalized
    //! @param ENSEMBLE ensemble for which exposures will be calculated
    //! @param SCHEME the scheme for this scoring function
    MutateDataSetPairwiseFilterExposure::MutateDataSetPairwiseFilterExposure
    (
      const double &EXPOSURE_CUTOFF,
      const assemble::ProteinEnsemble &ENSEMBLE,
      const std::string &SCHEME
    ) :
      m_ExposureCutoff( EXPOSURE_CUTOFF),
      m_ExposureMap(),
      m_Scheme( SCHEME)
    {
      score::DataSetPairwiseStructuralExposure::FillExposureMap( ENSEMBLE, m_ExposureMap);
    }

    //! @brief Clone function
    //! @return pointer to new MutateDataSetPairwiseFilterExposure
    MutateDataSetPairwiseFilterExposure *MutateDataSetPairwiseFilterExposure::Clone() const
    {
      return new MutateDataSetPairwiseFilterExposure( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateDataSetPairwiseFilterExposure::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &MutateDataSetPairwiseFilterExposure::GetScheme() const
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
    MutateDataSetPairwiseFilterExposure::operator()( const DataSetPairwise &DATA_SET) const
    {
      static util::ShPtr< DataSetPairwise> s_empty;

      if( DATA_SET.IsEmpty())
      {
        BCL_MessageStd( "MutateDataSetPairwiseFilterExposure data set is empty");
        return math::MutateResult< DataSetPairwise>( s_empty, *this);
      }

      // the dataset filtered
      util::ShPtr< DataSetPairwise> filtered_data( new DataSetPairwise());

      // iterate through the data set to filter out data pairs that don't meet exposure criteria
      for
      (
        DataSetPairwise::const_iterator data_itr( DATA_SET.Begin()), data_itr_end( DATA_SET.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // score of first data point
        const double score_first
        (
          score::DataSetPairwiseStructuralExposure::CalculateExposureScore
          (
            data_itr->First(), m_ExposureMap, m_ExposureCutoff
          )
        );

        // score of second data point
        const double score_second
        (
          score::DataSetPairwiseStructuralExposure::CalculateExposureScore
          (
            data_itr->Second(), m_ExposureMap, m_ExposureCutoff
          )
        );

        BCL_MessageDbg
        (
          "first score " + util::Format()( score_first) +
          " second score " + util::Format()( score_second) + " for " + data_itr->GetIdentification()
        );

        // true if the first and second data points are defined and score is zero meaning they meet the cutoff
        if( util::IsDefined( score_first) && !score_first && util::IsDefined( score_second) && !score_second)
        {
          // add data pair to filtered_data
          filtered_data->Insert( *data_itr);
        }
      }

      BCL_MessageDbg
      (
        "MutateDataSetPairwiseFilterExposure num unique locators " +
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
    std::istream &MutateDataSetPairwiseFilterExposure::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ExposureCutoff, ISTREAM);
      io::Serialize::Read( m_ExposureMap, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateDataSetPairwiseFilterExposure::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ExposureCutoff, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ExposureMap, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint
  
} // namespace bcl
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
#include "assemble/bcl_assemble_protein_model.h"
#include "restraint/bcl_restraint_data_set_pairwise.h"
#include "restraint/bcl_restraint_mutate_data_set_pairwise_filter_triangulation.h"
#include "score/bcl_score_data_set_pairwise_coordinate_triangulation.h"

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
    const util::SiPtr< const util::ObjectInterface> MutateDataSetPairwiseFilterTriangulation::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateDataSetPairwiseFilterTriangulation())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &MutateDataSetPairwiseFilterTriangulation::GetDefaultScheme()
    {
      static const std::string s_scheme( "triangulation");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor taking optional scheme
    //! @param SCHEME the scheme for this scoring function
    MutateDataSetPairwiseFilterTriangulation::MutateDataSetPairwiseFilterTriangulation( const std::string &SCHEME) :
      m_RadiusCutoff(),
      m_Ensemble(),
      m_SortedData(),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking members
    //! @param RADIUS_CUTOFF coordinates are considered far enough apart when above this distance
    //! @param ENSEMBLE ensemble for which exposures will be calculated
      //! @param SORTED_DATA data sorted in a desired order in which it will be filtered out
    //! @param SCHEME the scheme for this scoring function
    MutateDataSetPairwiseFilterTriangulation::MutateDataSetPairwiseFilterTriangulation
    (
      const double &RADIUS_CUTOFF,
      const util::ShPtr< assemble::ProteinEnsemble> &ENSEMBLE,
      const storage::List< DataPairwise> &SORTED_DATA,
      const std::string &SCHEME
    ) :
      m_RadiusCutoff( RADIUS_CUTOFF),
      m_Ensemble( ENSEMBLE),
      m_SortedData( SORTED_DATA),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateDataSetPairwiseFilterTriangulation
    MutateDataSetPairwiseFilterTriangulation *MutateDataSetPairwiseFilterTriangulation::Clone() const
    {
      return new MutateDataSetPairwiseFilterTriangulation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateDataSetPairwiseFilterTriangulation::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &MutateDataSetPairwiseFilterTriangulation::GetScheme() const
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
    MutateDataSetPairwiseFilterTriangulation::operator()( const DataSetPairwise &DATA_SET) const
    {
      static util::ShPtr< DataSetPairwise> s_empty;

      if( DATA_SET.IsEmpty())
      {
        BCL_MessageStd( "MutateDataSetPairwiseFilterTriangulation data set is empty");
        return math::MutateResult< DataSetPairwise>( s_empty, *this);
      }

      // the dataset filtered
      util::ShPtr< DataSetPairwise> filtered_data( new DataSetPairwise());

      storage::Set
      <
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, assemble::LocatorAtomCoordinatesInterface::PtrLessThan
      > used_atoms;

      // iterate through the sorted dataset
      for
      (
        storage::List< DataPairwise>::const_iterator data_itr( m_SortedData.Begin()), data_itr_end( m_SortedData.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // try to find the data pair in DATA
        DataSetPairwise::const_iterator found_data_itr( DATA_SET.Find( *data_itr));

        // true if the current data is not found in DATA
        if( found_data_itr == DATA_SET.End())
        {
          // go to next data in sorted list of possible data
          continue;
        }
        BCL_MessageDbg
        (
          "checking " + data_itr->First()->GetIdentification()   + " and " + data_itr->Second()->GetIdentification()
        );

        // check if the data pair is usable considering the atoms that have been used already
        const bool usable_a_a
        (
          score::DataSetPairwiseCoordinateTriangulation::UsableLocator
          (
            used_atoms, *m_Ensemble, found_data_itr->First(), m_RadiusCutoff
          )
        );
        const bool usable_a_b
        (
          score::DataSetPairwiseCoordinateTriangulation::UsableLocator
          (
            used_atoms, *m_Ensemble, found_data_itr->Second(), m_RadiusCutoff
          )
        );

        // true if either aa or ab are not usable
        if( !usable_a_a || !usable_a_b)
        {
          // go to next iteration
          continue;
        }

        // add data pair to filtered_data
        filtered_data->Insert( *found_data_itr);

        // iterate through dataset again
        for
        (
          storage::List< DataPairwise>::const_iterator data_itr_b( ++storage::List< DataPairwise>::const_iterator( data_itr));
          data_itr_b != data_itr_end;
          ++data_itr_b
        )
        {
          // try to find the data pair in DATA
          DataSetPairwise::const_iterator found_data_itr_b( DATA_SET.Find( *data_itr_b));

          // true if the current data is not found in DATA
          if( found_data_itr_b == DATA_SET.End())
          {
            // go to next data in sorted list of possible data
            continue;
          }

          // check if the data pair is usable considering the atoms that have been used already
          const bool usable_b_a
          (
            score::DataSetPairwiseCoordinateTriangulation::UsableLocator
            (
              used_atoms, *m_Ensemble, found_data_itr_b->First(), m_RadiusCutoff
            )
          );
          const bool usable_b_b
          (
            score::DataSetPairwiseCoordinateTriangulation::UsableLocator
            (
              used_atoms, *m_Ensemble, found_data_itr_b->Second(), m_RadiusCutoff
            )
          );

          // true if either ba or bb are not usable
          if( !usable_b_a || !usable_b_b)
          {
            // go to next iteration
            continue;
          }

          // add data pair to filtered_data
          filtered_data->Insert( *found_data_itr_b);
        }
        BCL_MessageDbg
        (
          "finished checking " + data_itr->First()->GetIdentification() +
          " and " + data_itr->Second()->GetIdentification()
        );
      }

      BCL_MessageDbg
      (
        "MutateDataSetPairwiseFilterTriangulation data set size " +
        util::Format()( filtered_data->GetSize())
      );

      return math::MutateResult< DataSetPairwise>( filtered_data, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateDataSetPairwiseFilterTriangulation::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_RadiusCutoff, ISTREAM);
      io::Serialize::Read( m_Ensemble, ISTREAM);
      io::Serialize::Read( m_SortedData, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateDataSetPairwiseFilterTriangulation::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_RadiusCutoff, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Ensemble, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SortedData, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint
  
} // namespace bcl
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
#include "restraint/bcl_restraint_mutate_data_set_pairwise_remove.h"

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
    const util::SiPtr< const util::ObjectInterface> MutateDataSetPairwiseRemove::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateDataSetPairwiseRemove())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &MutateDataSetPairwiseRemove::GetDefaultScheme()
    {
      static const std::string s_scheme( "remove");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from optional scheme and defaults number of elements to remove to 1
    //! @param SCHEME optional scheme
    MutateDataSetPairwiseRemove::MutateDataSetPairwiseRemove( const std::string &SCHEME) :
      m_SizeRange( 1, 1),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking members
    //! @param MIN the minimum possible number of elements desired to be removed
    //! @param MAX the maximum possible number of elements desired to be removed
    //! @param SCHEME scheme for mutate
    MutateDataSetPairwiseRemove::MutateDataSetPairwiseRemove( const size_t MIN, const size_t MAX, const std::string &SCHEME) :
      m_SizeRange( MIN, MAX),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateDataSetPairwiseRemove
    MutateDataSetPairwiseRemove *MutateDataSetPairwiseRemove::Clone() const
    {
      return new MutateDataSetPairwiseRemove( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateDataSetPairwiseRemove::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this mutate
    //! @return the scheme for this mutate
    const std::string &MutateDataSetPairwiseRemove::GetScheme() const
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
    //! @param DATA DataSetPairwise of interest that will be mutated
    //! @return MutateResult that results from mutating to the argument DATA_SET
    math::MutateResult< DataSetPairwise>
    MutateDataSetPairwiseRemove::operator()( const DataSetPairwise &DATA) const
    {
      // static empty model
      static util::ShPtr< DataSetPairwise> s_empty_data_set;

      if( DATA.IsEmpty())
      {
        // return skipped move
        return math::MutateResult< DataSetPairwise>( s_empty_data_set, *this);
      }

      util::ShPtr< DataSetPairwise> new_data_set( DATA.Clone());

      // get random number of elements to remove within range
      const size_t num_to_remove( random::GetGlobalRandom().Random( m_SizeRange));

      // remove as many elements as specified
      for( size_t remove( 0); remove < num_to_remove && ( num_to_remove - remove) <= new_data_set->GetSize(); ++remove)
      {
        DataSetPairwise::const_iterator remove_element
        (
          random::GetGlobalRandom().Iterator( new_data_set->Begin(), new_data_set->End(), new_data_set->GetSize())
        );

        new_data_set->Erase( remove_element);
      }

      return math::MutateResult< DataSetPairwise>( new_data_set, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateDataSetPairwiseRemove::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SizeRange, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateDataSetPairwiseRemove::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SizeRange, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint

} // namespace bcl
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
#include "restraint/bcl_restraint_mutate_transformation_matrix_3d_null.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_transformation_matrix_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MutateTransformationMatrix3DNull::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateTransformationMatrix3DNull())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateTransformationMatrix3DNull::MutateTransformationMatrix3DNull()
    {
    }

    //! @brief Clone is the virtual copy constructor
    //! @return new instance of this class
    MutateTransformationMatrix3DNull *MutateTransformationMatrix3DNull::Clone() const
    {
      return new MutateTransformationMatrix3DNull( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateTransformationMatrix3DNull::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator which mutates a math::TransformationMatrix3D - in this case it does nothing
    //! @param MATRIX the math::TransformationMatrix3D which is to be mutated
    //! @return returns a math::MutateResult with the mutated math::TransformationMatrix3D - in this case unchanged
    math::MutateResult< math::TransformationMatrix3D> MutateTransformationMatrix3DNull::operator()
    (
      const math::TransformationMatrix3D &MATRIX
    ) const
    {
      BCL_MessageDbg( "MutateTransformationMatrix3DNull::operator()\n");

      // util::ShPtr< math::TransformationMatrix3D> sp_matrix( MATRIX.Clone()); //< this does not work

      // create ShPtr to math::TransformationMatrix3D initialized with "MATRIX"
      util::ShPtr< math::TransformationMatrix3D> sp_matrix( new math::TransformationMatrix3D( MATRIX));

      // return a math::MutateResult containing the mutated (unchanged) "sp_matrix"
      return math::MutateResult< math::TransformationMatrix3D>( sp_matrix, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateTransformationMatrix3DNull::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &MutateTransformationMatrix3DNull::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_mutate_transformation_matrix_3d_rotate.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_transformation_matrix_3d.h"
#include "random/bcl_random_uniform_distribution.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MutateTransformationMatrix3DRotate::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateTransformationMatrix3DRotate())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateTransformationMatrix3DRotate::MutateTransformationMatrix3DRotate() :
      m_Axis(),
      m_MaxRotation(),
      m_MinRotation()
    {
    }

    //! @brief construct from an axis and an amount of rotation
    //! @param AXIS is the math::RotationAxis about which the rotation will be performed
    //! @param MAX_ROTATION is a double denoting the amount of of rotation
    MutateTransformationMatrix3DRotate::MutateTransformationMatrix3DRotate( const coord::Axis &AXIS, const double MAX_ROTATION, const double MIN_ROTATION) :
      m_Axis( AXIS),
      m_MaxRotation( MAX_ROTATION),
      m_MinRotation( MIN_ROTATION)
    {
    }

    //! virtual copy constructor
    MutateTransformationMatrix3DRotate *MutateTransformationMatrix3DRotate::Clone() const
    {
      return new MutateTransformationMatrix3DRotate( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateTransformationMatrix3DRotate::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
    math::MutateResult< math::TransformationMatrix3D> MutateTransformationMatrix3DRotate::operator()
    (
      const math::TransformationMatrix3D &MATRIX
    ) const
    {
      BCL_MessageDbg( "MutateTransformationMatrix3DRotate::operator()");
      BCL_MessageDbg( "MutateTransformationMatrix3DRotate:: m_Axis() " + util::Format()( m_Axis));
      BCL_MessageDbg( "MutateTransformationMatrix3DRotate:: m_MaxRotation() " + util::Format()( m_MaxRotation));
      BCL_MessageDbg( "MutateTransformationMatrix3DRotate:: m_MinRotation() " + util::Format()( m_MinRotation));

      // get the transformation matrix of "MATRIX" at the origin
      util::ShPtr< math::TransformationMatrix3D> sp_origin( new math::TransformationMatrix3D());

      // rotate "origin" as defined by "m_Axis" and "m_MaxRotation", "m_MinRotation"
      if( m_MaxRotation == m_MinRotation)
      {
        sp_origin->operator()( m_Axis, m_MaxRotation);
      }
      else
      {
        sp_origin->operator()( m_Axis, random::GetGlobalRandom().Random< double>( m_MinRotation, m_MaxRotation));
      }

      // translate "origin" back to starting position
      sp_origin->operator()( MATRIX);

      // return unchanged "MATRIX"
      return math::MutateResult< math::TransformationMatrix3D>( sp_origin, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateTransformationMatrix3DRotate::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Axis, ISTREAM);
      io::Serialize::Read( m_MaxRotation, ISTREAM);
      io::Serialize::Read( m_MinRotation, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &MutateTransformationMatrix3DRotate::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Axis, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MaxRotation, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MinRotation, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_noe_data.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_pick_sse_random.h"
#include "find/bcl_find_pick_criteria_wrapper.h"
#include "fold/bcl_fold_default_mutates.h"
#include "fold/bcl_fold_mutate_protein_model_sse_add.h"
#include "fold/bcl_fold_mutate_tree.h"
#include "fold/bcl_fold_placement_sse_distance_restraint.h"
#include "fold/bcl_fold_score_weight_set.h"
#include "nmr/bcl_nmr_star_noe_handler.h"
#include "score/bcl_score_restraint_atom_distance.h"
#include "score/bcl_score_restraint_noe_attraction.h"
#include "score/bcl_score_restraint_noe_knowledge_based.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    // initialize scores
    fold::Score NOEData::e_ScoreNOERestraint( fold::GetScores().e_Undefined);
    fold::Score NOEData::e_ScoreNOEPenalty( fold::GetScores().e_Undefined);

    //! the scores used internally
    const util::SiPtr< const score::RestraintAtomDistance> NOEData::s_ScoreNOERestraint
    (
      util::Enumerated< score::RestraintAtomDistance>::AddInstance
      (
        new score::RestraintAtomDistance
        (
          score::RestraintNoeKnowledgeBased(),
          1.0,
          "noe_restraint"
        )
      )
    );
    const util::SiPtr< const score::RestraintAtomDistance> NOEData::s_ScoreNOEPenalty
    (
      util::Enumerated< score::RestraintAtomDistance>::AddInstance
      (
        new score::RestraintAtomDistance
        (
          score::RestraintNoeAttraction(),
          1.0,
          "noe_penalty"
        )
      )
    );

    // initialize mutate
    fold::Mutate NOEData::e_MutateAddSSENOE( fold::GetMutates().e_Undefined);

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> NOEData::s_Instance
    (
      util::Enumerated< Interface>::AddInstance( new NOEData())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    NOEData::NOEData() :
      m_Handler( GetDefaultHandler()),
      m_Restraints()
    {
    }

    //! @brief constructor taking member variable parameters
    //! @param HANDLER the handler that will be used to read and write the restraints
    NOEData::NOEData( const HandlerBase< util::ShPtrVector< AtomDistance> > &HANDLER) :
      m_Handler( HANDLER),
      m_Restraints()
    {
    }

    //! @brief Clone function
    //! @return pointer to new NOEData
    NOEData *NOEData::Clone() const
    {
      return new NOEData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &NOEData::GetAlias() const
    {
      static const std::string s_name( "NOE");
      return s_name;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &NOEData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the default file extension
    //! @return the default file extension
    const std::string &NOEData::GetDefaultExtension() const
    {
      static const std::string s_extension( ".noe_star");
      return s_extension;
    }

    //! @brief get the default restraint handler
    //! @return the default restraint handler
    const nmr::StarNOEHandler &NOEData::GetDefaultHandler()
    {
      static const nmr::StarNOEHandler s_handler( ".noe_star");
      return s_handler;
    }

    //! @brief gives reference to the specific type of data this restraint uses
    //! @return gives reference to the specific type of data this restraint uses
    const util::ShPtrVector< AtomDistance> &NOEData::GetAtomDistanceRestraints() const
    {
      return *m_Restraints;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initialize the scores and add them to Scores enumerator
    void NOEData::InitializeScores()
    {
      if( !e_ScoreNOERestraint.IsDefined())
      {
        // read in restraints
        if( !m_Restraints.IsDefined() && m_Handler.IsDefined() && m_Handler->Exists())
        {
          m_Restraints = util::CloneToShPtr( m_Handler->ReadRestraintsFromFile());
        }
        util::ShPtr< score::RestraintAtomDistance> noe_restraint( s_ScoreNOERestraint->Clone());
        noe_restraint->SetRestraints( m_Restraints);
        util::ShPtr< score::RestraintAtomDistance> noe_penalty( s_ScoreNOEPenalty->Clone());
        noe_penalty->SetRestraints( m_Restraints);
        e_ScoreNOERestraint = fold::GetScores().AddScore( noe_restraint);
        e_ScoreNOEPenalty = fold::GetScores().AddScore( noe_penalty);
      }
    }

    //! @brief sets the weights of scores in a weight set
    //! @param SCORE_WEIGHT_SET the score weight set that will be modified
    void NOEData::ModifyScoreWeightSet( fold::ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
      SCORE_WEIGHT_SET.SetWeight( e_ScoreNOERestraint, 5);
      SCORE_WEIGHT_SET.SetWeight( e_ScoreNOEPenalty, 5);
    }

    //! @brief initialize the mutates and add them to Mutates enumerator
    void NOEData::InitializeMutates()
    {
      if( !e_MutateAddSSENOE.IsDefined())
      {
        // random picker
        const assemble::PickSSERandom picker_sse_random;

        // pool picker
        const find::PickCriteriaWrapper
        <
          util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::DomainInterface
        > picker_sse_pool( picker_sse_random);

        // read in restraints
        if( !m_Restraints.IsDefined() && m_Handler.IsDefined() && m_Handler->Exists())
        {
          m_Restraints = util::CloneToShPtr( m_Handler->ReadRestraintsFromFile());
        }

        // placement
        const fold::PlacementSSEDistanceRestraint placement( m_Restraints);

        // add SSE next to SSE weighted by # of restraints
        e_MutateAddSSENOE = fold::GetMutates().AddMutate
        (
          fold::MutateProteinModelSSEAdd
          (
            picker_sse_pool,
            placement,
            fold::MutateTree::GetMutateTypeName( fold::MutateTree::e_Add) + "_sse_next_to_sse_noe"
          )
        );
      }
    }

    //! @brief modify the mutate tree used
    //! @param MUTATE_TREE MutateTree to be modified
    void NOEData::ModifyMutateTree( fold::MutateTree &MUTATE_TREE) const
    {
      // get the probability of normal SSE add
      const double original_prob
      (
        MUTATE_TREE.GetMutateProbability( fold::DefaultMutates::GetInstance().e_AddSSENextToSSE)
      );

      // set the original probability to zero
      MUTATE_TREE.SetMutateProbability( fold::DefaultMutates::GetInstance().e_AddSSENextToSSE, 0.0);

      // set the NOE probability
      MUTATE_TREE.SetMutateProbability( e_MutateAddSSENOE, original_prob);
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void NOEData::MergeAndModifyMutateTree( fold::MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< fold::MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }
 
    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer NOEData::GetSerializer() const
    {
      io::Serializer serial;
      serial.SetClassDescription
      (
        "Nuclear overhauser effect (NOE) restraints; allows for placement of SSEs to satisfy NOE data."
      );
      serial.AddInitializer
      (
        "",
        "Handler for reading NOEs",
        io::Serialization::GetAgent( &m_Handler),
        util::Implementation< HandlerBase< util::ShPtrVector< AtomDistance> > >( GetDefaultHandler()).GetString()
      );
      serial.AddInitializer
      (
        "spin label length",
        "# of bonds the spin label is away from the backbone",
        io::Serialization::GetAgent( &score::RestraintNMRDistanceInterface::GetSpinLabelLength()),
        "6"
      );
      return serial;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads restraints formatted for this restraint type from an istream
    //! @return istream restraints formatted for this restraint type were read from
    std::istream &NOEData::ReadRestraints( std::istream &ISTREAM)
    {
      // read the restraints
      m_Restraints = util::CloneToShPtr( m_Handler->ReadRestraints( ISTREAM));

      // endStage
      return ISTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &NOEData::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Handler, ISTREAM);
      io::Serialize::Read( m_Restraints, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &NOEData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Handler, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Restraints, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_piesa.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> Piesa::s_Instance( new Piesa());

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Piesa::Piesa()
    {
    }

    //! @brief returns a pointer to a new Piesa
    //! @return pointer to a new Piesa
    Piesa *Piesa::Clone() const
    {
      return new Piesa( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &Piesa::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief simulates a PIESA spectrum from the given SSE
    //! @param SSE SSE for which to simulate the PIESA spectrum
    //! @param B0 direction of the external magnetic field
    //! @return PIESA spectrum for the given SSE
    util::ShPtr< Piesa> Piesa::Create( const assemble::SSE &SSE, const linal::Vector3D &B0)
    {
      // normalize the magnetic field axis
      linal::Vector3D b0( B0);
      b0.Normalize();

      // compute the spectral values for each residue
      const size_t num_aas( SSE.GetSize());
      for( size_t i( 0); i < num_aas; ++i)
      {
        // get the coordinates of the relevant atoms for the calculation
        const biol::AABase &aa( *SSE.GetAA( i));
        const biol::Atom &n( aa.GetAtom( biol::GetAtomTypes().N));
        const biol::Atom &ca( aa.GetAtom( biol::GetAtomTypes().CA));
        const biol::Atom &h( aa.GetAtom( biol::GetAtomTypes().H));
        const linal::Vector3D n_coord( n.GetCenter());
        const linal::Vector3D ca_coord( ca.GetCenter());
        const linal::Vector3D h_coord( h.GetCenter());

        // compute the axes of the peptide plane
        const linal::Vector3D x_axis( ( h_coord - n_coord).Normalize());
        const linal::Vector3D n_c( -( ca_coord - n_coord));
        const linal::Vector3D y_axis( ( n_c - ( x_axis * n_c) * x_axis).Normalize());
        const linal::Vector3D z_axis( linal::CrossProduct( x_axis, y_axis));

        // compute the projection angles with the chemical shift tensor
        const linal::Vector3D b0_proj( ( ( b0 * x_axis) * x_axis + ( b0 * y_axis) * y_axis).Normalize());
        const double alpha( std::acos( b0_proj * x_axis));
        const double beta( std::acos( b0 * z_axis));
        const double gamma( 17.0 / 180.0 * math::g_Pi);

        // get the other values necessary for the computation
        // const double gmr_h1( 267.513 * std::pow( 10.0, 6.0));
        // const double gmr_n15( -27.116 * std::pow( 10.0, 6.0));
        const storage::VectorND< 3, double> cst( GetCSTensor( aa.GetType()));

        // compute the chemical shift of the N15
        const double cs_n15
        (
          cst( 0) * std::pow( std::sin( beta), 2) * std::pow( std::sin( alpha - gamma), 2) +
          cst( 1) * std::pow( std::cos( beta), 2) +
          cst( 2) * std::pow( std::sin( beta), 2) * std::pow( std::cos( alpha - gamma), 2)
        );

        // compute the dipolar coupling between h1 and n15
        const double dpc_p
        (
          62462.0 * ( 3 * std::pow( std::sin( beta), 2.0) * std::pow( std::cos( alpha), 2.0) - 1.0) / 2.0
        );

        // BCL_MessageTop( util::Format()( i) + "\t" + util::Format()( cs_n15) + "\t" + util::Format()( dpc_p));
      }

      return util::ShPtr< Piesa>( new Piesa);
    }

    //! @brief returns the principal values of the chemical shift tensor for the given amino acid
    //! @param AA_TYPE type of the amino acid for which to return the principal values of the chemical shift tensor
    //! @return the principal values of the chemical shift tensor for the given amino acid
    storage::VectorND< 3, double> Piesa::GetCSTensor( const biol::AAType &AA_TYPE)
    {
      // return storage::VectorND< 3, double>( 64.0, 77.0, 217.0);
      if( AA_TYPE == biol::GetAATypes().GLY)
      {
        return storage::VectorND< 3, double>( 41.0, 64.0, 210.0);
      }
      else
      {
        return storage::VectorND< 3, double>( 64.0, 77.0, 217.0);
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads members from an input stream
    //! @param ISTREAM input stream to read members from
    //! @return returns the input stream
    std::istream &Piesa::Read( std::istream &ISTREAM)
    {
      // read members

      return ISTREAM;
    }

    //! @brief writes members into an output stream
    //! @param OSTREAM output stream to write members into
    //! @INDENT number of indentations to use
    //! @return returns the output stream
    std::ostream &Piesa::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_pofr_data.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_mutate_tree.h"
#include "fold/bcl_fold_score_weight_set.h"
#include "io/bcl_io_serialization.h"
#include "restraint/bcl_restraint_sas_debye.h"
#include "restraint/bcl_restraint_sas_density_data.h"
#include "restraint/bcl_restraint_sas_pofr.h"
#include "score/bcl_score_pofr.h"
#include "score/bcl_score_restraint_pofr.h"
#include "score/bcl_score_restraint_saxs.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    // initialize score
    fold::Score PofrData::e_ScorePofrRestraint( fold::GetScores().e_Undefined);

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PofrData::s_Instance
    (
      util::Enumerated< Interface>::AddInstance( new PofrData())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from extension
    //! @param EXTENSION the extension used to identify files containing data for this restraint
    PofrData::PofrData() :
      m_Data(),
      m_Extension( GetDefaultExtension())
    {
    }

    //! @brief Clone function
    //! @return pointer to new PofrData
    PofrData *PofrData::Clone() const
    {
      return new PofrData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &PofrData::GetAlias() const
    {
      static const std::string s_name( "POFR");
      return s_name;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PofrData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the default file extension
    //! @return the default file extension
    const std::string &PofrData::GetDefaultExtension() const
    {
      static const std::string s_extension( ".pofr");
      return s_extension;
    }

    //! @brief gives reference to the specific type of data this restraint uses
    //! @return gives reference to the specific type of data this restraint uses
    const SasDensityData &PofrData::GetDensityData() const
    {
      return *m_Data;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer PofrData::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Interface derived class implements functionality for a restraint based on PofR data");
      serializer.AddInitializer
      (
        "Data",
        "experimental PofR curve",
        io::Serialization::GetAgent( &m_Data)
      );
      serializer.AddInitializer(
        "Extension",
        "the extension used to identify files containing pofr data",
        io::Serialization::GetAgent( &m_Extension)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initialize the scores and add them to Scores enumerator
    void PofrData::InitializeScores()
    {
      if( !e_ScorePofrRestraint.IsDefined())
      {
        // read the restraints from the file, if they aren't already defined
        if( !m_Data.IsDefined())
        {
          // reset the data
          m_Data = util::Implementation< SasDensityData>( new SasDensityData());
          *m_Data = m_Data->ReadRestraintsFromFile();
        }

        SasPofR density_result;

        // Stores the experimental data in the Density Interface Class
        density_result.SetExperimentalDensity( m_Data);
        e_ScorePofrRestraint = fold::GetScores().AddScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::RestraintPofr( density_result, score::PofR())
          )
        );
      }
    }

    //! @brief sets the weights of scores in a weight set
    //! @param SCORE_WEIGHT_SET the score weight set that will be modified
    void PofrData::ModifyScoreWeightSet( fold::ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
      SCORE_WEIGHT_SET.SetWeight( e_ScorePofrRestraint, 5000);
    }

    //! @brief modify the mutate tree used
    //! @param MUTATE_TREE MutateTree to be modified
    void PofrData::ModifyMutateTree( fold::MutateTree &MUTATE_TREE) const
    {
      // set SSE remove to 0
      MUTATE_TREE.SetMutateTypeProbability( fold::MutateTree::e_Remove, 0);
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void PofrData::MergeAndModifyMutateTree( fold::MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< fold::MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }
 
  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads restraints formatted for this restraint type from an istream
    //! @return istream restraints formatted for this restraint type were read from
    std::istream &PofrData::ReadRestraints( std::istream &ISTREAM)
    {
      // reset the data
      m_Data = util::Implementation< SasDensityData>( new SasDensityData( m_Extension));

      // read from the stream
      m_Data->ReadFromDataFile( ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PofrData::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Data, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PofrData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Data, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_pre_data.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_score_weight_set.h"
#include "nmr/bcl_nmr_star_noe_handler.h"
#include "score/bcl_score_restraint_atom_distance.h"
#include "score/bcl_score_restraint_energy_well.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    // initialize scores
    fold::Score PREData::e_ScorePRERestraint( fold::GetScores().e_Undefined);

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PREData::s_Instance
    (
      util::Enumerated< Interface>::AddInstance( new PREData())
    );

    const util::SiPtr< const score::RestraintAtomDistance> PREData::s_Score
    (
      util::Enumerated< score::RestraintAtomDistance>::AddInstance
      (
        new score::RestraintAtomDistance
        (
          score::RestraintEnergyWell( score::RestraintEnergyWell::e_PRE),
          1.0,
          "pre_restraint"
        )
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor taking member variable parameters
    //! @param HANDLER the handler that will be used to read and write the restraints
    PREData::PREData() :
      m_Handler( GetDefaultHandler()),
      m_Restraints()
    {
    }

    //! @brief constructor taking member variable parameters
    //! @param HANDLER the handler that will be used to read and write the restraints
    PREData::PREData( const HandlerBase< util::ShPtrVector< AtomDistance> > &HANDLER) :
      m_Handler( HANDLER),
      m_Restraints()
    {
    }

    //! @brief Clone function
    //! @return pointer to new PREData
    PREData *PREData::Clone() const
    {
      return new PREData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &PREData::GetAlias() const
    {
      static const std::string s_name( "PRE");
      return s_name;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PREData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the default file extension
    //! @return the default file extension
    const std::string &PREData::GetDefaultExtension() const
    {
      static const std::string s_extension( ".pre_star");
      return s_extension;
    }

    //! @brief get the default restraint handler
    //! @return the default restraint handler
    const nmr::StarNOEHandler &PREData::GetDefaultHandler()
    {
      static const nmr::StarNOEHandler s_handler( ".pre_star");
      return s_handler;
    }

    //! @brief gives reference to the specific type of data this restraint uses
    //! @return gives reference to the specific type of data this restraint uses
    const util::ShPtrVector< AtomDistance> &PREData::GetAtomDistanceRestraints() const
    {
      return *m_Restraints;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initialize the scores and add them to Scores enumerator
    void PREData::InitializeScores()
    {
      if( !e_ScorePRERestraint.IsDefined())
      {
        // read in restraints
        if( !m_Restraints.IsDefined() && m_Handler.IsDefined() && m_Handler->Exists())
        {
          m_Restraints = util::CloneToShPtr( m_Handler->ReadRestraintsFromFile());
        }

        util::ShPtr< score::RestraintAtomDistance> score( s_Score->Clone());
        score->SetRestraints( m_Restraints);
        e_ScorePRERestraint = fold::GetScores().AddScore( score);
      }
    }

    //! @brief sets the weights of scores in a weight set
    //! @param SCORE_WEIGHT_SET the score weight set that will be modified
    void PREData::ModifyScoreWeightSet( fold::ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
      SCORE_WEIGHT_SET.SetWeight( e_ScorePRERestraint, 5);
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void PREData::MergeAndModifyMutateTree( fold::MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< fold::MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }
 
  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer PREData::GetSerializer() const
    {
      io::Serializer serial;
      serial.SetClassDescription( "Paramagenetic relaxation enhancement based restraints");
      serial.AddInitializer
      (
        "",
        "Handler for reading PREs",
        io::Serialization::GetAgent( &m_Handler),
        util::Implementation< HandlerBase< util::ShPtrVector< AtomDistance> > >( GetDefaultHandler()).GetString()
      );
      serial.AddInitializer
      (
        "spin label length",
        "# of bonds the spin label is away from the backbone",
        io::Serialization::GetAgent( &score::RestraintNMRDistanceInterface::GetSpinLabelLength()),
        "6"
      );
      return serial;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PREData::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Handler, ISTREAM);
      io::Serialize::Read( m_Restraints, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PREData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Handler, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Restraints, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_rdc_assignment.h"

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
    const util::SiPtr< const util::ObjectInterface> RDCAssignment::s_Instance
    (
      GetObjectInstances().AddInstance( new RDCAssignment())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RDCAssignment::RDCAssignment() :
      m_Data()
    {
    }

    //! @brief Clone function
    //! @return pointer to new RDCAssignment
    RDCAssignment *RDCAssignment::Clone() const
    {
      return new RDCAssignment( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RDCAssignment::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief add new data into the list
    //! @param FIRST_COORD first coordinate
    //! @param SECOND_COORD second coordinate
    //! @param VALUE experimental RDC value
    void RDCAssignment::PushBack
    (
      const linal::Vector3D &FIRST_COORD,
      const linal::Vector3D &SECOND_COORD,
      const double VALUE
    )
    {
      m_Data.PushBack
      (
        storage::Triplet< linal::Vector3D, linal::Vector3D, double>( FIRST_COORD, SECOND_COORD, VALUE)
      );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RDCAssignment::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Data, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &RDCAssignment::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Data, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
  
} // namespace bcl
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
#include "restraint/bcl_restraint_rdc.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix3x3.h"
#include "linal/bcl_linal_matrix_operations.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "restraint/bcl_restraint_rdc_assignment.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> RDC::s_Instance
    (
      GetObjectInstances().AddInstance( new RDC())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RDC::RDC() :
      m_Data()
    {
    }

    //! @brief Clone function
    //! @return pointer to new RDC
    RDC *RDC::Clone() const
    {
      return new RDC( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RDC::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief normalizes the RDC value to an NH value
    void RDC::NormalizetoNH()
    {
      // calculate the NH normalization factor
      static const double s_nh_value
      (
        CalculateNormalization
        (
          chemistry::GetElementTypes().e_Nitrogen,
          chemistry::GetElementTypes().e_Hydrogen,
          biol::GetAtomTypes().N->GetBondLength( biol::GetAtomTypes().H)
        )
      );

      // iterate through the rdc data
      for
      (
        storage::Vector< storage::Triplet< DataPairwise, double, double> >::iterator
          itr( m_Data.Begin()), itr_end( m_Data.End()); itr != itr_end; ++itr
      )
      {
        // update the rdc value
        itr->Third() /=
          CalculateNormalization
          (
            itr->First().First()->GetAtomType()->GetElementType(),  // first element type
            itr->First().Second()->GetAtomType()->GetElementType(), // second element type
            itr->Second()                                        // internuclear distance
          ) /
          s_nh_value;
      }
    }

    //! @brief adjust the RDC signs relative to NH (i.e. value is correct, but sign is flipped)
    void RDC::AdjustSigns()
    {
      // iterate through the rdc data
      for
      (
        storage::Vector< storage::Triplet< DataPairwise, double, double> >::iterator
          itr( m_Data.Begin()), itr_end( m_Data.End()); itr != itr_end; ++itr
      )
      {
        // if the two gyromagnetic ratios are positive when multiplied (NH will be negative)
        if
        (
          itr->First().First()->GetAtomType()->GetElementType()->GetProperty
          (
            chemistry::ElementTypeData::e_GyromagneticRatio
          ) *
          itr->First().Second()->GetAtomType()->GetElementType()->GetProperty
          (
            chemistry::ElementTypeData::e_GyromagneticRatio
          ) >
          0.0
        )
        {
          // switch the sign of the rdc value
          itr->Third() = -itr->Third();
        }
      }
    }

    //! @brief calculates an RDC value from two coordinates and a tensor
    //! @param COORDINATES_A first coordinates
    //! @param COORDINATES_B second coordinates
    //! @param TENSOR tensor to be applied
    //! @return calculated RDC value
    double RDC::CalculateValue
    (
      const linal::Vector3D &COORDINATES_A,
      const linal::Vector3D &COORDINATES_B,
      const linal::Matrix3x3< double> &TENSOR
    )
    {
      // calculate the vector
      linal::Vector3D vector_value( linal::UnitVector( COORDINATES_B, COORDINATES_A));

      // store the values as a matrix
      const linal::Matrix< double> coordinates( 3, 1, vector_value.Begin());

      // multiply the transposed coordinates by the tensor and then by themselves again in order
      // to obtain an RDC value
      return ( coordinates.Transposed() * TENSOR * coordinates)( 0, 0);
    }

    //! @brief add new data into the Vector
    //! @param LOCATOR_A first atom locator
    //! @param LOCATOR_B second atom locator
    //! @param INTERNUCLEAR_DISTANCE internuclear distance
    //! @param VALUE experimental RDC value
    void RDC::PushBack
    (
      const assemble::LocatorAtomCoordinatesInterface &LOCATOR_A,
      const assemble::LocatorAtomCoordinatesInterface &LOCATOR_B,
      const double INTERNUCLEAR_DISTANCE,
      const double VALUE
    )
    {
      m_Data.PushBack
      (
        storage::Triplet< DataPairwise, double, double>
        (
          DataPairwise
          (
            util::ShPtr< assemble::LocatorAtomCoordinatesInterface>( LOCATOR_A.Clone()),
            util::ShPtr< assemble::LocatorAtomCoordinatesInterface>( LOCATOR_B.Clone())
          ),
          INTERNUCLEAR_DISTANCE,
          VALUE
        )
      );
    }

    //! @brief shuffles the Vector of RDCs
    void RDC::Shuffle()
    {
      std::random_shuffle( m_Data.Begin(), m_Data.End());
    }

    //! @brief generates the assignment from the protein model
    //! @param PROTEIN_MODEL protein model to be used to generate the assignment
    //! @return assignment containing coordinates of located atoms and experimental distance
    RDCAssignment RDC::GenerateAssignment( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // initialize RDC assignment
      RDCAssignment assignment;

      // iterate through the rdc data
      for
      (
        storage::Vector< storage::Triplet< DataPairwise, double, double> >::const_iterator
          itr( m_Data.Begin()), itr_end( m_Data.End()); itr != itr_end; ++itr
      )
      {
        // get the atom coordinates
        const linal::Vector3D first_atom_coords( itr->First().First()->Locate( PROTEIN_MODEL));
        const linal::Vector3D second_atom_coords( itr->First().Second()->Locate( PROTEIN_MODEL));
        const double value( itr->Third());

        // if all the values are defined and the experimental value is non-zero
        if( first_atom_coords.IsDefined() && second_atom_coords.IsDefined() && util::IsDefined( value) && value != 0.0)
        {
          // add the information to the assignment
          assignment.PushBack( first_atom_coords, second_atom_coords, value);
        }
      }

      // end
      return assignment;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RDC::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Data, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &RDC::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Data, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculate the normalization factor for the given element types and bond length
    //! @param ELEMENT_TYPE_A first element type
    //! @param ELEMENT_TYPE_B first element type
    //! @param BOND_LENGTH bond length
    //! @return the normalization factor
    double RDC::CalculateNormalization
    (
      const chemistry::ElementType &ELEMENT_TYPE_A,
      const chemistry::ElementType &ELEMENT_TYPE_B,
      const double BOND_LENGTH
    )
    {
      // calculate the normalization factor
      return ELEMENT_TYPE_A->GetProperty( chemistry::ElementTypeData::e_GyromagneticRatio) *
        ELEMENT_TYPE_B->GetProperty( chemistry::ElementTypeData::e_GyromagneticRatio) /
        math::Pow( BOND_LENGTH, 3.0);
    }

  } // namespace restraint

} // namespace bcl
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
#include "restraint/bcl_restraint_rdc_data.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_score_weight_set.h"
#include "nmr/bcl_nmr_residual_dipolar_coupling_least_square_deviation.h"
#include "nmr/bcl_nmr_star_rdc_handler.h"
#include "score/bcl_score_residual_dipolar_coupling_q_value.h"
#include "score/bcl_score_restraint_residual_dipolar_coupling.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    // initialize score
    fold::Score RDCData::e_ScoreRDCRestraint( fold::GetScores().e_Undefined);

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> RDCData::s_Instance
    (
      util::Enumerated< Interface>::AddInstance( new RDCData())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RDCData::RDCData() :
      m_Restraints(),
      m_Handler( GetDefaultHandler())
    {
    }

    //! @brief constructor taking member variable parameters
    //! @param HANDLER the handler that will be used to read and write the restraints
    RDCData::RDCData( const HandlerBase< RDC> &HANDLER) :
      m_Restraints(),
      m_Handler( HANDLER)
    {
    }

    //! @brief Clone function
    //! @return pointer to new RDCData
    RDCData *RDCData::Clone() const
    {
      return new RDCData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &RDCData::GetAlias() const
    {
      static const std::string s_name( "RDC");
      return s_name;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RDCData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the default file extension
    //! @return the default file extension
    const std::string &RDCData::GetDefaultExtension() const
    {
      static const std::string s_extension( ".rdc_star");
      return s_extension;
    }

    //! @brief get the default restraint handler
    //! @return the default restraint handler
    const nmr::StarRDCHandler &RDCData::GetDefaultHandler()
    {
      static const nmr::StarRDCHandler s_handler( ".rdc_star");
      return s_handler;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initialize the scores and add them to Scores enumerator
    void RDCData::InitializeScores()
    {
      if( !e_ScoreRDCRestraint.IsDefined())
      {
        if( !m_Restraints.IsDefined() && m_Handler.IsDefined() && m_Handler->Exists())
        {
          m_Restraints = util::CloneToShPtr( m_Handler->ReadRestraintsFromFile());
        }
        else
        {
          *m_Restraints = m_Handler->ReadRestraintsFromFile();
        }

        e_ScoreRDCRestraint = fold::GetScores().AddScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::RestraintResidualDipolarCoupling
            (
              m_Restraints,
              nmr::ResidualDipolarCouplingLeastSquareDeviation(),
              score::ResidualDipolarCouplingQValue()
            )
          )
        );
      }
    }

    //! @brief sets the weights of scores in a weight set
    //! @param SCORE_WEIGHT_SET the score weight set that will be modified
    void RDCData::ModifyScoreWeightSet( fold::ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
      SCORE_WEIGHT_SET.SetWeight( e_ScoreRDCRestraint, 5);
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void RDCData::MergeAndModifyMutateTree( fold::MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< fold::MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }
 
  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RDCData::GetSerializer() const
    {
      io::Serializer serial;
      serial.SetClassDescription
      (
        "Residual dipolar coupling restraints"
      );
      serial.AddInitializer
      (
        "",
        "Handler for reading RDCs",
        io::Serialization::GetAgent( &m_Handler),
        util::Implementation< HandlerBase< RDC> >( GetDefaultHandler()).GetString()
      );
      return serial;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RDCData::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Restraints, ISTREAM);
      io::Serialize::Read( m_Handler, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &RDCData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Restraints, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Handler, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_sas_analysis.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_cubic_spline_damped.h"
#include "math/bcl_math_statistics.h"
#include "restraint/bcl_restraint_sas_debye_interface.h"
#include "restraint/bcl_restraint_sas_experimental_and_calculated_density.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  /////////////////
  // data access //
  /////////////////

    //! @brief Compute the max dimension of the macromolecule from the p(r) curve
    //! @returns maximum pairwise distance in the macromolecule
    double SasAnalysis::ComputeDmax( const SasDensityData &EXPERIMENTAL_DATA)
    {
      // To get the last element of the Vector go the the end using the End function and move back 1 position
      storage::Vector< SasDistanceDensityPoint>::const_iterator data_itr( EXPERIMENTAL_DATA.End());
      data_itr--;

      return data_itr->GetRvalue();
    }

    //! @brief find maximum density point of the p(r) curve
    //! @returns the maximum value
    double SasAnalysis::FindDensitymax( const SasDensityData &EXPERIMENTAL_DATA)
    {
      double max_density( EXPERIMENTAL_DATA.Begin()->GetDensity());
      for
      (
        storage::Vector< SasDistanceDensityPoint>::const_iterator
          data_itr( EXPERIMENTAL_DATA.Begin()),
          data_itr_end( EXPERIMENTAL_DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        if( data_itr->GetDensity() > max_density)
        {
          max_density = data_itr->GetDensity();
        }
      }
      return max_density;
    }

    //! @brief find maximum density point of the p(r) curve
    //! @returns the maximum value
    double SasAnalysis::FindxDensitymax( const SasDensityData &EXPERIMENTAL_DATA)
    {
      double max_density( EXPERIMENTAL_DATA.Begin()->GetDensity());
      double max_r( EXPERIMENTAL_DATA.Begin()->GetRvalue());
      for
      (
        storage::Vector< SasDistanceDensityPoint>::const_iterator
          data_itr( EXPERIMENTAL_DATA.Begin()),
          data_itr_end( EXPERIMENTAL_DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        if( data_itr->GetDensity() > max_density)
        {
          max_density = data_itr->GetDensity();
          max_r = data_itr->GetRvalue();
        }
      }
      return max_r;
    }

    //! @brief scale the p(r) function to desired range
    //! @param EXPERIMENTAL_DATA - the data to scale
    //! @param SCALING_FACTOR - the value to scale all of the Intensity values by
    //! @returns scaled SasDensityData Object
    SasDensityData SasAnalysis::ScaleDensityData
    (
      const SasDensityData &EXPERIMENTAL_DATA,
      const double &SCALING_FACTOR
    )
    {
      SasDensityData scaled_data;

      scaled_data.SetBinSize( EXPERIMENTAL_DATA.GetBinSize());
      scaled_data.SetBinNumber( EXPERIMENTAL_DATA.GetBinNumber());
      scaled_data.SetDmax   ( EXPERIMENTAL_DATA.GetDmax());
      for
      (
        storage::Vector< SasDistanceDensityPoint>::const_iterator
          data_itr( EXPERIMENTAL_DATA.Begin()),
          data_itr_end( EXPERIMENTAL_DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        double r_value( data_itr->GetRvalue());
        double density( data_itr->GetDensity() * SCALING_FACTOR);
        double error( data_itr->GetError() * SCALING_FACTOR);

        SasDistanceDensityPoint scaled_point( r_value, density, error);
        scaled_data.PushBackDensity( scaled_point);
      }

      scaled_data.SetHmax( scaled_data.ComputeHmax());
      scaled_data.SetHxmax( scaled_data.ComputeHxmax());

      return scaled_data;
    }

    //! @brief shift the experimental p(r) function to overlap calculated p(r) function by aligning max peak
    //! @param EXPERIMENTAL_DATA - the data to scale
    //! @param SHIFT - the value to shift all of the r values by
    //! @returns scaled SasDensityData Object
    SasDensityData SasAnalysis::ShiftDensity( const SasDensityData &CALCULATED_DATA, const double &SHIFT)
    {
      SasDensityData shift_data;

      shift_data.SetBinSize( CALCULATED_DATA.GetBinSize());
      shift_data.SetBinNumber( CALCULATED_DATA.GetBinNumber());
      shift_data.SetDmax   ( CALCULATED_DATA.GetDmax());

      for
      (
        storage::Vector< SasDistanceDensityPoint>::const_iterator
          data_itr( CALCULATED_DATA.Begin()),
          data_itr_end( CALCULATED_DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        double r_shift_value( data_itr->GetRvalue() + SHIFT);
        double density( data_itr->GetDensity());
        double error( data_itr->GetError());

        SasDistanceDensityPoint shift_point( r_shift_value, density, error);
        shift_data.PushBackDensity( shift_point);
      }
      shift_data.SetHmax( shift_data.ComputeHmax());
      shift_data.SetHxmax( shift_data.ComputeHxmax());

      return shift_data;
    }

    //! @brief compute course integral of rectangles of width binsize
    //! @param DATA - the data to operate on
    //! @returns courseintegral value
    double SasAnalysis::ComputeCourseIntegral( const SasDensityData &DATA)
    {

      double integral( 0.0);
      double bin_size( DATA.GetBinSize());

      for
      (
        storage::Vector< SasDistanceDensityPoint>::const_iterator
          data_itr( DATA.Begin()),
          data_itr_end( DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        double density( data_itr->GetDensity());
        integral += density * bin_size;
      }

      return integral;
    }

    //! @brief function to calculate cumulative integral score for SAS pofr curves
    //! @param SAS_DATA experimental and calculated density distribution data
    //! @return integral score
    double SasAnalysis::CalculatePofRIntegralScore( const SasExperimentalAndCalculatedDensity &SAS_DATA)
    {
      SasExperimentalAndCalculatedDensity normalized_data( SAS_DATA);

      normalized_data.ScaleCalculatedDensity( 1 / SAS_DATA.GetCalculatedDensity().GetHmax());
      normalized_data.ScaleExperimentalDensity( 1 / SAS_DATA.GetExperimentalDensity().GetHmax());

      double cal_integral( ComputeCourseIntegral( normalized_data.GetCalculatedDensity()));
      double exp_integral( ComputeCourseIntegral( normalized_data.GetExperimentalDensity()));

      return exp_integral - cal_integral;
    }

    //! @brief function to calculate the excess integral of the calculated and experimental pofr curves
    //! @param SAS_DATA experimental and calculated density distribution data
    //! @return excess integral score
    double SasAnalysis::CalculatePofRExcessIntegralScore
    (
      const SasExperimentalAndCalculatedDensity &SAS_DATA
    )
    {
      double excess_integral( 0.0);

      SasExperimentalAndCalculatedDensity normalized_data( SAS_DATA);

      normalized_data.ScaleCalculatedDensity( 1 / SAS_DATA.GetCalculatedDensity().GetHmax());
      normalized_data.ScaleExperimentalDensity( 1 / SAS_DATA.GetExperimentalDensity().GetHmax());

      // Shift the calculated data
      normalized_data.ShiftDensity();

      for
      (
        storage::Vector< SasDistanceDensityPoint>::const_iterator
         cal_data_itr( normalized_data.GetCalculatedDensity().Begin()),
         exp_data_itr( normalized_data.GetExperimentalDensity().Begin()),
         cal_data_itr_end( normalized_data.GetCalculatedDensity().End());
        cal_data_itr != cal_data_itr_end;
        ++cal_data_itr, ++exp_data_itr
      )
      {
        double cal_density( cal_data_itr->GetDensity());
        double exp_density( exp_data_itr->GetDensity());

        if( cal_density > exp_density)
        {
          excess_integral += cal_density - exp_density;
        }
      }

      return excess_integral;
    }

    //! @brief function to calculate the amount of oscillation in the pofr curves
    //! @param SAS_DATA experimental and calculated density distribution data
    //! @return oscillation score
    double SasAnalysis::CalculatePofROscillationScore( const SasDensityData &DATA)
    {
      // A measure of smoothness is provided by the ratio ( ||P'|| / ||P||) / ( pi/( dmax - dmin))

      // Get the length of the distance distribution curve
      double delta( DATA.GetDmax());

      // Scale the Density Profiles
      SasDensityData scaled_data( ScaleDensityData( DATA, 1 / DATA.GetHmax()));

      // Compute the derivatives of the scaled data
      SasDensityData derivative( PofrDerivative( scaled_data));

      // Compute the Norm of the scaled pofr data
      storage::List< double> density( ConvertDensityDataToList( scaled_data));
      double density_norm( math::Statistics::Norm( density.Begin(), density.End()));

      // Compute the Norm of the derivative pofr data
      storage::List< double> der_density( ConvertDensityDataToList( derivative));
      double density_der_norm( math::Statistics::Norm( der_density.Begin(), der_density.End()));

      // Compute the oscillation score
      double oscillation_score( ( density_der_norm / density_norm) / ( math::g_Pi / delta));

      return oscillation_score;
    }

    //! @brief Compute the max momentum transfer vector value from the i(q) curve
    //! @returns maximum momentum transfer vector
    double SasAnalysis::ComputeQmax( const SasScatteringData &EXPERIMENTAL_DATA)
    {
      // To get the last element of the Vector go the the end using the End function and move back 1 position
      storage::Vector< SasScatteringPoint>::const_iterator data_itr( EXPERIMENTAL_DATA.End());
      --data_itr;

      return data_itr->GetQvalue();
    }

    //! @brief find the minimum intensity value in a SasScatteringData Object
    //! @param SAXS_DATA_OBJECT - the data to search
    //! @return the minimum intensity value
    double SasAnalysis::MinIntensity( const SasScatteringData &SAXS_DATA_OBJECT)
    {
      double min_intensity( SAXS_DATA_OBJECT.GetScatteringData().Begin()->GetIntensity());
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
          data_itr( SAXS_DATA_OBJECT.Begin()),
          data_itr_end( SAXS_DATA_OBJECT.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        if( data_itr->GetIntensity() < min_intensity)
        {
          min_intensity = data_itr->GetIntensity();
        }
      }
      return min_intensity;
    }

    //! @brief find the maximum intensity value in a SasScatteringData Object
    //! @param SAXS_DATA_OBJECT - the data to search
    //! @return the maximum intensity value
    double SasAnalysis::MaxIntensity( const SasScatteringData &SAXS_DATA_OBJECT)
    {
      double max_intensity( SAXS_DATA_OBJECT.GetScatteringData().Begin()->GetIntensity());
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
          data_itr( SAXS_DATA_OBJECT.Begin()),
          data_itr_end( SAXS_DATA_OBJECT.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        if( data_itr->GetIntensity() > max_intensity)
        {
          max_intensity = data_itr->GetIntensity();
        }
      }
      return max_intensity;
    }

    //! @brief scale the calculated intensity to align the experimental and calculated curves
    //! @param EXPERIMENTAL_DATA - the data to scale
    //! @param SCALING_FACTOR - the value to scale all of the Intensity values by
    //! @returns scaled SasScatteringData Object
    SasScatteringData SasAnalysis::ScaleData
    (
      const SasScatteringData &EXPERIMENTAL_DATA,
      const double &SCALING_FACTOR
    )
    {
      SasScatteringData scaled_data;
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
          data_itr( EXPERIMENTAL_DATA.Begin()),
          data_itr_end( EXPERIMENTAL_DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        double q_value( data_itr->GetQvalue());
        double intensity( data_itr->GetIntensity() * SCALING_FACTOR);
        double error( data_itr->GetError() * SCALING_FACTOR);

        SasScatteringPoint scaled_point( q_value, intensity, error);
        scaled_data.PushBackScattering( scaled_point);
      }
      return scaled_data;
    }

    //! @brief move the data to all positive numbers while maintaining identical morphology
    //! @brief This move will not adjust errors
    //! @param SAXS_DATA_OBJECT - the data to move
    //! @param MIN_VALUE - the value to base the move from
    //! @return data set with all positive values
    SasScatteringData SasAnalysis::SlideData( const SasScatteringData &SAXS_DATA_OBJECT, const double &MIN_VALUE)
    {
      SasScatteringData slid_data;
      double offset( math::Absolute( MIN_VALUE) + 1);

      BCL_MessageStd( "offset: " + util::Format()( offset));

      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
          data_itr( SAXS_DATA_OBJECT.Begin()),
          data_itr_end( SAXS_DATA_OBJECT.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        double q_value( data_itr->GetQvalue());
        double intensity( data_itr->GetIntensity() + offset);
        double error( data_itr->GetError());

        SasScatteringPoint scaled_point( q_value, intensity, error);
        slid_data.PushBackScattering( scaled_point);
      }
      return slid_data;
    }

    //! @brief take the data to log base 10
    SasScatteringData SasAnalysis::Log10( const SasScatteringData &EXPERIMENTAL_DATA)
    {
      SasScatteringData log_base10_data;

      // iterate over data to get values
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
        itr( EXPERIMENTAL_DATA.Begin()),
        itr_end( EXPERIMENTAL_DATA.End());
        itr != itr_end;
        ++itr
      )
      {
        // Store intensity
        double q_value( itr->GetQvalue());
        double intensity( log10( itr->GetIntensity()));
        double error( itr->GetError());

        // transform the error to a logaithmic scale.  If the error is greater than the intensity, only compute the
        // transformation based off the upper value ( should not happen),
        // otherwise take the average of the upper and lower values.

        if( itr->GetError() >= itr->GetIntensity())
        {
          error = log10( ( itr->GetIntensity() + itr->GetError()) / itr->GetIntensity());
        }
        else
        {
          error = 0.5 * log10( ( itr->GetIntensity() + itr->GetError()) / ( itr->GetIntensity() - itr->GetError()));
        }
        SasScatteringPoint log_point( q_value, intensity, error);
        log_base10_data.PushBackScattering( log_point);
      }
      return log_base10_data;
    }

    //! @brief transform the data from log base 10 scale to absolute scale
    //! @param DATA_OBJECT - the data to compute the log10 values of
    //! @returns data on absolute scale
    SasScatteringData SasAnalysis::LogtoAbsolute( const SasScatteringData &DATA_OBJECT)
    {
      SasScatteringData absolute_data;

      // iterate over data to get values
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
        itr( DATA_OBJECT.Begin()),
        itr_end( DATA_OBJECT.End());
        itr != itr_end;
        ++itr
      )
      {
        // Store intensity
        double q_value( itr->GetQvalue());
        double intensity( pow( 10, itr->GetIntensity()));
        double error( itr->GetError());

        double transform_error_1( error / 0.5);
        double transform_error_2( pow( 10, transform_error_1));

        double numerator( intensity * ( transform_error_2 - 1));
        double denominator( transform_error_2 + 1);

        double final_error( numerator / denominator);

        SasScatteringPoint absolute_point( q_value, intensity, final_error);
        absolute_data.PushBackScattering( absolute_point);
      }
      return absolute_data;
    }

    //! @brief take the pofr derivative of the data
    SasDensityData SasAnalysis::PofrDerivative( const SasDensityData &DATA)
    {

      const size_t data_size( DATA.GetDensitySize());

      // initialize math vectors with size of data set
      linal::Vector< double> data_values( data_size);

      // calculate delta
      storage::Vector< SasDistanceDensityPoint>::const_iterator delta_itr( DATA.Begin());
      storage::Vector< SasDistanceDensityPoint>::const_iterator delta_next( delta_itr);
      ++delta_next;

      const double delta( delta_next->GetRvalue() - delta_itr->GetRvalue());

      // populate math vectors with values from SAXS_DATA
      size_t pos( 0);

      for
      (
        storage::Vector< SasDistanceDensityPoint>::const_iterator
         data_itr( DATA.Begin()),
         data_itr_end( DATA.End());
        data_itr != data_itr_end;
        ++data_itr, ++pos
      )
      {
        data_values( pos) = data_itr->GetDensity();
      }

      // Create CublicSpline objects for the dataset
      math::CubicSplineDamped data_function;

      // Train the data set
      data_function.Train
      (
        delta_itr->GetRvalue(),
        delta,
        data_values
      );

      SasDensityData derivative_data;

      // Use the Cubic Spline Functions to Calculate the Derivative at provided point
      // Push the result into the Storage Containers
      for
      (
        storage::Vector< SasDistanceDensityPoint>::const_iterator
         data_itr( DATA.Begin()),
         data_itr_end( DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        double r_value( data_itr->GetRvalue());
        double density( data_function.dF( data_itr->GetRvalue()));
        double error( data_itr->GetError());

        SasDistanceDensityPoint derivative_point( r_value, density, error);
        derivative_data.PushBackDensity( derivative_point);
      }
      return derivative_data;
    }

    //! @brief take the derivative of the data using variable delta
    SasScatteringData SasAnalysis::Derivative( const SasScatteringData &EXPERIMENTAL_DATA)
    {
      const size_t data_size( EXPERIMENTAL_DATA.GetScatteringSize());

      // initialize math vectors with size of data set
      linal::Vector< double> x_values( data_size);
      linal::Vector< double> data_values( data_size);

      // populate math vectors with values from SAXS_DATA
      size_t pos( 0);

      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
          data_itr( EXPERIMENTAL_DATA.Begin()),
          data_itr_end( EXPERIMENTAL_DATA.End());
        data_itr != data_itr_end;
        ++data_itr, ++pos
      )
      {
        x_values( pos) = data_itr->GetQvalue();
        data_values( pos) = data_itr->GetIntensity();
      }

      //BCL_MessageStd( " data_values: " + util::Format()( data_values));

      // Create CublicSpline objects for the dataset
      math::CubicSplineDamped data_function;

      // Train the data set
      data_function.Train( x_values, data_values);

      SasScatteringData derivative_data;

      // Use the Cubic Spline Functions to Calculate the Derivative at provided point
      // Push the result into the Storage Containers
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
          data_itr( EXPERIMENTAL_DATA.Begin()),
          data_itr_end( EXPERIMENTAL_DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        double q_value( data_itr->GetQvalue());
        double intensity( data_function.dF( data_itr->GetQvalue()));
        double error( data_itr->GetError());

        //BCL_MessageStd( " derivative vd q_value: " + util::Format()( q_value));

        SasScatteringPoint derivative_point( q_value, intensity, error);
        derivative_data.PushBackScattering( derivative_point);
      }

      //BCL_MessageStd( " derivative vd values are: " + util::Format()( derivative_data));
      return derivative_data;
    }

    //! @brief Compute Experimental error from Simulated SAXS curves.  example: Crysol file
    //! @param SIMULATED_DATA - the data to simulate the experimental error for
    //! @return data set with simulated experimental error
    SasScatteringData SasAnalysis::AddErrors( const SasScatteringData &SIMULATED_DATA)
    {
      SasScatteringData simulated_errors;

      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
          data_itr( SIMULATED_DATA.Begin()),
          data_itr_end( SIMULATED_DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        double q_value( data_itr->GetQvalue());
        double intensity( data_itr->GetIntensity());
        double error( ComputeError( q_value, intensity));

        SasScatteringPoint simulated_error_point( q_value, intensity, error);
        simulated_errors.PushBackScattering( simulated_error_point);
      }
      return simulated_errors;
    }

    //! @brief Compute Experimental error for given Data point
    //! @param Q - momentum transfer vector
    //! @param INTENISTY - Intensity value for given momentum transfer vector
    //! @return data set with simulated experimental error
    double SasAnalysis::ComputeError( const double &Q, const double &INTENSITY)
    {
      double poisson_noise( math::Absolute( random::GetGlobalRandom().RandomPoisson( 10) / 10.0 - 1.0) + 1);
      double error( 0.15 * INTENSITY * ( Q + 0.001) * poisson_noise);
      return error;
    }

    //! @brief Set experimental error to defined value
    //! @param SIMULATED_DATA - the data to set the experimental error for
    //! @param ERROR_VALUE - the value to set all the experimental errors to
    //! @return data set with experimental errors set to the error value
    SasScatteringData SasAnalysis::SetErrors( const SasScatteringData &SIMULATED_DATA, const double &ERROR_VALUE)
    {
      SasScatteringData set_errors;
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
          data_itr( SIMULATED_DATA.Begin()),
          data_itr_end( SIMULATED_DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        double q_value( data_itr->GetQvalue());
        double intensity( data_itr->GetIntensity());
        double error( ERROR_VALUE);

        SasScatteringPoint set_error_point( q_value, intensity, error);
        set_errors.PushBackScattering( set_error_point);
      }
      return set_errors;
    }

    //! @brief function to calculate scaling factor for different SAXS intensity curves.  This version reads the error
    //! @brief associated with each q value.  The default error value is set to 1.0.
    //! @param SAXS_DATA experimental and calculated saxs data
    //! @return scaling factor
    double SasAnalysis::CalculateScalingWeight
    (
      const SasExperimentalAndCalculatedData &SAXS_DATA,
      const bool &USE_ERRORS
    )
    {
      // initialize sum
      double numerator( 0.0);
      double denominator( 0.0);

      // iterate over experimental and calculated data to get values
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
          exp_data_itr( SAXS_DATA.GetExperimentalData().Begin()),
          cal_data_itr( SAXS_DATA.GetCalculatedData().Begin()),
          exp_data_itr_end( SAXS_DATA.GetExperimentalData().End());
        exp_data_itr != exp_data_itr_end;
        ++exp_data_itr, ++cal_data_itr
      )
      {

        if( USE_ERRORS)
        {
          // if both I values are defined
          numerator   += exp_data_itr->GetIntensity() * cal_data_itr->GetIntensity() / math::Sqr( exp_data_itr->GetError());
          denominator += math::Sqr( cal_data_itr->GetIntensity()) / math::Sqr( exp_data_itr->GetError());
        }
        else
        {
          numerator   += exp_data_itr->GetIntensity() * cal_data_itr->GetIntensity();
          denominator += math::Sqr( cal_data_itr->GetIntensity());
        }
      }

      // endscaled_data.GetCalculatedDensity()
      return denominator == 0.0 ? 0.0 : numerator / denominator;
    }

    //! @brief function to calculate scaling factor for different SAXS intensity curves
    //! @param SAXS_DATA experimental and calculated saxs data
    //! @return scaling factor
    double SasAnalysis::CalculateStovgaardScalingWeight( const SasExperimentalAndCalculatedData &SAXS_DATA)
    {
      // initialize sum
      double numerator( 0.0);
      double denominator( 0.0);
      double alpha( 0.15);
      double beta( 0.30);
      double sigma( 0.0);
      double sqrSigma( 0.0);

      // iterate over experimental and calculated data to get values
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
         exp_data_itr( SAXS_DATA.GetExperimentalData().Begin()),
         cal_data_itr( SAXS_DATA.GetCalculatedData().Begin()),
         exp_data_itr_end( SAXS_DATA.GetExperimentalData().End());
        exp_data_itr != exp_data_itr_end;
        ++exp_data_itr, ++cal_data_itr
      )
      {
        // if both I values are defined
        sigma = ( exp_data_itr->GetIntensity() * ( exp_data_itr->GetQvalue() + alpha) * beta);
        sqrSigma = math::Sqr( sigma);

        numerator   += ( exp_data_itr->GetIntensity() * cal_data_itr->GetIntensity()) / sqrSigma;
        denominator += ( math::Sqr( cal_data_itr->GetIntensity())) / sqrSigma;
      }

      return denominator == 0.0 ? 0.0 : numerator / denominator;
    }

    //! @brief function to normalize the experimental and calculated curves individually
    //! @param SAXS_DATA scattering profile to normalize.  Range is between 0 and 1
    //! @return vector of normalized intensities
    SasScatteringData SasAnalysis::NormalizeData( const SasScatteringData &SAXS_DATA)
    {
      SasScatteringData normalized_data( SAXS_DATA);

      // Create list of doubles to hold raw intensity data
      storage::List< double> intensity( ConvertIntensityDataToList( SAXS_DATA));

      // Get the Norm of the experimental data to divide the errors by
      double exp_norm( math::Statistics::Norm( intensity.Begin(), intensity.End()));

      // Normalize the data set
      math::Statistics::Normalize( intensity.Begin(), intensity.End());

      // Write the Normalized data to the ScatteringData Object and scale the error by the norm value
      for
      (
        storage::Vector< SasScatteringPoint>::iterator
         norm_data_itr( normalized_data.Begin()),
         norm_data_itr_end( normalized_data.End());
         norm_data_itr != norm_data_itr_end;
        ++norm_data_itr
      )
      {
        const double normalized_value( intensity.FirstElement());
        norm_data_itr->SetIntensity( normalized_value);
        norm_data_itr->SetError( norm_data_itr->GetError() / exp_norm);
        intensity.PopFront();
      }

      //return SasScatteringData;
      return normalized_data;
    }

    //! @brief function to convert Scattering Intensities to a list of doubles
    //! @param SAXS_DATA input scattering profile to normalize.
    //! @return list of intensities of type double
    storage::List< double> SasAnalysis::ConvertIntensityDataToList( const SasScatteringData &SAXS_DATA)
    {
      // Create list of doubles to hold raw data
      storage::List< double> intensity;
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
        data_itr( SAXS_DATA.Begin()),
        data_itr_end( SAXS_DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        intensity.PushBack( data_itr->GetIntensity());
      }

      return intensity;
    }

    //! @brief function to convert Scattering Intensities to a list of doubles
    //! @param SAXS_DATA input scattering profile to normalize.
    //! @return list of intensities of type double
    storage::List< double> SasAnalysis::ConvertDensityDataToList( const SasDensityData &POFR_DATA)
    {
      // Create list of doubles to hold raw data
      storage::List< double> density;
      for
      (
        storage::Vector< SasDistanceDensityPoint>::const_iterator
        data_itr( POFR_DATA.Begin()),
        data_itr_end( POFR_DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        density.PushBack( data_itr->GetDensity());
      }

      return density;
    }

    //! @brief function to create Debye Implementation (GPU or CPU) based on parameters
    //! @param APPROXIMATE_LOOPS set to true to approximate loop regions
    //! @param APPROXIMATE_SIDE_CHAINS set to true to approximate side chains
    //! @param C1 Value of Exluded Volume adjustment parameter
    //! @param C2 Value of Hydration Shell adjustment parameter
    //! @param USE_CPU set to true to only use CPU implementation
    //! @return either CPU or GPU debye implementation based on provided parameters
    util::Implementation< SasDebyeInterface> SasAnalysis::SetDebyeImplementation
    (
      const bool &APPROXIMATE_LOOPS,
      const bool &APPROXIMATE_SIDE_CHAINS,
      const double &C1,
      const double &C2,
      const bool &USE_CPU,
      const bool &USE_SANS,
      const double &DEUTERIUM_EXCHANGE_PARAMETER
    )
    {
      // Setup Commandline Strings for either the opencl or non-opencl version of the code

      // First set up variables
      std::string opencl_parameters;
      std::string parameters;

      // Need to update opencl for SANS

      opencl_parameters =
        "OpenCLSaxsDebye(consider loops="
        + util::Format()( APPROXIMATE_LOOPS)
        + ", analytic=0, excluded volume="
        + util::Format()( C1)
        + ", hydration shell="
        + util::Format()( C2)
        + ", approximate_sidechains="
        + util::Format()( APPROXIMATE_SIDE_CHAINS)
        + " )";

      parameters =
        "SasDebye(consider loops="
        + util::Format()( APPROXIMATE_LOOPS)
        + ", analytic=0, excluded volume="
        + util::Format()( C1)
        + ", hydration shell="
        + util::Format()( C2)
        + ", approximate_sidechains="
        + util::Format()( APPROXIMATE_SIDE_CHAINS)
        + ", use_sans="
        + util::Format()( USE_SANS)
        + ", deuterium_percentage="
        + util::Format()( DEUTERIUM_EXCHANGE_PARAMETER)
        + " )";

      // Try to use OpenCL to compute the curves with the provided parameters, if that fails us the non-openCL version
      std::stringstream err_stream;
      util::Implementation< SasDebyeInterface> sas( opencl_parameters, err_stream);

      if( !sas.IsDefined())
      {
        // use the non-opencl version
        sas = parameters;
      }

      if( USE_CPU || USE_SANS)
      {
        sas = parameters;
      }

      // return the desired implementation
      return sas;
    }
  } // namespace restraint
} // namespace bcl
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

// include header of this class
#include "restraint/bcl_restraint_sas_data_parameters.h"
// includes from bcl - sorted alphabetically
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SasDataParameters::s_Instance
    (
      GetObjectInstances().AddInstance( new SasDataParameters())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SasDataParameters::SasDataParameters() :
      m_UseSansImplementation( false),
      m_Qvalue( util::GetUndefinedDouble()),
      m_Sasa( util::GetUndefinedDouble()),
      m_ExcludedVolume( util::GetUndefinedDouble()),
      m_HydrationShell( util::GetUndefinedDouble()),
      m_DeuteriumExchangeRate( util::GetUndefinedDouble())
    {
    }

    SasDataParameters::SasDataParameters( const double &Q_VALUE) :
        m_UseSansImplementation( false),
        m_Qvalue( Q_VALUE),
        m_Sasa( util::GetUndefinedDouble()),
        m_ExcludedVolume( util::GetUndefinedDouble()),
        m_HydrationShell( util::GetUndefinedDouble()),
        m_DeuteriumExchangeRate( util::GetUndefinedDouble())
    {
    }

    //! @brief constructor from given input data
    //! @param QVALUE the scattering angle from SAXS curve ( x-axis)
    //! @param INTENSITY the wavelength intensity for a given scattering angle ( y-axis)
    //! @param ERROR - the variation in measurement for a given intensity
    SasDataParameters::SasDataParameters
    (
       const bool &SANS_IMPLEMENTATION,
       const double &Q_VALUE,
       const double &SASA_VALUE,
       const double &EXCLUDED_VOLUME,
       const double &HYDRATION_SHELL,
       const double &DEUTERIUM_EXCHANGE_RATE
    ) :
       m_UseSansImplementation( SANS_IMPLEMENTATION),
       m_Qvalue( Q_VALUE),
       m_Sasa( SASA_VALUE),
       m_ExcludedVolume( EXCLUDED_VOLUME),
       m_HydrationShell( HYDRATION_SHELL),
       m_DeuteriumExchangeRate( DEUTERIUM_EXCHANGE_RATE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SaxsScatteringPoint
    SasDataParameters *SasDataParameters::Clone() const
    {
      return new SasDataParameters( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SasDataParameters::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SasDataParameters::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_UseSansImplementation, ISTREAM);
      io::Serialize::Read( m_Qvalue, ISTREAM);
      io::Serialize::Read( m_Sasa, ISTREAM);
      io::Serialize::Read( m_ExcludedVolume, ISTREAM);
      io::Serialize::Read( m_HydrationShell, ISTREAM);
      io::Serialize::Read( m_DeuteriumExchangeRate, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SasDataParameters::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_UseSansImplementation, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Qvalue, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Sasa, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ExcludedVolume, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_HydrationShell, OSTREAM, INDENT)<< '\n';
      io::Serialize::Write( m_DeuteriumExchangeRate, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_sas_debye.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_ca_cb.h"
#include "biol/bcl_biol_atom_group_types.h"
#include "biol/bcl_biol_sasa_data.h"
#include "fold/bcl_fold_add_parabolic_loops.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_sum_function.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> SasDebye::s_Instance
    (
      util::Enumerated< SasDebyeInterface>::AddInstance( new SasDebye())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Constructor
    //! @param LOOPS bool value to represent loops that are not present in the protein model
    //! @param USE_REGULA_FALSI_APPROXIMATION whether to determine the norm factor with regula falsi (true) or
    //!        pythagorean approximation (false)
    //! @param EXCLUDED_VOLUME_PARAMETER C1 tuning parameter for data fitting
    //! @param HYDRATION_SHELL_PARAMETER C2 tuning parameter for data fitting
    //! @param SIDE_CHAIN_APPROXIMATION - true to approximate side chains, false otherwise
    //! @param REDUCED_EXPERIMENTAL_DATA - small data set of experimental data based on Shannon Sampling
    SasDebye::SasDebye
    (
      const bool LOOPS,
      const bool USE_REGULA_FALSI_APPROXIMATION,
      double EXCLUDED_VOLUME_PARAMETER,
      double HYDRATION_SHELL_PARAMETER,
      const bool SIDE_CHAIN_APPROXIMATION,
      const bool USE_SANS,
      double DEUTERIUM_EXCHANGE_PARAMETER,
      const util::ShPtr< storage::Vector< SasScatteringPoint> > REDUCED_EXPERIMENTAL_DATA
    )
    :
      m_ShouldApproximateLoops( LOOPS),
      m_DetermineAnalyticNormFactor( USE_REGULA_FALSI_APPROXIMATION),
      m_ExcludedVolumeParameter( EXCLUDED_VOLUME_PARAMETER),
      m_HydrationShellParameter( HYDRATION_SHELL_PARAMETER),
      m_ShouldApproximateSideChains( SIDE_CHAIN_APPROXIMATION),
      m_UseSans( USE_SANS),
      m_DeuteriumExchangeParameter( DEUTERIUM_EXCHANGE_PARAMETER)
    {
      SasDebyeInterface::SetReducedExperimentalData( REDUCED_EXPERIMENTAL_DATA);
    }

    //! @brief Clone function
    //! @return pointer to new SasDebye function
    SasDebye *SasDebye::Clone() const
    {
      return new SasDebye( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &SasDebye::GetAlias() const
    {
      static const std::string s_Name( "SasDebye");
      return s_Name;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SasDebye::GetClassIdentifier() const
    {
      // Get BCL Standardized Class name
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Get the positions and form factor functions for all the atoms in the protein model
    //! @param MODEL the protein model of interest
    //! @return a vector containing pairs of position and form factor function
    storage::Vector
    <
      storage::Triplet< linal::Vector3D, util::ShPtr< math::FunctionInterfaceSerializable< SasDataParameters, double> >, double>
    > SasDebye::GetAtomsAndFormFactors( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      m_CompleteAllAtoms.Reset();
      m_AtomicCoordinates.Reset();
      m_AtomSasa.Reset();

      BCL_MessageDbg( " Inside Get Atoms and Form Factors Function: ");

      // initialize vector to hold atom coordinates and Form Factor F(q) function
      storage::Vector
      <
        storage::Triplet
        <
          linal::Vector3D,
          util::ShPtr< math::FunctionInterfaceSerializable< SasDataParameters, double> >,
          double
        >
      > all_atoms;

      size_t all_atoms_size( all_atoms.GetSize());
      BCL_MessageDbg( "initial all_atoms size: " + util::Format()( all_atoms_size));

      size_t data_member_size( m_AtomicCoordinates.GetSize());
      BCL_MessageDbg( "initial Data member size: " + util::Format()( data_member_size) + "\n");

      // Get Structure Factors for Loop coordinate approximation.  The form factors for all atoms in the residue will
      // be summed together
      if( m_ShouldApproximateLoops)
      {
        BCL_MessageDbg( " Inside the Loop approximation section: ");

        // storage for loop coordinates, if m_ShouldApproximateLoops was given
        const storage::Vector< storage::Pair< biol::AAType, linal::Vector3D> >
          loop_coordinates_aa_types( fold::AddParabolicLoops( m_DetermineAnalyticNormFactor).GetLoopCoordinates( PROTEIN_MODEL));
        for
        (
          storage::Vector< storage::Pair< biol::AAType, linal::Vector3D> >::const_iterator
            itr_loop_coords_aa_types( loop_coordinates_aa_types.Begin()),
            itr_loop_coords_aa_types_end( loop_coordinates_aa_types.End());
          itr_loop_coords_aa_types != itr_loop_coords_aa_types_end;
          ++itr_loop_coords_aa_types
        )
        {
          all_atoms.PushBack
          (
            storage::Triplet< linal::Vector3D, util::ShPtr< math::FunctionInterfaceSerializable< SasDataParameters, double> >, double>
            (
              itr_loop_coords_aa_types->Second(),
              itr_loop_coords_aa_types->First()->GetStructureFactor(),
              0.0
            )
          );

          // Expand the residue to its atom types
          const storage::Set< biol::AtomType> &atoms_in_residue( itr_loop_coords_aa_types->First()->GetAllowedHeavyAtomTypes());
          for
           (
             storage::Set< biol::AtomType>::const_iterator
               itr_atom_type( atoms_in_residue.Begin()),
               itr_atom_type_end( atoms_in_residue.End());
             itr_atom_type != itr_atom_type_end;
             ++itr_atom_type
           )
          {

            biol::AAType aa_type( itr_loop_coords_aa_types->First());
            biol::AtomType atom_type( *itr_atom_type);

            m_CompleteAllAtoms.PushBack( biol::GetAtomGroupTypes().GetTypeString( aa_type, atom_type));

            m_AtomicCoordinates.PushBack( itr_loop_coords_aa_types->Second());
            m_AtomSasa.PushBack( 0.0);
          }

        } // end loop estimation
      }

      all_atoms_size = all_atoms.GetSize();
      BCL_MessageDbg( "after loop approximation all_atoms size: " + util::Format()( all_atoms_size));

      data_member_size = m_AtomicCoordinates.GetSize();
      BCL_MessageDbg( "after loop approximation Data member size: " + util::Format()( data_member_size) + "\n");

      size_t atom_number( 0);
      size_t sse_count( 0);
      size_t residue_count( 0);

      // iterate over the chains
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
          chain_itr( PROTEIN_MODEL.GetChains().Begin()), chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end; ++chain_itr
      )
      {
        // get all the structured sses
        const util::SiPtrVector< const assemble::SSE> structured_sses
        (
          ( *chain_itr)->GetSSEs( storage::Set< biol::SSType>( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND))
        );

        // if structured_sses are empty
        if( structured_sses.IsEmpty() && m_ShouldApproximateLoops)
        {
          // warn user and return empty vector
          BCL_Message( util::Message::e_Standard, "No structured SSEs found in protein model");
          continue;
        }

        // iterate over all the SSEs including coil (depending on min_sse_size parameter)
        for
        (
          storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
            sse_itr( ( *chain_itr)->GetData().Begin()), sse_itr_end( ( *chain_itr)->GetData().End());
          sse_itr != sse_itr_end; ++sse_itr
        )
        {

           BCL_MessageDbg( "sse_count: " + util::Format()( sse_count));

          // iterate over the Amino Acids in a given SSE
          for
          (
            util::ShPtrVector< biol::AABase>::const_iterator aa_itr( ( *sse_itr)->GetData().Begin()),
              aa_itr_end( ( *sse_itr)->GetData().End());
            aa_itr != aa_itr_end; ++aa_itr
          )
          {
            ++residue_count;
            std::string amino_acid( ( *aa_itr)->GetType()->GetName());
            //BCL_MessageStd( "Residue: " + util::Format()( amino_acid));

            // initialize sum function beta carbon and for missing atom structure factors if specified
            util::ShPtr< math::SumFunction< SasDataParameters, double> > cb_atom_structure_factors
            (
              new math::SumFunction< SasDataParameters, double>()
            );

            // get the atom types for this residue
            // Be sure to get the Heavy Atom Types ( all atoms except hydrogen).  This set contains all of the residues
            // That a given amino acid contains.  Not the atoms actually present in the input PDB file
            const storage::Set< biol::AtomType> &complete_atom_types( ( *aa_itr)->GetType()->GetAllowedHeavyAtomTypes());

            //BCL_MessageStd( "AtomTypes: " + util::Format()( complete_atom_types));

            // get the first side chain atom
            const biol::Atom &first_side_chain_atom( ( *aa_itr)->GetFirstSidechainAtom());

            // if you do not have experimental Sasa Data then atom_sasa will always be zero.
            double atom_sasa( 0.0);
            double beta_sasa( 0.0);

            // iterate over the atoms in a given residue
            for
            (
              storage::Set< biol::AtomType>::const_iterator atom_itr( complete_atom_types.Begin()),
                atom_itr_end( complete_atom_types.End());
              atom_itr != atom_itr_end; ++atom_itr
            )
            {
              atom_number++;

              // get the atom type
              const biol::AtomType atom_type_name( ( *aa_itr)->GetAtom( *atom_itr).GetType());

              // get the atom coordinates
              const linal::Vector3D &atom_coords( ( *aa_itr)->GetAtom( *atom_itr).GetCoordinates());

              // get the SASA data
              util::ShPtr< biol::SasaData> sp_sasa_data
              (
                PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Sasa)
              );

              if( atom_coords.IsDefined())
              {
                std::string atom( ( *atom_itr)->GetName());
              }

              // Using Solvent Excluded Surface Area
              // Make sure the atom coordinates are defined
              if( sp_sasa_data.IsDefined() && atom_coords.IsDefined())
              {
                // Get the pdbid for the atom type
                const size_t atom_id( ( *aa_itr)->GetAtom( *atom_itr).GetPdbID());
                BCL_MessageDbg( "atom_id: " + util::Format()( atom_id));

                atom_sasa = sp_sasa_data->GetData()[ atom_id - 1]->GetSolventAccessibleSurface();
                BCL_MessageDbg( "atom_sasa: " + util::Format()( atom_sasa));
              }

              // if the coordinates are defined and this is not CB
              if( atom_coords.IsDefined() && *atom_itr != first_side_chain_atom.GetType())
              {

                const std::string atomtypecheck( biol::GetAtomGroupTypes().GetTypeString( ( *aa_itr)->GetType(), *atom_itr));
                //BCL_MessageStd( "After map AtomTypes: " + util::Format()( atomtypecheck));

                // pushback into vector
                all_atoms.PushBack
                (
                  storage::Triplet< linal::Vector3D, util::ShPtr< math::FunctionInterfaceSerializable< SasDataParameters, double> >, double>
                  (
                    atom_coords,
                    util::CloneToShPtr( *biol::GetAtomGroupTypes().GetType( ( *aa_itr)->GetType(), *atom_itr)),
                    atom_sasa
                  )
                );

                m_CompleteAllAtoms.PushBack
                (
                  biol::GetAtomGroupTypes().GetTypeString( ( *aa_itr)->GetType(), *atom_itr)
                );

                m_AtomicCoordinates.PushBack( atom_coords);

                ( sp_sasa_data.IsDefined()) ? m_AtomSasa.PushBack( atom_sasa / 100) : m_AtomSasa.PushBack( 0.0);

                all_atoms_size = all_atoms.GetSize();
                BCL_MessageDbg( "Line 323 all_atoms size: " + util::Format()( all_atoms_size));

                data_member_size = m_AtomicCoordinates.GetSize();
                BCL_MessageDbg( "Line 326 Data member size: " + util::Format()( data_member_size) + "\n");
              }
              else if( first_side_chain_atom.GetCoordinates().IsDefined())
              {
                // add cb atom to structure factor calculated
                if( *atom_itr == first_side_chain_atom.GetType())
                {
                  *cb_atom_structure_factors = *biol::GetAtomGroupTypes().GetType( ( *aa_itr)->GetType(), *atom_itr);

                  m_CompleteAllAtoms.PushBack
                  (
                    biol::GetAtomGroupTypes().GetTypeString( ( *aa_itr)->GetType(), *atom_itr)
                  );

                  m_AtomicCoordinates.PushBack( first_side_chain_atom.GetCoordinates());
                  linal::Vector3D tempcoords( first_side_chain_atom.GetCoordinates());

                  ( sp_sasa_data.IsDefined()) ? m_AtomSasa.PushBack( atom_sasa / 100) : m_AtomSasa.PushBack( 0.0);

                  // Store Atom Sasa at this iteration to insert into CB later
                  beta_sasa = atom_sasa;

                  all_atoms_size = all_atoms.GetSize();
                  BCL_MessageDbg( "Line 349 all_atoms size: " + util::Format()( all_atoms_size));

                  data_member_size = m_AtomicCoordinates.GetSize();
                  BCL_MessageDbg( "Line 352 Data member size: " + util::Format()( data_member_size) + "\n:");
                }

                // add missing residues to structure factor calculation for CB
                if( m_ShouldApproximateSideChains && *atom_itr != first_side_chain_atom.GetType())
                {
                  *cb_atom_structure_factors += *biol::GetAtomGroupTypes().GetType( ( *aa_itr)->GetType(), *atom_itr);

                  // This is the case where the atom coordinates are not defined and the form factors are summed on
                  // the CB position. In this situation the SASA value will be 0

                  m_CompleteAllAtoms.PushBack
                  (

                    biol::GetAtomGroupTypes().GetTypeString( ( *aa_itr)->GetType(), *atom_itr)
                  );

                  m_AtomicCoordinates.PushBack( first_side_chain_atom.GetCoordinates());
                  linal::Vector3D tempcoords( first_side_chain_atom.GetCoordinates());

                  m_AtomSasa.PushBack( 0.0);

                  all_atoms_size = all_atoms.GetSize();
                  BCL_MessageDbg( "Line 373 all_atoms size: " + util::Format()( all_atoms_size));

                  data_member_size = m_AtomicCoordinates.GetSize();
                  BCL_MessageDbg( "Line 376 Data member size: " + util::Format()( data_member_size));
                }
              }
            }

            // if cb is defined push back the cb point that will either be the CB atom alone, or the CB atom with
            // the other side chain atoms form factor contributions summed together
            if( !cb_atom_structure_factors->GetFunction().IsEmpty())
            {
              all_atoms.PushBack
              (
                storage::Triplet< linal::Vector3D, util::ShPtr< math::FunctionInterfaceSerializable< SasDataParameters, double> >, double>
                (
                  first_side_chain_atom.GetCoordinates(),
                  cb_atom_structure_factors,
                  beta_sasa
                )
              );

              all_atoms_size = all_atoms.GetSize();
              BCL_MessageDbg
              (
                "Added Cb Atom to all_atoms: " + util::Format()( all_atoms_size) +
                "\n ------------------------------------------------------------------------- \n"
              );

            } // Atom Type Iteration
          } // AA iteration
        } // SSE iteration
      } // chain iteration

        all_atoms_size = all_atoms.GetSize();
        BCL_MessageDbg( "Line 400 all_atoms size: " + util::Format()( all_atoms_size));

        data_member_size = m_AtomicCoordinates.GetSize();
        BCL_MessageDbg
        ( "Line 403 Data member size: " + util::Format()( data_member_size) + "\n");

      return all_atoms;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes a proteinModel as input and calculates an intensity using the debye formula
    //! @param PROTEIN_MODEL
    //! @return the intensity for a given q value
    SasExperimentalAndCalculatedData SasDebye::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {

//      //********************************** Block to circumvent SAXS Calculation *********************************
//
//      // Read Processed Data File and return result
//      util::ShPtr< restraint::SasExperimentalAndCalculatedData> sp_precalculated_data( new restraint::SasExperimentalAndCalculatedData);
//
//      io::IFStream read;
//      io::File::MustOpenIFStream( read, "scaled_y_axis.data");
//      sp_precalculated_data->ReadFromDataFile( read);
//      io::File::CloseClearFStream( read);
//
//      return *sp_precalculated_data;
//
//      //*************************************** End of Temp Block

      // get the experimental SAXS data
      util::ShPtr< SasScatteringData> sp_experimental_data;

      if( !this->GetReducedExperimentalData().IsDefined())
      {
        // get the experimental SAXS data
        sp_experimental_data = this->GetExperimentalData();

        // Verify you have experimental Saxs Data
        if( !sp_experimental_data.IsDefined())
        {
          // warn user and return empty data
          BCL_Message( util::Message::e_Critical, "No experimental SAXS data found, returning empty data");
          return SasExperimentalAndCalculatedData();
        }
      }
      else
      {
        SasScatteringData experimental_data;

        for
        (
          storage::Vector< SasScatteringPoint>::const_iterator
          exp_data_itr( this->GetReducedExperimentalData()->Begin()),
          exp_data_itr_end( this->GetReducedExperimentalData()->End());
          exp_data_itr != exp_data_itr_end;
          ++exp_data_itr
        )
        {
          experimental_data.PushBackScattering( *exp_data_itr);
        }

        // use the Clone to ShPtr to create a shared pointer to experimental data
        util::ShPtr< SasScatteringData> sp_reduced_data( util::CloneToShPtr( experimental_data));

        // set sp_experimental_data to the reduced data set
        sp_experimental_data = sp_reduced_data;
      }

      storage::Vector
      <
        storage::Triplet< linal::Vector3D, util::ShPtr< math::FunctionInterfaceSerializable< SasDataParameters, double> >, double>
      > all_atoms( GetAtomsAndFormFactors( PROTEIN_MODEL));

      // create object to hold calculated data
      SasScatteringData calculated_data;
      calculated_data.AllocateScatteringMemory( sp_experimental_data->GetScatteringData().GetSize());

      const size_t number_of_atoms( all_atoms.GetSize());
      linal::Matrix< double> distance_matrix( number_of_atoms, number_of_atoms);
      size_t row( 0);

      // Compute the Distance Matrix

      for
      (
        storage::Vector
        <
          storage::Triplet< linal::Vector3D, util::ShPtr< math::FunctionInterfaceSerializable< SasDataParameters, double> >, double>
        >::const_iterator atom_itr_a( all_atoms.Begin()), atom_itr_end( all_atoms.End());
        atom_itr_a != atom_itr_end; ++atom_itr_a, ++row
      )
      {
        size_t col( row + 1);
        for
        (
          storage::Vector
          <
            storage::Triplet< linal::Vector3D, util::ShPtr< math::FunctionInterfaceSerializable< SasDataParameters, double> >, double>
          >::const_iterator atom_itr_b( atom_itr_a + 1);
          atom_itr_b != atom_itr_end; ++atom_itr_b, ++col
        )
        {
          distance_matrix( row, col) = linal::Distance( atom_itr_a->First(), atom_itr_b->First());
        }
      }

      // iterate over experimental data to get q-values
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
          data_itr( sp_experimental_data->GetScatteringData().Begin()),
          data_itr_end( sp_experimental_data->GetScatteringData().End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // variable to hold q-value
        const double &q( data_itr->GetQvalue());

        // setup variable to hold intensity
        double intensity( 0.0);
        storage::Vector< double> residues;
        residues.AllocateMemory( number_of_atoms);

        for
        (
          storage::Vector
          <
            storage::Triplet< linal::Vector3D, util::ShPtr< math::FunctionInterfaceSerializable< SasDataParameters, double> >, double>
          >::const_iterator atom_itr_a( all_atoms.Begin()), atom_itr_end( all_atoms.End());
          atom_itr_a != atom_itr_end; ++atom_itr_a
        )
        {
          // restraint::SasDataParameters data( q, atom_itr_a->Third(), m_ExcludedVolumeParameter, m_HydrationShellParameter);
          SasDataParameters data
          (
            m_UseSans,
            q,
            atom_itr_a->Third(),
            m_ExcludedVolumeParameter,
            m_HydrationShellParameter,
            m_DeuteriumExchangeParameter
          );

          // calculate form factor for given q value for residue type
          residues.PushBack( atom_itr_a->Second()->operator()( data));
        }

        // iterate over pairs in the vector
        for( size_t row( 0); row < number_of_atoms; ++row)
        {
          // Get residue a for the debye calculation
          const double residue_a( residues( row));

          // calculate I using cb-cb distance (rij) and form factors from the amino acid types
          // I(q) = (sum i=1 to M)(sum j=1 to M) Fi(q)*Fj(q) * sin(q*rij)/(q*rij)
          // which can be rewritten as
          // I(q) = (sum i=1 to M) Fi(q)^2 +  2 * Fi(q) * (sum j=i+1 to M) Fj(q) * sin(q*rij)/(q*rij)
          // here we will call (sum j=i+1 to M) Fj(q) * sin(q*rij)/(q*rij) inner sum
          double inner_sum( 0.0);

          // iterate over upper triangle of atom-atom distance pairs
          for( size_t col( row + 1); col < number_of_atoms; ++col)
          {

            // const double cb_distance( linal::Distance( atom_itr_a->First(), atom_itr_b->First()));

            // calculate form factor for given q value for residue type
            //const double residue_b( atom_itr_b->Second()->operator()( q));
            const double residue_b( residues( col));
            const double &current_distance = distance_matrix( row, col);

            // calculate x = q * current_distance
            const double x( q * current_distance);

            // Sum the intensity over the protein
            // There is a nuance at point q = 0.  The limit of the function sin(x)/x is 1 as x approaches 0.
            // Therefore, at q==0 we only sum by residue_b.  The second term goes to 1.
            if( q == 0.0)
            {
              inner_sum += residue_b;
            }
            else
            {
              inner_sum += residue_b * ( sin( x)) / x;

            }
          }
          intensity += 2.0 * residue_a * inner_sum + math::Sqr( residue_a);
        }

        // push back into storage::Vector< storage::VectorND< 3, double> >
        // first -> q, second -> I, third -> error
        calculated_data.PushBackScattering( SasScatteringPoint( q, intensity, 0.0));

      }
      SasExperimentalAndCalculatedData saxs_data( *sp_experimental_data, calculated_data);

      return saxs_data;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SasDebye::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "performs saxs debye calculation"
      );

      parameters.AddInitializer
      (
        "consider loops",
        "should loops be considered",
        io::Serialization::GetAgent( &m_ShouldApproximateLoops),
        "0"
      );
      parameters.AddInitializer
      (
        "analytic",
        "whether to determine the norm factor with regula falsi (1) or pythagorean approximation (0)",
        io::Serialization::GetAgent( &m_DetermineAnalyticNormFactor),
        "0"
      );

      parameters.AddInitializer
      (
        "excluded volume",
        "tuning parameter for excluded solvent volume",
        io::Serialization::GetAgent( &m_ExcludedVolumeParameter),
        "1.0"
      );

      parameters.AddInitializer
      (
        "hydration shell",
        "tuning parameter for hydration shell",
        io::Serialization::GetAgent( &m_HydrationShellParameter),
        "0.0"
      );

      parameters.AddOptionalInitializer
      (
        "approximate_sidechains",
        "sum up form factor contribution on cb position",
        io::Serialization::GetAgent( &m_ShouldApproximateSideChains)
      );

      parameters.AddOptionalInitializer
      (
        "use_sans",
        "use sans implementation of debye formula",
        io::Serialization::GetAgent( &m_UseSans)
      );

      parameters.AddInitializer
      (
        "deuterium_percentage",
        "percentage of deuterium comprising solvent",
        io::Serialization::GetAgent( &m_DeuteriumExchangeParameter),
        "0.0"
      );

      return parameters;
    }

  } // namespace restraint
} // namespace bcl
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
#include "linal/bcl_linal_vector.h"
#include "math/bcl_math_histogram.h"
#include "restraint/bcl_restraint_sas_analysis.h"
#include "restraint/bcl_restraint_sas_density_data.h"

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
    const util::SiPtr< const util::ObjectInterface> SasDensityData::s_Instance
    (
      util::Enumerated< HandlerBase< SasDensityData> >::AddInstance
      (
        new SasDensityData()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SasDensityData::SasDensityData( const std::string &EXTENSION) :
      HandlerBase< SasDensityData>( EXTENSION),
      m_DensityDistributionData(),
      m_BinSize( 0),
      m_BinNumber( 0),
      m_Dmax( 0),
      m_Hmax( 0),
      m_Hxmax( 0)
    {
    }

    //! @brief constructor from given input data
    SasDensityData::SasDensityData
    (
      const storage::Vector< SasDistanceDensityPoint> &DENSITY_DATA
    ) :
      HandlerBase< SasDensityData>( ".pofr"),
      m_DensityDistributionData( DENSITY_DATA),
      m_Dmax( 0),
      m_Hmax( 0),
      m_Hxmax( 0)
    {
      storage::Vector< SasDistanceDensityPoint>::const_iterator data_itr( DENSITY_DATA.Begin());
      data_itr++;

      m_BinSize = data_itr->GetRvalue();
      m_BinNumber = DENSITY_DATA.GetSize();
      m_Dmax = SasAnalysis::ComputeDmax( *this);
      m_Hmax = SasAnalysis::FindDensitymax( *this);
      m_Hxmax = SasAnalysis::FindxDensitymax( *this);
    }

    //! @brief constructor from a histogram
    SasDensityData::SasDensityData
    (
      const math::Histogram &DENSITY_HISTOGRAM,
      const double &DMAX
    ) :
      m_Dmax( DMAX)
    {
      linal::Vector< double> binning( DENSITY_HISTOGRAM.GetBinning());
      linal::Vector< double> counts( DENSITY_HISTOGRAM.GetHistogram());

      double h_max( 0.0);
      double hx_max( 0.0);

      // Convert Histogram to SasDensityData
      for
      (
        const double *x( binning.Begin()), *x_end( binning.End()), *y( counts.Begin()), *y_end( counts.End());
        x != x_end && y != y_end;
        ++x, ++y
      )
      {
        // 0.5 was subtracted from x to align the bin on the left boundary to match experimental data
        m_DensityDistributionData.PushBack( SasDistanceDensityPoint( *x - 0.5, *y, util::GetUndefined< double>()));
        if( *y > h_max)
        {
          h_max = *y;
          hx_max = *x - 0.5;
        }
      }

      m_BinSize = DENSITY_HISTOGRAM.GetBinSize();
      m_BinNumber = DENSITY_HISTOGRAM.GetNumberOfBins();
      m_Hmax = h_max;
      m_Hxmax = hx_max;
    }

    //! @brief Clone function
    //! @return pointer to new SasDensityData
    SasDensityData *SasDensityData::Clone() const
    {
      return new SasDensityData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SasDensityData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief pushback function to add object to Dataset vector
    //! @param DATAPOINT_OBJECT  //! @brief pushback function to add object to m_Data vectorDataPoint values R, P(R), and Error
    void SasDensityData::PushBackDensity( const SasDistanceDensityPoint &DENSITY_POINT_OBJECT)
    {
      m_DensityDistributionData.PushBack( DENSITY_POINT_OBJECT);
    }

    //! @param VALUE to set binsize to
    void SasDensityData::SetBinSize( const double &BIN_SIZE)
    {
      m_BinSize = BIN_SIZE;
    }

    //! @param VALUE to set binsize to
    void SasDensityData::SetBinNumber( const size_t &BIN_NUMBER)
    {
      m_BinNumber = BIN_NUMBER;
    }

    //! @param VALUE to set binsize to
    void SasDensityData::SetDmax( const double &DMAX)
    {
      m_Dmax = DMAX;
    }

    //! @param VALUE to set Hmax to
    void SasDensityData::SetHmax( const double &HMAX)
    {
      m_Hmax = HMAX;
    }

    //! @param VALUE to set Hmax to
    void SasDensityData::SetHxmax( const double &HXMAX)
    {
      m_Hxmax = HXMAX;
    }

    const double SasDensityData::ComputeHmax() const
    {
      return SasAnalysis::FindDensitymax( *this);
    }

    const double SasDensityData::ComputeHxmax() const
    {
      return SasAnalysis::FindxDensitymax( *this);
    }

  /////////////////
  // operations  //
  /////////////////

    //! @brief preallocate Density Data memory
    //! @param SIZE size to preallocate
    void SasDensityData::AllocateDensityMemory( const size_t &SIZE)
    {
      m_DensityDistributionData.AllocateMemory( SIZE);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read experimental data from BCLModel
    //! @param ISTREAM input data stream
    void SasDensityData::ReadBCLModel( std::istream &ISTREAM)
    {
      // String to hold line data
      std::string read_line;

      // skip the header line
      std::getline( ISTREAM, read_line);

      // iterate over each line in the file
      while( std::getline( ISTREAM, read_line).good())
      {
        // split the string
        const storage::Vector< std::string> split_line( util::SplitString( read_line, " "));

        double p( util::ConvertStringToNumericalValue< double>( split_line( 0)));
        double density( util::ConvertStringToNumericalValue< double>( split_line( 1)));
        double error( util::ConvertStringToNumericalValue< double>( split_line( 2)));

        // push back the data p, p(r), Error(p)
        m_DensityDistributionData.PushBack( SasDistanceDensityPoint( p, density, error));
      } // close while loop
    } // close function

    //! @brief read experimental data from Gnom
    //! @param ISTREAM input data stream
    void SasDensityData::ReadGnomData( std::istream &ISTREAM)
    {

      std::string read_line;

      // read first line of input file
      std::getline( ISTREAM, read_line);

      // split the string based on white space as a delimeter
      storage::Vector< std::string> Line( util::SplitString( read_line, " "));

      // get the size of the first line
      size_t line_size( Line.GetSize());

      // paramaters to read through gnom filetype
      bool readDDF( false);

      // while the end of the file is not reached
      while( std::getline( ISTREAM, read_line).good())
      {
        // split the string
        const storage::Vector< std::string> split_line( util::SplitString( read_line, " "));

        line_size = split_line.GetSize();
        if( line_size == 0)
        {
          continue;
        }
        if( split_line( 0) == "R" && split_line( 1) == "P(R)")
        {
          readDDF = true;
        }
        else if( split_line( 0) == "Reciprocal")
        {
          readDDF = false;
        }
        else if( readDDF && line_size == 3)
        {
          // push back the data
          m_DensityDistributionData.PushBack
          (
            SasDistanceDensityPoint
            (
              util::ConvertStringToNumericalValue< double>( split_line( 0)),  // R-value
              util::ConvertStringToNumericalValue< double>( split_line( 1)),  // P(R)-value
              util::ConvertStringToNumericalValue< double>( split_line( 2))   // error
            )
          );
        }
      } // close while loop

      storage::Vector< SasDistanceDensityPoint>::const_iterator data_itr( m_DensityDistributionData.Begin());
      data_itr++;

      m_BinSize = data_itr->GetRvalue();
      m_BinNumber = m_DensityDistributionData.GetSize();
      m_Dmax = SasAnalysis::ComputeDmax( *this);
      m_Hmax = SasAnalysis::FindDensitymax( *this);
      m_Hxmax = SasAnalysis::FindxDensitymax( *this);
    } // close the function

    //! @brief reads saxs restraints from an input stream
    //! @brief input stream to read the restraints from
    //! @return the read in restraints
    SasDensityData SasDensityData::ReadRestraints( std::istream &ISTREAM) const
    {
      SasDensityData data;
      data.ReadFromDataFile( ISTREAM);
      return data;
    }

    //! @brief reads in the member data from a formatted file containing 3 columns:
    //! @brief scattering angle q ( 4*Pi*sin(theta)/lambda) where lambda is measured in Angstroms, I is the intensity
    //! @brief at a given a value, E is the experimental error.
    //! @brief Crysol generated files must used the paramater /dro 0.0.  The algorithm does not support adding the
    //! @brief hydration layer around the molecule.
    //! @param ISTREAM input stream
    //! @param FORMAT the file format to use for reading
    //! @return istream which was read from
    std::istream &SasDensityData::ReadFromDataFile( std::istream &ISTREAM)
    {
      // build a string to hold the line information
      std::string read_line;

      enum Types { gnom, unknown};
      Types filetype( unknown);

      // read first line of input file
      std::getline( ISTREAM, read_line);

      // split the string based on white space as a delimeter
      storage::Vector< std::string> Line( util::SplitString( read_line, " "));

      // get the size of the first line
      size_t line_size( Line.GetSize());

      size_t number_blank_lines( 0);

      // ignore empty lines at the first of the file
      while( line_size <= 1)
      {
        number_blank_lines++;
        std::getline( ISTREAM, read_line);
        Line = util::SplitString( read_line, " ");
        line_size = Line.GetSize();

        if( number_blank_lines == 1000)
        {
          BCL_Exit( "Can not read file passed to SasDensityData::ReadFromDataFile, greater than 1000 blank lines", -1);
        }

      }

      // Verify the file is a GNOM file
      if( line_size >= 5 && Line( 1) == "G" && Line( 2) == "N" && Line( 3) == "O" && Line( 4) == "M")
      {
        BCL_MessageStd( "Filetype is gnom:");
        filetype = gnom;
      }
      else
      {
        BCL_Exit( "Could not read file; last read line was: " + read_line, -1);
      }

      switch( filetype)
      {
        case gnom:
          ReadGnomData( ISTREAM);
          break;
        case unknown:
        default:
          break;
      }

      size_t size( 0);

      // make sure there are r values for analysis
      size = m_DensityDistributionData.GetSize();

      BCL_Assert( size != 0, "The number of Q values are: " + util::Format()( size));

      // return the stream
      return ISTREAM;
    }

    //! @brief writes out the member data from a formatted file containing 3 columns
    //! @param OSTREAM output stream
    //! @return ostream which was written to
    std::ostream &SasDensityData::WriteToDataFile( std::ostream &OSTREAM) const
    {
      // initialize format
      const util::Format format;

      OSTREAM << format( " BCL DENSITY INPUT PARAMETERS: RValue Density Error") << '\n';

      // iterate over m_DensityDistributionData
      for
      (
        storage::Vector< SasDistanceDensityPoint>::const_iterator density_itr
        (
          m_DensityDistributionData.Begin()),
         density_itr_end( m_DensityDistributionData.End()
        );
        density_itr != density_itr_end; ++density_itr
      )
      {
        // write the data
        OSTREAM << format( density_itr->GetRvalue()) << ' '
               << format( density_itr->GetDensity()) << ' '
               << format( density_itr->GetError()) << '\n';
      }

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SasDensityData::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_DensityDistributionData, ISTREAM);
      io::Serialize::Read( m_BinSize, ISTREAM);
      io::Serialize::Read( m_BinNumber, ISTREAM);
      io::Serialize::Read( m_Dmax, ISTREAM);
      io::Serialize::Read( m_Hmax, ISTREAM);
      io::Serialize::Read( m_Hxmax, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SasDensityData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_DensityDistributionData, OSTREAM, INDENT);
      io::Serialize::Write( m_BinSize, OSTREAM, INDENT);
      io::Serialize::Write( m_BinNumber, OSTREAM, INDENT);
      io::Serialize::Write( m_Dmax, OSTREAM, INDENT);
      io::Serialize::Write( m_Hmax, OSTREAM, INDENT);
      io::Serialize::Write( m_Hxmax, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_sas_distance_density_point.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SasDistanceDensityPoint::s_Instance
    (
      GetObjectInstances().AddInstance( new SasDistanceDensityPoint())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SasDistanceDensityPoint::SasDistanceDensityPoint() :
      m_Rvalue( util::GetUndefinedDouble()),
      m_Density( util::GetUndefinedDouble()),
      m_Error( util::GetUndefinedDouble())
    {
    }

    //! @brief constructor from given input data
    //! @param QVALUE the scattering angle from SAXS curve ( x-axis)
    //! @param INTENSITY the wavelength intensity for a given scattering angle ( y-axis)
    //! @param ERROR - the variation in measurement for a given intensity
    SasDistanceDensityPoint::SasDistanceDensityPoint
    (
      const double RVALUE,
      const double DENSITY,
      const double MEASUREMENT_ERROR
    ) :
       m_Rvalue( RVALUE),
       m_Density( DENSITY),
       m_Error( MEASUREMENT_ERROR)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SasDistanceDensityPoint
    SasDistanceDensityPoint *SasDistanceDensityPoint::Clone() const
    {
      return new SasDistanceDensityPoint( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SasDistanceDensityPoint::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SasDistanceDensityPoint::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Rvalue, ISTREAM);
      io::Serialize::Read( m_Density, ISTREAM);
      io::Serialize::Read( m_Error, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SasDistanceDensityPoint::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Rvalue, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Density, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Error, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_sas_experimental_and_calculated_data.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "math/bcl_math.h"
#include "restraint/bcl_restraint_sas_analysis.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SasExperimentalAndCalculatedData::s_Instance
    (
      GetObjectInstances().AddInstance( new SasExperimentalAndCalculatedData())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SasExperimentalAndCalculatedData::SasExperimentalAndCalculatedData() :
      m_ExperimentalData(),
      m_CalculatedData()
    {
    }

    //! @brief constructor from given input data
    SasExperimentalAndCalculatedData::SasExperimentalAndCalculatedData
    (
      const SasScatteringData &EXPERIMENTAL,
      const SasScatteringData &CALCULATED
    ) :
      m_ExperimentalData( EXPERIMENTAL),
      m_CalculatedData( CALCULATED)
    {
      BCL_Assert
      (
        EXPERIMENTAL.GetScatteringData().GetSize() == CALCULATED.GetScatteringData().GetSize(),
        "Experimental and calculated SAXS data must have the same size"
      );
    }

    //! @brief Clone function
    //! @return pointer to new SasExperimentalAndCalculatedData
    SasExperimentalAndCalculatedData *SasExperimentalAndCalculatedData::Clone() const
    {
      return new SasExperimentalAndCalculatedData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SasExperimentalAndCalculatedData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  /////////////////
  // operations  //
  /////////////////

    //! @brief scale the calculated and experimental data
    //! @param SCALING_FACTOR - the value to scale all of the Intensity values by
    void SasExperimentalAndCalculatedData::ScaleData( const double &SCALING_FACTOR)
    {
      // Scale Experimental and Calculated Data
      m_ExperimentalData = SasAnalysis::ScaleData( m_ExperimentalData, SCALING_FACTOR);
      m_CalculatedData   = SasAnalysis::ScaleData( m_CalculatedData, SCALING_FACTOR);
    }

    //! @brief scale the calculated intensity to align the experimental and calculated curves
    //! @param SCALING_FACTOR - the value to scale all of the Intensity values by
    void SasExperimentalAndCalculatedData::ScaleCalculatedData( const double &SCALING_FACTOR)
    {
      m_CalculatedData = SasAnalysis::ScaleData( m_CalculatedData, SCALING_FACTOR);
    }

    //! @brief Set the Experimental Error Values of current object to the values of passed object
    //! @param DATA_OBJECT - the object with the error values to copy
    void SasExperimentalAndCalculatedData::SetExperimentalData( const SasScatteringData &DATA_OBJECT)
    {
      m_ExperimentalData = DATA_OBJECT;
    }

    //! @brief move the data to all positive numbers while maintaining identical morphology
    void SasExperimentalAndCalculatedData::SlideData( const double &Y_MIN)
    {
      m_ExperimentalData = SasAnalysis::SlideData( m_ExperimentalData, Y_MIN);
      m_CalculatedData = SasAnalysis::SlideData( m_CalculatedData, Y_MIN);
    }

    //! @brief Set the scale for the experimental and calculated intensities to the specified boundary
    void SasExperimentalAndCalculatedData::SetYScale( const double &Y_MAX)
    {

      double max_intensity( SasAnalysis::MaxIntensity( m_ExperimentalData));

      double min_intensity
      (
        std::min
        (
          SasAnalysis::MinIntensity( m_ExperimentalData), SasAnalysis::MinIntensity( m_CalculatedData)
        )
      );

      if( min_intensity < 0)
      {
        SlideData( min_intensity);
        max_intensity = SasAnalysis::MaxIntensity( m_ExperimentalData);
      }
      const double scale_factor( Y_MAX / max_intensity);
      ScaleData( scale_factor);
    }

    //! @brief Normalize the experimental data set by its largest intensity and then Normalize the calculated data
    //! @brief by its largest intensity value
    void SasExperimentalAndCalculatedData::NormalizeData()
    {
      // Normalize the 2 data sets
      m_ExperimentalData = SasAnalysis::NormalizeData( m_ExperimentalData);
      m_CalculatedData =   SasAnalysis::NormalizeData( m_CalculatedData);
    }

    //! @brief take the data to log base 10
    void SasExperimentalAndCalculatedData::Log10()
    {
      m_ExperimentalData = SasAnalysis::Log10( m_ExperimentalData);
      m_CalculatedData   = SasAnalysis::Log10( m_CalculatedData);
    }

    //! @brief transform the log10 data to the absolute scale
    void SasExperimentalAndCalculatedData::LogtoAbsolute()
    {
      m_ExperimentalData = SasAnalysis::LogtoAbsolute( m_ExperimentalData);
      m_CalculatedData = SasAnalysis::LogtoAbsolute( m_CalculatedData);
    }

    //! @brief take the derivative of the data
    void SasExperimentalAndCalculatedData::Derivative()
    {
      m_ExperimentalData = SasAnalysis::Derivative( m_ExperimentalData);
      m_CalculatedData =   SasAnalysis::Derivative( m_CalculatedData);
    }

    //! @brief compute the chi score between the experimental and calculated curves
    //! @return the chi score between the experimental and calculated curves
    double SasExperimentalAndCalculatedData::ComputeScoringFunction( bool USE_ERRORS) const
    {
      double sum( 0.0);
      // iterate over experimental and calculated data to get values
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
          exp_data_itr( m_ExperimentalData.Begin()),
          cal_data_itr( m_CalculatedData.Begin()),
          exp_data_itr_end( m_ExperimentalData.End());
        exp_data_itr != exp_data_itr_end;
        ++exp_data_itr, ++cal_data_itr
      )
      {
        if( USE_ERRORS)
        {
          sum += math::Sqr( ( exp_data_itr->GetIntensity() - cal_data_itr->GetIntensity()) / exp_data_itr->GetError());
        }
        else
        {
          sum += math::Sqr( exp_data_itr->GetIntensity() - cal_data_itr->GetIntensity());
        }
      }

      return math::Sqrt( sum / double( GetExperimentalData().GetScatteringData().GetSize()));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to std::ostream in four columns : Q-value, experimental intensity, error and calculated intensity
    //! @param OSTREAM outputstream to write to
    //! @return outputstream which was written to
    std::ostream &SasExperimentalAndCalculatedData::WriteToGnuplot( std::ostream &OSTREAM) const
    {
      OSTREAM << "Q_Value Experimental_Intensity Experimental_Error Computed_Intensity" << '\n';

      // iterate over experimental and calculated data to get values
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
          exp_data_itr( GetExperimentalData().Begin()),
          cal_data_itr( GetCalculatedData().Begin()),
          exp_data_itr_end( GetExperimentalData().End());
        exp_data_itr != exp_data_itr_end;
        ++exp_data_itr, ++cal_data_itr
      )
      {
        // write the data to the ofstream
        OSTREAM << exp_data_itr->GetQvalue()    << ' '
                << exp_data_itr->GetIntensity() << ' '
                << exp_data_itr->GetError()     << ' '
                << cal_data_itr->GetIntensity() << '\n';
      }
      return OSTREAM;
    }

    //! @brief write to filename in three columns : Q-value, calculated intensity, experimental error
    //! @param filename to write to
    std::ostream &SasExperimentalAndCalculatedData::WriteToGnomeFormat( std::ostream &OSTREAM) const
    {
      // iterate over experimental and calculated data to get values
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
          exp_data_itr( GetExperimentalData().Begin()),
          cal_data_itr( GetCalculatedData().Begin()),
          exp_data_itr_end( GetExperimentalData().End());
        exp_data_itr != exp_data_itr_end;
        ++exp_data_itr, ++cal_data_itr
      )
      {
        // write the data to the ofstream
        OSTREAM << cal_data_itr->GetQvalue()    << ' '
                << cal_data_itr->GetIntensity() << ' '
                << exp_data_itr->GetError()     << '\n';
      }
      return OSTREAM;
    }

    //! @brief write to filename in four columns : Q-value, experimental intensity, error and calculated intensity
    //! @param filename to write to
    void SasExperimentalAndCalculatedData::WriteToGnomeFileName( const std::string &FILENAME) const
    {
      io::OFStream write;
      io::File::MustOpenOFStream( write, FILENAME);

      // write the data
      WriteToGnomeFormat( write);

      // close the ofstream
      io::File::CloseClearFStream( write);
    }

    //! @brief write to filename in four columns : Q-value, experimental intensity, error and calculated intensity
    //! @param filename to write to
    void SasExperimentalAndCalculatedData::WriteToGnuplotFileName( const std::string &FILENAME) const
    {
      io::OFStream write;
      io::File::MustOpenOFStream( write, FILENAME);

      // write the data
      WriteToGnuplot( write);

      // close the ofstream
      io::File::CloseClearFStream( write);
    }

    //! @brief read in the member data from a pre-formatted file
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SasExperimentalAndCalculatedData::ReadFromDataFile( std::istream &ISTREAM)
    {
      // build a string to hold the line information
      std::string read_line;

      // read first line of input file
      std::getline( ISTREAM, read_line);

      // split the string based on white space as a delimeter
      storage::Vector< std::string> line( util::SplitString( read_line, " "));

      // get the size of the first line
      size_t line_size( line.GetSize());

      // ignore empty lines at the first of the file
      while( line_size < 1)
      {
        std::getline( ISTREAM, read_line);
        line = util::SplitString( read_line, " ");
        line_size = line.GetSize();
      }

      bool is_correct_file( false);

      // verify bcl file type
      if(
          line( 0) == "Q_Value" &&
          line( 1) == "Experimental_Intensity" &&
          line( 2) == "Experimental_Error" &&
          line( 3) == "Computed_Intensity"
        )
      {
        is_correct_file = true;
      }

      BCL_Assert( is_correct_file == true, "Incorrect Input File Type" + util::Format()( line));

      // iterate over each line in the file
      while( std::getline( ISTREAM, read_line).good())
      {
        // split the string
        const storage::Vector< std::string> split_line( util::SplitString( read_line, " "));

        double q( util::ConvertStringToNumericalValue< double>( split_line( 0)));
        double exp_intensity( util::ConvertStringToNumericalValue< double>( split_line( 1)));
        double error( util::ConvertStringToNumericalValue< double>( split_line( 2)));
        double cal_intensity( util::ConvertStringToNumericalValue< double>( split_line( 3)));

        m_ExperimentalData.PushBackScattering( SasScatteringPoint( q, exp_intensity, error));
        m_CalculatedData.PushBackScattering( SasScatteringPoint( q, cal_intensity, 0.0));
      }

      size_t experimental_size( m_ExperimentalData.GetScatteringSize());
      size_t calculated_size( m_CalculatedData.GetScatteringSize());

      BCL_Assert( experimental_size != 0 && experimental_size == calculated_size, " Data was incorrectly read");

      return ISTREAM;
    }

    //! @brief read fit file format from CRYSOL
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SasExperimentalAndCalculatedData::ReadFromCrysolFitFile
    ( std::istream &ISTREAM, const double &FIRST_EXPERIMENTAL_POINT)
    {

      // build a string to hold the line information
      std::string read_line;

      // read first line of input file
      std::getline( ISTREAM, read_line);

      // split the string based on white space as a delimeter
      storage::Vector< std::string> line( util::SplitString( read_line, " "));

      // get the size of the first line
      size_t line_size( line.GetSize());

      // ignore empty lines at the first of the file
      while( line_size < 1)
      {
        std::getline( ISTREAM, read_line);
        line = util::SplitString( read_line, " ");
        line_size = line.GetSize();
      }

      // Switch to control when to start reading data
      bool use_data( false);

      // iterate over each line in the file
      while( std::getline( ISTREAM, read_line).good())
      {
        // split the string
        const storage::Vector< std::string> split_line( util::SplitString( read_line, " "));

        double q( util::ConvertStringToNumericalValue< double>( split_line( 0)));
        double exp_intensity( util::ConvertStringToNumericalValue< double>( split_line( 1)));
        double error( util::GetUndefined< double>());
        double cal_intensity( util::ConvertStringToNumericalValue< double>( split_line( 2)));

        // Do not use extrapolated data from Crysol. Only use the data provided by experiment
        // Use within tolerance function to account for numerical drift between platforms
        if( math::EqualWithinAbsoluteTolerance( q, FIRST_EXPERIMENTAL_POINT, 0.001))
        {
          use_data = true;
        }

        if( use_data)
        {
          m_ExperimentalData.PushBackScattering( SasScatteringPoint( q, exp_intensity, error));
          m_CalculatedData.PushBackScattering( SasScatteringPoint( q, cal_intensity, error));
        }
      }

      size_t experimental_size( m_ExperimentalData.GetScatteringSize());
      size_t calculated_size( m_CalculatedData.GetScatteringSize());

      BCL_Assert( experimental_size != 0 && experimental_size == calculated_size, " Data was incorrectly read");
      return ISTREAM;
    }

    //! @brief read fit file format from CRYSOL
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SasExperimentalAndCalculatedData::ReadFromFoxsFile( std::istream &ISTREAM)
    {
      // build a string to hold the line information
      std::string read_line;

      // read first line of input file
      std::getline( ISTREAM, read_line);

      // split the string based on white space as a delimeter
      storage::Vector< std::string> line( util::SplitString( read_line, " "));

      // get the size of the first line
      size_t line_size( line.GetSize());

      // ignore empty lines at the first of the file
      while( line_size < 1)
      {
        std::getline( ISTREAM, read_line);
        line = util::SplitString( read_line, " ");
        line_size = line.GetSize();
      }

      // iterate over each line in the file
      while( std::getline( ISTREAM, read_line).good())
      {
        // split the string
        const storage::Vector< std::string> split_line( util::SplitString( read_line, " "));

        std::string first_value( split_line( 0));

        // Do not use header information from foxs. Only use the data provided by experiment
        if( first_value.compare( "#") == 0)
        {
          continue;
        }

        double q( util::ConvertStringToNumericalValue< double>( split_line( 0)));
        double exp_intensity( util::ConvertStringToNumericalValue< double>( split_line( 1)));
        double error( util::GetUndefined< double>());
        double cal_intensity( util::ConvertStringToNumericalValue< double>( split_line( 2)));

        m_ExperimentalData.PushBackScattering( SasScatteringPoint( q, exp_intensity, error));
        m_CalculatedData.PushBackScattering( SasScatteringPoint( q, cal_intensity, error));
      }

      size_t experimental_size( m_ExperimentalData.GetScatteringSize());
      size_t calculated_size( m_CalculatedData.GetScatteringSize());

      BCL_Assert( experimental_size != 0 && experimental_size == calculated_size, " Data was incorrectly read");

      return ISTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SasExperimentalAndCalculatedData::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ExperimentalData, ISTREAM);
      io::Serialize::Read( m_CalculatedData, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SasExperimentalAndCalculatedData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ExperimentalData, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_CalculatedData, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_sas_experimental_and_calculated_density.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "math/bcl_math.h"
#include "math/bcl_math_cubic_spline.h"
#include "math/bcl_math_cubic_spline_damped.h"
#include "restraint/bcl_restraint_sas_analysis.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SasExperimentalAndCalculatedDensity::s_Instance
    (
      GetObjectInstances().AddInstance( new SasExperimentalAndCalculatedDensity())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SasExperimentalAndCalculatedDensity::SasExperimentalAndCalculatedDensity() :
      m_ExperimentalDensity(),
      m_CalculatedDensity()
    {
    }

    //!@brief copy constructor
    SasExperimentalAndCalculatedDensity::SasExperimentalAndCalculatedDensity( const SasExperimentalAndCalculatedDensity &RHS)
    {
      m_ExperimentalDensity = RHS.GetExperimentalDensity();
      m_ExperimentalDmax  = RHS.GetExperimentalDensity().GetDmax();
      m_CalculatedDensity = RHS.GetCalculatedDensity();
      m_CalculatedDmax = RHS.GetCalculatedDensity().GetDmax();
    }

    //! @brief constructor from given input data
    SasExperimentalAndCalculatedDensity::SasExperimentalAndCalculatedDensity
    (
      const SasDensityData &EXPERIMENTAL,
      const SasDensityData &CALCULATED
    ) :
      m_ExperimentalDensity( EXPERIMENTAL),
      m_CalculatedDensity( CALCULATED)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SasExperimentalAndCalculatedDensity
    SasExperimentalAndCalculatedDensity *SasExperimentalAndCalculatedDensity::Clone() const
    {
      return new SasExperimentalAndCalculatedDensity( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SasExperimentalAndCalculatedDensity::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  /////////////////
  // operations  //
  /////////////////

    //! @brief scale the calculated and experimental data
    //! @param SCALING_FACTOR - the value to scale all of the Intensity values by
    void SasExperimentalAndCalculatedDensity::ScaleExperimentalDensity( const double &EXP_SCALING_FACTOR)
    {
      // Scale Experimental and Calculated Data
      m_ExperimentalDensity = SasAnalysis::ScaleDensityData( m_ExperimentalDensity, EXP_SCALING_FACTOR);
    }

    //! @brief scale the calculated and experimental data
    //! @param SCALING_FACTOR - the value to scale all of the Intensity values by
    void SasExperimentalAndCalculatedDensity::ScaleCalculatedDensity( const double &CAL_SCALING_FACTOR)
    {
      // Scale Experimental and Calculated Data
      m_CalculatedDensity = SasAnalysis::ScaleDensityData( m_CalculatedDensity, CAL_SCALING_FACTOR);
    }

    //! @brief align the max peak of the experimental and calculated curves
    void SasExperimentalAndCalculatedDensity::ShiftDensity()
    {
      double cal_hx_max( m_CalculatedDensity.GetHxmax());
      double exp_hx_max( m_ExperimentalDensity.GetHxmax());

      double shift( exp_hx_max - cal_hx_max);

      m_CalculatedDensity = SasAnalysis::ShiftDensity( m_CalculatedDensity, shift);

      // Compute spline for experimental Density
      const size_t data_size( m_ExperimentalDensity.GetDensitySize());

      // initialize math vectors with size of data set
      linal::Vector< double> data_values( data_size);

      const double delta( m_ExperimentalDensity.GetBinSize());

      // populate math vectors with values from SAXS_DATA
      size_t pos( 0);

      for
      (
        storage::Vector< SasDistanceDensityPoint>::const_iterator
         data_itr( m_ExperimentalDensity.Begin()),
         data_itr_end( m_ExperimentalDensity.End());
        data_itr != data_itr_end;
        ++data_itr, ++pos
      )
      {
        data_values( pos) = data_itr->GetDensity();
      }

      // Create CublicSpline objects for the dataset
      math::CubicSplineDamped data_function;

      // Train the data set
      data_function.Train
      (
        m_ExperimentalDensity.Begin()->GetRvalue(),
        delta,
        data_values
      );

      // iterate over calculated data using the splines from the experimental data
      for
      (
        storage::Vector< SasDistanceDensityPoint>::iterator
         cal_data_itr( m_CalculatedDensity.Begin()),
         exp_data_itr( m_ExperimentalDensity.Begin()),
         cal_data_itr_end( m_CalculatedDensity.End());
        cal_data_itr != cal_data_itr_end;
        ++cal_data_itr, ++exp_data_itr
      )
      {
        if( cal_data_itr->GetRvalue() < m_ExperimentalDensity.GetDmax())
        {
          exp_data_itr->SetDensity( data_function( cal_data_itr->GetRvalue()));
          //exp_data_itr->SetRvalue ( cal_data_itr->GetRvalue());
        }
        else
        {
          exp_data_itr->SetDensity( 0.0);
          //exp_data_itr->SetRvalue ( cal_data_itr->GetRvalue());
        }
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to filename in three columns : R-value, density, experimental error
    //! @param filename to write to
    std::ostream &SasExperimentalAndCalculatedDensity::WriteToOstream( std::ostream &OSTREAM) const
    {
      OSTREAM << "R_Value Experimental_Density Experimental_Error Computed_Density" << '\n';

      // iterate over experimental and calculated data to get values
      for
      (
        storage::Vector< SasDistanceDensityPoint>::const_iterator
          exp_data_itr( GetExperimentalDensity().Begin()),
          cal_data_itr( GetCalculatedDensity().Begin()),
          exp_data_itr_end( GetExperimentalDensity().End());
        exp_data_itr != exp_data_itr_end;
        ++exp_data_itr, ++cal_data_itr
      )
      {
        // write the data to the ofstream
        OSTREAM << exp_data_itr->GetRvalue()    << ' '
                << exp_data_itr->GetDensity()   << ' '
                << exp_data_itr->GetError()     << ' '
                << cal_data_itr->GetDensity()   << '\n';
      }
      return OSTREAM;
    }

    //! @brief write to filename in three columns : R-value, experimental density, and calculated density
    //! @param filename to write to
    void SasExperimentalAndCalculatedDensity::WriteToGnuplotFileName( const std::string &FILENAME) const
    {
      io::OFStream write;
      io::File::MustOpenOFStream( write, FILENAME);

      // write the data
      WriteToOstream( write);

      // close the ofstream
      io::File::CloseClearFStream( write);
    }

    //! @brief write to filename in four columns : Q-value, experimental intensity, error and calculated intensity
    //! @param filename to write to
    void SasExperimentalAndCalculatedDensity::WriteToFileName( const std::string &FILENAME) const
    {
      io::OFStream write;
      io::File::MustOpenOFStream( write, FILENAME);

      // write the data
      WriteToOstream( write);

      // close the ofstream
      io::File::CloseClearFStream( write);
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SasExperimentalAndCalculatedDensity::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ExperimentalDensity, ISTREAM);
      io::Serialize::Read( m_CalculatedDensity, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SasExperimentalAndCalculatedDensity::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ExperimentalDensity, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_CalculatedDensity, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_sas_optimization.h"

// includes from bcl - sorted alphabetically
#include "restraint/bcl_restraint_sas_analysis.h"
#include "restraint/bcl_restraint_sas_debye_interface.h"
#include "restraint/bcl_restraint_saxs_opti_result.h"
#include "util/bcl_util_implementation.h"

using bcl::util::Message;

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

    const float SasOptimization::default_C1Min      = 0.80;
    const float SasOptimization::default_C1Max      = 1.20;
    const float SasOptimization::default_C2Min      = 0.00;
    const float SasOptimization::default_C2Max      = 4.00;
    const float SasOptimization::default_C1StepSize = 0.005;
    const float SasOptimization::default_C2StepSize = 0.100;

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SasOptimization::s_Instance
    (
      GetObjectInstances().AddInstance( new SasOptimization())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SasOptimization::SasOptimization() :
      m_c1Min( default_C1Min),
      m_c1Max( default_C1Max),
      m_c2Min( default_C2Min),
      m_c2Max( default_C2Max),
      m_c1StepSize( default_C1StepSize),
      m_c2StepSize( default_C2StepSize),
      m_scoreType( score::SasType::e_chi),
      m_useErrors( false),
      m_numEvaluations( 0),
      m_searchGrid( ( size_t) 41, ( size_t) 81, -100.0),
      m_ApproximateSideChains( false),
      m_Sans( false),
      m_DeuteriumExchangeParameter( 0.0)
    {

      storage::Vector< SasTransformation::TransformationTypeEnum> default_transform_vector;
      default_transform_vector.PushBack( SasTransformation::e_None);
      SasTransformation default_transform( default_transform_vector, false, m_useErrors, 1.0);

      m_transformation = default_transform;
    }

    //! @brief Constructor from members
    //! @param C1_MIN - minimum C1 value
    //! @param C1_MAX - maximum C1 value
    //! @param C2_MIN - minimum C2 value
    //! @param C2_MAX - maximum C2 value
    //! @param C1_STEPSIZE - step size for C1 parameter
    //! @param C2_STEPSIZE - step size for C2 parameter
    SasOptimization::SasOptimization
    (
      const float &C1_MIN,
      const float &C1_MAX,
      const float &C2_MIN,
      const float &C2_MAX,
      const float &C1_STEPSIZE,
      const float &C2_STEPSIZE,
      const score::SasType::ScoreFunctionEnum &SCORE_FUNCTION,
      const bool &USE_ERRORS,
      const SasTransformation &TRANSFORMATION_TYPES,
      const bool &APPROXIMATE_SIDE_CHAINS,
      const bool &HARDWARE_TYPE,
      const bool &SAS_TYPE,
      const float &DEUTERIUM_EXCHANGE_PARAMETER
    )
    :
      m_c1Min( C1_MIN),
      m_c1Max( C1_MAX),
      m_c2Min( C2_MIN),
      m_c2Max( C2_MAX),
      m_c1StepSize( C1_STEPSIZE),
      m_c2StepSize( C2_STEPSIZE),
      m_scoreType( SCORE_FUNCTION),
      m_useErrors( USE_ERRORS),
      m_numEvaluations( 0),
      m_searchGrid
      (
        (size_t)( (( m_c2Max - m_c2Min) / m_c2StepSize) + 1),
        (size_t)( (( m_c1Max - m_c1Min) / m_c1StepSize) + 1),
        ( float) - 100.0
      ),
      m_transformation( TRANSFORMATION_TYPES),
      m_ApproximateSideChains( APPROXIMATE_SIDE_CHAINS),
      m_Cpu( HARDWARE_TYPE),
      m_Sans( SAS_TYPE),
      m_DeuteriumExchangeParameter( DEUTERIUM_EXCHANGE_PARAMETER)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SasOptimization
    SasOptimization *SasOptimization::Clone() const
    {
      return new SasOptimization( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SasOptimization::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief overloaded () operator to optimize c1 and c2 parameters
    //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
    //! @return returns storage vector of 3 values: c1, c2, and min chi score
    SaxsOptiResult SasOptimization::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      SaxsOptiResult result( GridWalk( PROTEIN_MODEL));

      BCL_MessageStd( "Transformation values:" + util::Format()( m_transformation));

      return result;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns lowest chi score from adusting c1 and c2 parameters over a grid of 3200 points
    //! @brief searches the grid by a line search using quadratic interpolation along each row of the grid
    //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
    //! @return returns storage vector of 3 values: c1, c2, and min chi score
    SaxsOptiResult SasOptimization::GridWalk( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      storage::Vector< float> chi_array;
      storage::Vector< float> c1_array;
      storage::Vector< float> c2_array;

      m_numEvaluations = 0;

      for( float c2( m_c2Min); c2 <= m_c2Max; c2 += m_c2StepSize)
      {
        SaxsOptiResult result( GoldenSection( m_c1Min, m_c1Max, c2, PROTEIN_MODEL));
        //SaxsOptiResult result( QuadraticInterpolation( m_c1Min, m_c1Max, c2, PROTEIN_MODEL));
        c1_array.PushBack( result.GetC1());
        c2_array.PushBack( result.GetC2());
        chi_array.PushBack( result.GetChi());
      }

      BCL_MessageStd( "chi_array: " + util::Format()( chi_array));

      //Get index of smallest element in list
      size_t index( GetMinimumIndex( chi_array));

      float final_chi( chi_array( index));
      float final_c1( c1_array( index));
      float final_c2( c2_array( index));

      SaxsOptiResult result( final_c1, final_c2, final_chi);

      return result;
    }

    //! @brief returns lowest chi score from adusting c1 and c2 parameters over a grid of 3200 points
    //! @brief for each column of the grid each row is tested until the function at a given point increases
    //! @brief this method is faster than a complete grid search, but still slow to find a final answer
    //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
    //! @return returns storage vector of 3 values: c1, c2, and min chi score
    SaxsOptiResult SasOptimization::SmartGridSearch( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // initialize variables for grid search
      float chi( 0.0);
      float lowest_chi( 10000.0);
      float lowest_c1( 0.0);
      float lowest_c2( 0.0);

      float previous( 10000.0);
      float current( 0.0);
      size_t count( 0);

      // instead of running 1 calculation, I want to run a grid search on the c1 and c2 parameters
      // c1 = [0.80 - 1.20] step = 0.005   80 total points
      // c2 = [0.0 -  4.0] step = 0.1      40 total points
      // 3200 or less total evaluations for the optimimal c1 and c2 combination to minimize chi

      for( float c1( m_c1Min); c1 <= m_c1Max; c1 += m_c1StepSize)
      {
        for( float c2( m_c2Min); c2 <= m_c2Max; c2 += m_c2StepSize)
        {
          chi = EvaluateFunction( c1, c2, PROTEIN_MODEL);

          if( chi < lowest_chi)
          {
            lowest_chi = chi;
            lowest_c1 = c1;
            lowest_c2 = c2;
          }

          current = chi;
          if( current < previous)
          {
            previous = current;
          }
          else
          {
            previous = 10000.0;
            current = 0.0;
            break;
          }
          //BCL_Message( util::Message::e_Standard, "iteration: " + util::Format()( count));
          ++count;
        }
      }

      SaxsOptiResult result( lowest_c1, lowest_c2, lowest_chi);
      return result;
    }

    //! @brief one Dimensional unimodal derivative free optimization method using Quadratic Interpolation
    //! @param MIN_C1 - Minimum C1 value
    //! @param MAX_C1 - Maximum C1 value
    //! @param C2_VALUE - Constant row of the matrix
    //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
    //! @return returns storage vector of 3 values: c1, c2, and min chi score
    SaxsOptiResult SasOptimization::QuadraticInterpolation
    (
      float MIN_C1,
      float MAX_C1,
      float C2_VALUE,
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      float a, b, c, d;
      float function_at_a, function_at_b, function_at_c, function_at_d;

      a = MIN_C1;
      b = MAX_C1;
      c = (MIN_C1 + MAX_C1) / 2;

      function_at_a = EvaluateFunction( a, C2_VALUE, PROTEIN_MODEL);
      function_at_b = EvaluateFunction( b, C2_VALUE, PROTEIN_MODEL);
      function_at_c = EvaluateFunction( c, C2_VALUE, PROTEIN_MODEL);

      d = Quadratic( a, function_at_a, b, function_at_b, c, function_at_c);

      // BCL_Message( util::Message::e_Standard, " d is: " + util::Format()( d));

      // handle the case where d is outside of the interval [a, b]
      if( !InsideQuadraticInterval( a, b, d))
      {
        SaxsOptiResult result;

        if( function_at_a < function_at_b)
        {
          result.SetC1( a);
          result.SetC2( C2_VALUE);
          result.SetChi( function_at_a);
        }
        else
        {
          result.SetC1( b);
          result.SetC2( C2_VALUE);
          result.SetChi( function_at_b);
        }
        return result;
      }

      // handle the case where c = d
      if( c == d)
      {
        c = c + 0.001;
      }

      function_at_d = EvaluateFunction( d, C2_VALUE, PROTEIN_MODEL);

      float prev_d( d);
      do
      {
        prev_d = d;
        // f increases on [C,D]
        if( c < d && function_at_c < function_at_d)
        {
          // a = a
          // c = c
          b = d;
          function_at_b = function_at_d;
          d = Quadratic( a, function_at_a, b, function_at_b, c, function_at_c);

          function_at_d = EvaluateFunction( d, C2_VALUE, PROTEIN_MODEL);
        }

        // f decreases on [D,C]
        else if( c > d && function_at_c < function_at_d)
        {
          a = d;
          function_at_a = function_at_d;
          // c = c
          // b = b;
          d = Quadratic( a, function_at_a, b, function_at_b, c, function_at_c);

          function_at_d = EvaluateFunction( d, C2_VALUE, PROTEIN_MODEL);
        }
        else if( c < d && function_at_c > function_at_d)
        {
          a = c;
          function_at_a = function_at_c;

          c = d;
          function_at_c = function_at_d;

          // b = b;
          d = Quadratic( a, function_at_a, b, function_at_b, c, function_at_c);

          function_at_d = EvaluateFunction( d, C2_VALUE, PROTEIN_MODEL);
        }
        else
        {
          // a = a;

          b = c;
          function_at_b = function_at_c;

          c = d;
          function_at_c = function_at_d;

          // b = b;
          d = Quadratic( a, function_at_a, b, function_at_b, c, function_at_c);

          function_at_d = EvaluateFunction( d, C2_VALUE, PROTEIN_MODEL);
        }

      } while( InsideQuadraticInterval( MIN_C1, MAX_C1, d) && ( d - prev_d >= 0.002));

      SaxsOptiResult result;

      if( function_at_c < function_at_d)
      {
        result.SetC1( c);
        result.SetC2( C2_VALUE);
        result.SetChi( function_at_c);
      }
      else
      {
        result.SetC1( d);
        result.SetC2( C2_VALUE);
        result.SetChi( function_at_d);

      }
      return result;
    }

    //! @brief one Dimensional unimodal derivative free optimization method using GoldenSection Method
    //! @param MIN_C1 - Minimum C1 value
    //! @param MAX_C1 - Maximum C1 value
    //! @param C2_VALUE - Constant row of the matrix
    //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
    //! @return returns storage vector of 3 values: c1, c2, and min chi score
    SaxsOptiResult SasOptimization::GoldenSection
    (
      float MIN_C1,
      float MAX_C1,
      float C2_VALUE,
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // This is an implementation of the Golden Section algorithm
      // first initialize all variables;

      // golden ratio (sqrt(5) - 1) / 2
      float theta( 0.618033989);

      float a( MIN_C1);
      float function_at_a( EvaluateFunction( a, C2_VALUE, PROTEIN_MODEL));

      float b( MAX_C1);
      float function_at_b( EvaluateFunction( b, C2_VALUE, PROTEIN_MODEL));

      float c( theta * a + ( 1 - theta) * b);
      float function_at_c( EvaluateFunction( c, C2_VALUE, PROTEIN_MODEL));

      float d( ( 1 - theta) * a + theta * b);
      float function_at_d( EvaluateFunction( d, C2_VALUE, PROTEIN_MODEL));

      storage::Vector< float> storage_function_values
      (
        storage::Vector< float>::Create( function_at_a, function_at_b, function_at_c, function_at_d)
      );

      float prev_min_value( GetMinimumValue( storage_function_values));
      float current_min_value( prev_min_value);

      do
      {
        prev_min_value = current_min_value;

        if( function_at_c < function_at_d)
        {
          // a remains the same
          b = d;
          function_at_b = function_at_d;

          d = c;
          function_at_d = function_at_c;

          c =    theta * a + ( 1 - theta)  * b;
          function_at_c = EvaluateFunction( c, C2_VALUE, PROTEIN_MODEL);

          storage::Vector< float> c_less_than_d
          (
            storage::Vector< float>::Create( function_at_a, function_at_b, function_at_c, function_at_d)
          );

          current_min_value = GetMinimumValue( c_less_than_d);
        }
        else
        {
          // b remains the same
          a = c;
          function_at_a = function_at_c;

          c = d;
          function_at_c = function_at_d;

          d = ( 1 - theta) * a + theta * b;
          function_at_d = EvaluateFunction( d, C2_VALUE, PROTEIN_MODEL);

          storage::Vector< float> c_greater_equal_d
          (
            storage::Vector< float>::Create( function_at_a, function_at_b, function_at_c, function_at_d)
          );

          current_min_value = GetMinimumValue( c_greater_equal_d);
        }
      } while( math::Absolute( double( c - a)) > ( 2 * m_c1StepSize));

      //while ( math::Absolute( double( current_min_value - prev_min_value)) > m_c1StepSize);

      storage::Vector< SaxsOptiResult> result_vector;
      storage::Vector< float> values;

      SaxsOptiResult result_a( a, C2_VALUE, function_at_a);
      result_vector.PushBack( result_a);
      values.PushBack( function_at_a);

      SaxsOptiResult result_b( b, C2_VALUE, function_at_b);
      result_vector.PushBack( result_b);
      values.PushBack( function_at_b);

      SaxsOptiResult result_c( c, C2_VALUE, function_at_c);
      result_vector.PushBack( result_c);
      values.PushBack( function_at_c);

      SaxsOptiResult result_d( d, C2_VALUE, function_at_d);
      result_vector.PushBack( result_d);
      values.PushBack( function_at_d);

      size_t index( GetMinimumIndex( values));

      return result_vector( index);
    }

    //! @brief returns lowest chi score from a search pattern starting at quadratic interpolation line search
    //! @param RADIUS, search radius
    //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
    //! @return SaxsOptiResult: c1, c2, and min chi score
    SaxsOptiResult SasOptimization::HookJeeves
    (
      linal::VectorND< float, 2> INITIAL_POINT,
      float RADIUS,
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      //restraint::SaxsOptiResult start_point( QuadraticInterpolation( m_c1Min, m_c1Max, m_c2Min, PROTEIN_MODEL));

      linal::VectorND< float, 2> grid_point( DiscretePoint( INITIAL_POINT));

      storage::Vector< linal::VectorND< float, 2> > x_coordinates;
      x_coordinates.PushBack( grid_point);
      linal::VectorND< float, 2> x;
      float h( RADIUS);
      size_t k( 0);

      float epsilon( 1.0);
      while( h >= epsilon)
      {
        bool found( false);

        while( found == false && h >= epsilon)
        {
          x = Explore( x_coordinates( k), h, PROTEIN_MODEL);
          if( x != x_coordinates( k))
          {
            x_coordinates.PushBack( x);
            k = k + 1;
            found = true;
          }
          else
          {
            h = h / 2;
          }
        }
      }

      linal::VectorND< float, 2> final_x = x_coordinates.LastElement();

      float min_value( EvaluatePoint( final_x, PROTEIN_MODEL));

      SaxsOptiResult result( final_x.First(), final_x.Last(), min_value);
      return result;
    }

    // helper functions

    //! @brief Convert C1 coordinates to Matrix position
    //! @param C1 - value of contrast density parameter
    //! @return matrix location of C1 value
    int SasOptimization::C1ParamtoMatrix( float C1) const
    {
      float x( ( C1 - m_c1Min) / m_c1StepSize);
      return x > 0 ? ( x + 0.5) : ( x - 0.5);
    }

    //! @brief Convert C2 coordinates to Matrix position
    //! @param C2 - value of contrast density parameter
    //! @return matrix location of C2 value
    int SasOptimization::C2ParamtoMatrix( float C2) const
    {
      float x( ( C2 - m_c2Min) / m_c2StepSize);
      return x > 0 ? ( x + 0.5) : ( x - 0.5);
    }

    //! @brief takes input point from line search method and rounds it to a valid start point on the grid
    //! @param INITIAL_POINT point from line search that is not on a grid point
    //! @return point on Grid
    linal::VectorND< float, 2> SasOptimization::DiscretePoint( const linal::VectorND< float, 2> INITIAL_POINT)const
    {
      float c1( INITIAL_POINT.First());
      float c2( INITIAL_POINT.Last());

      int c1_int( c1 * 100);
      int c2_int( c2 * 10);

      float c1_grid( (float)c1_int / 100);
      float c2_grid( (float)c2_int / 10);
      linal::VectorND< float, 2> final_point;
      final_point(0) = c1_grid;
      final_point(1) = c2_grid;
      return final_point;
    }

    //! @brief compute chi for given parameters and protein model
    //! @param C1 - electron density adjustment parameter for form factors
    //! @param C2 - hydration shell adjustmet parameter for form factors
    //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
    //! @return chi score between experimental and computed saxs profiles
    float SasOptimization::EvaluateFunction( float C1, float C2, const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      m_numEvaluations++;
      BCL_Message
      (
        util::Message::e_Standard,
           " evaluation number: " + util::Format()( m_numEvaluations)
      );

      // Setup Commandline Strings for either the opencl or non-opencl version of the code
      util::Implementation< SasDebyeInterface> sas
      (
        SasAnalysis::SetDebyeImplementation
        ( false, false, C1, C2, m_Cpu, m_Sans, m_DeuteriumExchangeParameter)
      );

      sas->SetExperimentalData( m_ExpData);

      SasExperimentalAndCalculatedData data_sets( sas->operator()( PROTEIN_MODEL));
      SasExperimentalAndCalculatedData transformed_data( m_transformation( data_sets));

      float score = ( float) score::SasType( m_useErrors, m_scoreType)( transformed_data);

      BCL_MessageStd( "C1: " + util::Format()( C1) + " C2: " + util::Format()( C2) + " Score: " + util::Format()( score));

      return score;
    }

    //! @brief returns lowest chi score all points around the center point based on Radius
    //! @param POINT, center point consisting of c1 and c2
    //! @param RADIUS, search radius
    //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
    //! @return linal::VectorND< float, 2>: c1, c2 coordinates of min function value
    linal::VectorND< float, 2> SasOptimization::Explore
    (
      const linal::VectorND< float, 2> &POINT,
      float RADIUS,
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // These are the direction unit vectors
      linal::VectorND< float, 2> e1;
      e1( 0) = m_c1StepSize;
      e1( 1) = 0;

      linal::VectorND< float, 2> e2;
      e2( 0) = 0;
      e2( 1) = m_c2StepSize;

      // These are the diagonal directions
      linal::VectorND< float, 2> e3;
      e3( 0) = m_c1StepSize;
      e3( 1) = m_c2StepSize;

      linal::VectorND< float, 2> e4;
      e4( 0) = -m_c1StepSize;
      e4( 1) = -m_c2StepSize;

      linal::VectorND< float, 2> e5;
      e5( 0) = m_c1StepSize;
      e5( 1) = -m_c2StepSize;

      linal::VectorND< float, 2> e6;
      e6( 0) = -m_c1StepSize;
      e6( 1) = m_c2StepSize;

      float h( RADIUS);

      linal::VectorND< float, 2> x0 = POINT;

      // pick initial evaluation points
      linal::VectorND< float, 2> x2, x3, x4, x5, x6, x7, x8, x9;
      x2 = x0 + h * e1;
      x3 = x0 + h * e2;
      x4 = x0 - h * e1;
      x5 = x0 - h * e2;

      x6 = x0 + h * e3;
      x7 = x0 + h * e4;
      x8 = x0 + h * e5;
      x9 = x0 + h * e6;

      storage::Vector< linal::VectorND< float, 2> > search_points;
      search_points.PushBack( x0);
      search_points.PushBack( x2);
      search_points.PushBack( x3);
      search_points.PushBack( x4);
      search_points.PushBack( x5);
      search_points.PushBack( x6);
      search_points.PushBack( x7);
      search_points.PushBack( x8);
      search_points.PushBack( x9);

      storage::Vector< float> function_values;

      for
      (
        storage::Vector< linal::VectorND< float, 2> >::const_iterator
          itr_search_points( search_points.Begin()),
          itr_search_points_end( search_points.End());
        itr_search_points != itr_search_points_end;
        ++itr_search_points
      )
      {
        if( Inbounds( *itr_search_points))
        {
          float value = EvaluatePoint( *itr_search_points, PROTEIN_MODEL);
          function_values.PushBack( value);
        }
        else
        {
          function_values.PushBack( 99999);
        }
      }

      size_t min_index( GetMinimumIndex( function_values));

      return search_points( min_index);
    }

    //! @brief returns the function value of a given c1 and c2 combination.  Results are stored in m_searchGrid to
    //! @brief prevent duplication of work
    //! @param POINT, center point consisting of c1 and c2
    //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
    //! @return function value of given c1 and c2 parameter
    float SasOptimization::EvaluatePoint
    (
      const linal::VectorND< float, 2> &POINT,
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      float function_at_x = util::GetUndefined< float>();

      float c1( POINT.First());
      float c2( POINT.Last());

      int c1_matrix = C1ParamtoMatrix( c1);
      int c2_matrix = C2ParamtoMatrix( c2);

      if( m_searchGrid( c2_matrix, c1_matrix) == -100.0)
      {
        function_at_x = EvaluateFunction( c1, c2, PROTEIN_MODEL);
        m_searchGrid( c2_matrix, c1_matrix) = function_at_x;
      }
      else
      {
        function_at_x = m_searchGrid( c2_matrix, c1_matrix);
      }

      return function_at_x;
    }

    //! @brief find index of minimum value in storage vector
    //! @param ARRAY storage vector to search
    //! @return index of minimum value in array
    size_t SasOptimization::GetMinimumIndex( const storage::Vector< float> &ARRAY) const
    {
      size_t size( ARRAY.GetSize());
      size_t index( 0);

      for( size_t i( 1); i < size; ++i)
      {
        if( ARRAY( i) < ARRAY( index))
        {
          index = i;
        }
      }

      return index;
    }

    //! @brief return minimum value in storage vector
    //! @param ARRAY storage vector to search
    //! @return minimum value in array
    float SasOptimization::GetMinimumValue( const storage::Vector< float> &ARRAY) const
    {
      linal::Vector< float> linal_function_values( ARRAY);
      return linal_function_values.Min();
    }

    //! @brief Ensure that the proposed point to evaluate is inside the searchGrid
    //! @param POINT, center point consisting of c1 and c2
    //! @return true if inbounds
    bool SasOptimization::Inbounds( linal::VectorND< float, 2> POINT) const
    {
      bool inbounds( true);

      if( POINT.First() < m_c1Min)
      {
        inbounds = false;
      }
      else if( POINT.First() > m_c1Max)
      {
        inbounds = false;
      }

      else if( POINT.Last() < m_c2Min)
      {
        inbounds = false;
      }
      else if( POINT.Last() > m_c2Max)
      {
        inbounds = false;
      }

      return inbounds;
    }

    //! @brief Ensure that the proposed point to evaluate (D) is inside the specified interval
    //! @param A - Left Boundary Point
    //! @param B - Right Boundary Point
    //! @param TEST_POINT, point to evaluate
    //! @return true if inbounds
    bool SasOptimization::InsideQuadraticInterval( float A, float B, float TEST_POINT) const
    {
      bool inbounds( true);
      if( TEST_POINT < A)
      {
        inbounds = false;
      }
      else if( TEST_POINT > B)
      {
        inbounds = false;
      }
      return inbounds;
    }

    //! @brief computes a quadratic model based on 3 values and thier function
    //! @param A - First point
    //! @param FA - function at first point
    //! @param B - Second point
    //! @param FB - function at second point
    //! @param C - Third point
    //! @param FC - function at third point
    float SasOptimization::Quadratic( float A, float FA, float B, float FB, float C, float FC) const
    {
      float calculation
      (
        0.5
        * ( FA * ( B*B - C*C ) - FC * ( B*B - A*A) + FB * ( C*C - A*A))
        / ( FA * ( B - C)      - FC * ( B - A)     + FB * ( C - A))
       );

      return calculation;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SasOptimization::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_c1Min, ISTREAM);
      io::Serialize::Read( m_c1Max, ISTREAM);
      io::Serialize::Read( m_c2Min, ISTREAM);
      io::Serialize::Read( m_c2Max, ISTREAM);
      io::Serialize::Read( m_c1StepSize, ISTREAM);
      io::Serialize::Read( m_c2StepSize, ISTREAM);
      io::Serialize::Read( m_numEvaluations, ISTREAM);
      io::Serialize::Read( m_searchGrid, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SasOptimization::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_c1Min, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_c1Max, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_c2Min, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_c2Max, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_c1StepSize, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_c2StepSize, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_numEvaluations, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_searchGrid, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_sas_pofr.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_ca_cb.h"
#include "biol/bcl_biol_atom_group_types.h"
#include "biol/bcl_biol_sasa_data.h"
#include "fold/bcl_fold_add_parabolic_loops.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_histogram.h"
#include "math/bcl_math_sum_function.h"
#include "restraint/bcl_restraint_sas_experimental_and_calculated_density.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> SasPofR::s_Instance
    (
      util::Enumerated< SasPofRInterface>::AddInstance( new SasPofR())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new SasPofR function
    SasPofR *SasPofR::Clone() const
    {
      return new SasPofR( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &SasPofR::GetAlias() const
    {
      static const std::string s_Name( "SasPofR");
      return s_Name;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SasPofR::GetClassIdentifier() const
    {
      // Get BCL Standardized Class name
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Get the positions for all the atoms in the protein model
    //! @param MODEL the protein model of interest
    //! @return a vector containing positions
    storage::Vector< linal::Vector3D> SasPofR::GetAtoms( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      BCL_MessageDbg( " Inside Get all atom coordinates from the protein model: ");

      // initialize vector to hold atom coordinates
      storage::Vector< linal::Vector3D> all_atoms;

      size_t atom_number( 0);

      // iterate over the chains
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
         chain_itr( PROTEIN_MODEL.GetChains().Begin()), chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end; ++chain_itr
      )
      {
        // get all the structured sses
        const util::SiPtrVector< const assemble::SSE> structured_sses
        (
          ( *chain_itr)->GetSSEs( storage::Set< biol::SSType>( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND))
        );

        // if structured_sses are empty
        if( structured_sses.IsEmpty())
        {
          // warn user and return empty vector
          BCL_Message( util::Message::e_Standard, "No structured SSEs found in protein model");
          continue;
        }

        // iterate over all the SSEs including coil (depending on min_sse_size parameter)
        for
        (
          storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
           sse_itr( ( *chain_itr)->GetData().Begin()), sse_itr_end( ( *chain_itr)->GetData().End());
          sse_itr != sse_itr_end; ++sse_itr
        )
        {
          // iterate over the Amino Acids in a given SSE
          for
          (
            util::ShPtrVector< biol::AABase>::const_iterator aa_itr( ( *sse_itr)->GetData().Begin()),
             aa_itr_end( ( *sse_itr)->GetData().End());
            aa_itr != aa_itr_end; ++aa_itr
          )
          {
            std::string amino_acid( ( *aa_itr)->GetType()->GetName());
            BCL_MessageDbg( "Residue: " + util::Format()( amino_acid));

            // get the atom types for this residue
            // Be sure to get the Heavy Atom Types ( all atoms except hydrogen).  This set contains all of the residues
            // That a given amino acid contains.  Not the atoms actually present in the input PDB file
            const storage::Set< biol::AtomType> &complete_atom_types( ( *aa_itr)->GetType()->GetAllowedHeavyAtomTypes());

            BCL_MessageDbg( "AtomTypes: " + util::Format()( complete_atom_types));

            // iterate over the atoms in a given residue
            for
            (
              storage::Set< biol::AtomType>::const_iterator atom_itr( complete_atom_types.Begin()),
               atom_itr_end( complete_atom_types.End());
              atom_itr != atom_itr_end; ++atom_itr
            )
            {
              atom_number++;

              // get the atom coordinates
              const linal::Vector3D &atom_coords( ( *aa_itr)->GetAtom( *atom_itr).GetCoordinates());

              BCL_MessageDbg( "Atom Coords: " + util::Format()( atom_coords));

              if( atom_coords.IsDefined())
              {
                std::string atom( ( *atom_itr)->GetName());
              }

              // if the coordinates are defined
              if( atom_coords.IsDefined())
              {
                // pushback into vector
                all_atoms.PushBack( atom_coords);
              }
            } // Atom Type Iteration
          } // AA iteration
        } // SSE iteration
      } // chain iteration
      return all_atoms;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes a proteinModel as input and calculates an intensity using the debye formula
    //! @param PROTEIN_MODEL
    //! @return the intensity for a given q value
    SasExperimentalAndCalculatedDensity SasPofR::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // This holds the Euclidean coordinates of the atoms in the model
      storage::Vector< linal::Vector3D> all_atoms( GetAtoms( PROTEIN_MODEL));

      BCL_MessageDbg( " All atoms: " + util::Format()( all_atoms));

      double bin_size( this->GetExperimentalDensity()->GetBinSize());
      size_t bin_number( this->GetExperimentalDensity()->GetBinNumber());

      // construct Histogram from starting bin, bin size and number of bins
      math::Histogram distance_histogram( 0.0, bin_size, bin_number);

      size_t row( 0);
      double dmax( 0);

      // Compute the upper triangle of pairwise distances and increment the historgram
      for
      (
        storage::Vector< linal::Vector3D>::const_iterator atom_itr_a( all_atoms.Begin()), atom_itr_end( all_atoms.End());
        atom_itr_a != atom_itr_end; ++atom_itr_a, ++row
      )
      {
        size_t col( row + 1);
        for
        (
          storage::Vector< linal::Vector3D>::const_iterator atom_itr_b( atom_itr_a + 1);
          atom_itr_b != atom_itr_end; ++atom_itr_b, ++col
        )
        {

          double distance = linal::Distance( *atom_itr_a, *atom_itr_b);
          distance_histogram.PushBack( distance);
          if( dmax < distance)
          {
            dmax = distance;
          }

          // Take care of the lower half of the matrix
          distance_histogram.PushBack( distance);
        }
      }

      // create object to hold calculated data
      // The computed dmax is set here, the experimental dmax is read in at object creation
      SasDensityData calculated_data( distance_histogram, dmax);

      SasExperimentalAndCalculatedDensity pofr_data( *this->GetExperimentalDensity(), calculated_data);

      BCL_MessageDbg( "data: " + util::Format()( pofr_data));

      double cal_dmax( pofr_data.GetCalculatedDensity().GetDmax());
      double exp_dmax( pofr_data.GetExperimentalDensity().GetDmax());

      BCL_MessageDbg( "cal_dmax: " + util::Format()( cal_dmax));
      BCL_MessageDbg( "exp_dmax: " + util::Format()( exp_dmax));

      double cal_hmax( pofr_data.GetCalculatedDensity().GetHmax());
      double exp_hmax( pofr_data.GetExperimentalDensity().GetHmax());

      BCL_MessageDbg( "cal_hmax: " + util::Format()( cal_hmax));
      BCL_MessageDbg( "exp_hmax: " + util::Format()( exp_hmax));

      return pofr_data;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SasPofR::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "performs pairwise atomic distance calculation"
      );

      return parameters;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_sas_scattering_data.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "math/bcl_math.h"
#include "restraint/bcl_restraint_sas_analysis.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SasScatteringData::s_Instance
    (
      util::Enumerated< HandlerBase< SasScatteringData> >::AddInstance
      (
        new SasScatteringData()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SasScatteringData::SasScatteringData( const std::string &EXTENSION) :
      HandlerBase< SasScatteringData>( EXTENSION),
      m_Data()
    {
    }

    //! @brief constructor from given input data
    SasScatteringData::SasScatteringData
    (
      const storage::Vector< SasScatteringPoint> &INIT_DATA
    ) :
      HandlerBase< SasScatteringData>( ".saxs"),
      m_Data( INIT_DATA)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SasScatteringData
    SasScatteringData *SasScatteringData::Clone() const
    {
      return new SasScatteringData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SasScatteringData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns short name for the class used over the command line
    //! @return short name for the class used over the command line
    const std::string &SasScatteringData::GetAlias() const
    {
      // format is automatically determined
      static const std::string s_name( "Detect");
      return s_name;
    }

    //! @brief pushback function to add object to Dataset vector
    //! @param DATAPOINT_OBJECT DataPoint values Q, I, and Error
    void SasScatteringData::PushBackScattering( const SasScatteringPoint &DATAPOINT_OBJECT)
    {
      m_Data.PushBack( DATAPOINT_OBJECT);
    }

  /////////////////
  // operations  //
  /////////////////

    //! @brief preallocate Scattering Data memory
    //! @param SIZE size to preallocate
    void SasScatteringData::AllocateScatteringMemory( const size_t &SIZE)
    {
      m_Data.AllocateMemory( SIZE);
    }

    //! @brief bool test for error values
    //! @return true if error is defined for all values of the dataset
    const bool SasScatteringData::IsErrorDefined() const
    {
      // iterate over m_Data
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator itr( m_Data.Begin()), itr_end( m_Data.End());
         itr != itr_end; ++itr
      )
      {
        // if error is defined do not terminate
        if( !itr->IsErrorDefined())
        {
          // if error is not defined return false immediately
          return false;
        }
      }

      // error is defined for all cases, return true
      return true;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read experimental data from BCL
    //! @param ISTREAM input data stream
    void SasScatteringData::ReadBCLProcessedData( std::istream &ISTREAM)
    {
      // String to hold line data
      std::string read_line;

      // iterate over each line in the file
      while( std::getline( ISTREAM, read_line).good())
      {
        // split the string
        const storage::Vector< std::string> split_line( util::SplitString( read_line, " "));

        double q( util::ConvertStringToNumericalValue< double>( split_line( 0)));
        double intensity( util::ConvertStringToNumericalValue< double>( split_line( 1)));
        double error( util::ConvertStringToNumericalValue< double>( split_line( 2)));

        // push back the data
        m_Data.PushBack( SasScatteringPoint( q, intensity, error));
      } // close while loop
    } // close function

    //! @brief read experimental data from Crysol. See Crysol documentation for column definitions
    //! @param ISTREAM input data stream
    void SasScatteringData::ReadCrysolData( std::istream &ISTREAM)
    {
      // String to hold line data
      std::string read_line;

      // iterate over each line in the file
      while( std::getline( ISTREAM, read_line).good())
      {
        // split the string
        const storage::Vector< std::string> split_line( util::SplitString( read_line, " "));

        double q( util::ConvertStringToNumericalValue< double>( split_line( 0)));
        double intensity( util::ConvertStringToNumericalValue< double>( split_line( 1)));

        m_Data.PushBack( SasScatteringPoint( q, intensity, SasAnalysis::ComputeError( q, intensity)));
      } // close while loop

    } // close function

    //! @brief read experimental data
    //! @param ISTREAM input data stream
    void SasScatteringData::ReadExperimentalData( std::istream &ISTREAM)
    {
      //String to hold line data
      std::string read_line;
      while( std::getline( ISTREAM, read_line).good())
      {
        // split the string
        const storage::Vector< std::string> split_line( util::SplitString( read_line, " "));

        double intensity( util::ConvertStringToNumericalValue< double>( split_line( 1)));
        if( intensity > 0)
        {
          // push back the data
          m_Data.PushBack
          (
            SasScatteringPoint
            (
              util::ConvertStringToNumericalValue< double>( split_line( 0)),    // Q-value
              util::ConvertStringToNumericalValue< double>( split_line( 1)),    // I-value for intensity in solution
              util::ConvertStringToNumericalValue< double>( split_line( 2))     // error
            )
          );
        }
        else
        {
           double value( util::ConvertStringToNumericalValue< double>( split_line( 1)));
           BCL_MessageStd( "Improper intensity value: " + util::Format()( value));
           BCL_Exit( "Intensity values cannot be zero or negative please examine input experimental data: ", -1);
        }
      } // close while loop
    } // close function

    //! @brief read experimental data from Gnom
    //! @param ISTREAM input data stream
    void SasScatteringData::ReadGnomData( std::istream &ISTREAM, size_t &INTENSITY_COLUMN, bool USE_EXTRAPOLATION)
    {

      std::string read_line;

      // read first line of input file
      std::getline( ISTREAM, read_line);

      // split the string based on white space as a delimeter
      storage::Vector< std::string> Line( util::SplitString( read_line, " "));

      // get the size of the first line
      size_t line_size( Line.GetSize());

      // paramaters to read through gnom filetype
      bool readFlag( false);

      // while the end of the file is not reached
      while( std::getline( ISTREAM, read_line).good())
      {
        // split the string
        const storage::Vector< std::string> split_line( util::SplitString( read_line, " "));
        if
        (
          split_line.GetSize() >= 7 &&
          split_line( 0) == "S"     &&
          split_line( 1) == "J"     &&
          split_line( 2) == "EXP"   &&
          split_line( 3) == "ERROR" &&
          split_line( 4) == "J"     &&
          split_line( 5) == "REG"   &&
          split_line( 6) == "I"
        )
        {
          readFlag = true;
        }
        else if( readFlag)
        {
          // get size of the vector
          line_size = split_line.GetSize();

          if( line_size > 1 && !util::IsNumerical( split_line( 0)))
          {
            readFlag = false;
          }
          else if( line_size == 2 && USE_EXTRAPOLATION)
          {
            // push back the extrapolated data
             m_Data.PushBack
             (
               SasScatteringPoint
               (
                 util::ConvertStringToNumericalValue< double>( split_line( 0)),  // Q-value
                 util::ConvertStringToNumericalValue< double>( split_line( 1)),  // I-value
                 util::GetUndefined< double>()                                   // error
               )
             );
          }
          else if( line_size == 5)
          {
            // push back the data
            m_Data.PushBack
            (
              SasScatteringPoint
              (
                util::ConvertStringToNumericalValue< double>( split_line( 0)),                 // Q-value
                util::ConvertStringToNumericalValue< double>( split_line( INTENSITY_COLUMN)),  // I-value
                util::ConvertStringToNumericalValue< double>( split_line( 2))                  // error
              )
            );
          }
        }
      } // close while loop
    } // close the function

    //! @brief reads in the member data from a formatted file containing 3 columns:
    //! @brief scattering angle q ( 4*Pi*sin(theta)/lambda) where lambda is measured in Angstroms, I is the intensity
    //! @brief at a given a value, E is the experimental error.
    //! @brief Crysol generated files must used the paramater /dro 0.0.  The algorithm does not support adding the
    //! @brief hydration layer around the molecule.
    //! @param ISTREAM input stream
    //! @param FORMAT the file format to use for reading
    //! @return istream which was read from
    std::istream &SasScatteringData::ReadFromDataFile( std::istream &ISTREAM)
    {
      // build a string to hold the line information
      std::string read_line;

      enum Types { bcl_processed_data, crysol, experimental, gnom, unknown};
      Types filetype( unknown);

      // read first line of input file
      std::getline( ISTREAM, read_line);

      // split the string based on white space as a delimeter
      storage::Vector< std::string> Line( util::SplitString( read_line, " "));

      // get the size of the first line
      size_t line_size( Line.GetSize());

      // ignore empty lines at the first of the file
      while( line_size <= 1)
      {
        std::getline( ISTREAM, read_line);
        Line = util::SplitString( read_line, " ");
        line_size = Line.GetSize();
      }

      // verify bclfile type
      if(
          line_size == 7 &&
          Line( 0) == "BCL" &&
          Line( 1) == "SCATTERING" &&
          Line( 2) == "INPUT" &&
          Line( 3) == "PARAMETERS:" &&
          Line( 4) == "QValue" &&
          Line( 5) == "Intensity" &&
          Line( 6) == "Error"
        )
      {
        BCL_MessageStd( "bcl_processed_data");
        filetype = bcl_processed_data;
      }

      else if(
          line_size == 4 &&
          Line( 0) == "Q_Value" &&
          Line( 1) == "Experimental_Intensity" &&
          Line( 2) == "Experimental_Error" &&
          Line( 3) == "Computed_Intensity"
        )
      {
        BCL_MessageStd( "bcl_processed_data");
        filetype = bcl_processed_data;
      }

      // verify crysol file type
      else if( Line( 0) == "Dif/Atom/Shape/Bord")
      {
        filetype = crysol;
      }

      // Verify the file is a GNOM file
      else if( line_size >= 5 && Line( 1) == "G" && Line( 2) == "N" && Line( 3) == "O" && Line( 4) == "M")
      {
        filetype = gnom;
      }

      // experimental data file type
      else if( line_size == 3)
      {
        filetype = experimental;
        double intensity( util::ConvertStringToNumericalValue< double>( Line( 1)));
        if( intensity > 0)
        {
          // push back the data
          m_Data.PushBack
          (
            SasScatteringPoint
            (
               util::ConvertStringToNumericalValue< double>( Line( 0)),    // Q-value
               util::ConvertStringToNumericalValue< double>( Line( 1)),    // I-value for intensity in solution
               util::ConvertStringToNumericalValue< double>( Line( 2))     // error
            )
          );
        }
        else
        {
          BCL_MessageStd( "Skipped negative or zero intensity value");
        }

      }
      else
      {
        BCL_Exit( "Could not read file; last read line was: " + read_line, -1);
      }

      size_t experimental_data_column( 1);
      switch( filetype)
      {
        case bcl_processed_data:
          ReadBCLProcessedData( ISTREAM);
          break;
        case crysol:
          ReadCrysolData( ISTREAM);
          break;
        case experimental:
          ReadExperimentalData( ISTREAM);
          break;
        case gnom:
          ReadGnomData( ISTREAM, experimental_data_column);
          break;
        case unknown:
        default:
          break;
      }

      size_t size( 0);

      // make sure there are q values for analysis
      size = m_Data.GetSize();

      BCL_Assert( size != 0, "The number of Q values are: " + util::Format()( size));

      // return the stream
      return ISTREAM;
    }

    //! @brief reads saxs restraints from an input stream
    //! @brief input stream to read the restraints from
    //! @return the read in restraints
    SasScatteringData SasScatteringData::ReadRestraints( std::istream &ISTREAM) const
    {
      SasScatteringData data;
      data.ReadFromDataFile( ISTREAM);
      return data;
    }

    //! @brief reads the computed SAXS profile from the fit density curve from gnom
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SasScatteringData::ReadFitFromGnom( std::istream &ISTREAM)
    {
      size_t fit_data_column( 4);
      bool use_extrapolation( true);
      ReadGnomData( ISTREAM, fit_data_column, use_extrapolation);

      size_t size( 0);
      size = m_Data.GetSize();
      BCL_Assert( size != 0, "The number of Q values are: " + util::Format()( size));

      // return the stream
      return ISTREAM;
    }

    //! @brief writes out the member data from a formatted file containing 3 columns
    //! @param OSTREAM output stream
    //! @return ostream which was written to
    std::ostream &SasScatteringData::WriteToDataFile( std::ostream &OSTREAM, bool HEADER) const
    {
      // initialize format
      const util::Format format;

      if( HEADER)
      {
        OSTREAM << format( "BCL SCATTERING INPUT PARAMETERS: QValue Intensity Error") << '\n';
      }

      // iterate over m_Data
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator itr( m_Data.Begin()), itr_end( m_Data.End());
         itr != itr_end; ++itr
      )
      {
        // write the data
        OSTREAM << format( itr->GetQvalue()) << ' '
                << format( itr->GetIntensity()) << ' '
                << format( itr->GetError()) << '\n';
      }

      // end
      return OSTREAM;
    }

    //! @brief write to filename in four columns : Q-value, experimental intensity, error and calculated intensity
    //! @param filename to write to
    void SasScatteringData::WriteToDataFileName( const std::string &FILENAME, const bool &HEADER) const
    {
      io::OFStream write;
      io::File::MustOpenOFStream( write, FILENAME);

      // write the data
      WriteToDataFile( write, HEADER);

      // close the ofstream
      io::File::CloseClearFStream( write);
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SasScatteringData::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Data, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SasScatteringData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Data, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_sas_scattering_point.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SasScatteringPoint::s_Instance
    (
      GetObjectInstances().AddInstance( new SasScatteringPoint())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SasScatteringPoint::SasScatteringPoint() :
      m_Qvalue( util::GetUndefinedDouble()),
      m_Intensity( util::GetUndefinedDouble()),
      m_Error( util::GetUndefinedDouble())
    {
    }

    //! @brief constructor from given input data
    //! @param QVALUE the scattering angle from SAXS curve ( x-axis)
    //! @param INTENSITY the wavelength intensity for a given scattering angle ( y-axis)
    //! @param ERROR - the variation in measurement for a given intensity
    SasScatteringPoint::SasScatteringPoint
    (
      const double QVALUE,
      const double INTENSITY,
      const double MEASUREMENT_ERROR
    ) :
       m_Qvalue( QVALUE),
       m_Intensity( INTENSITY),
       m_Error( MEASUREMENT_ERROR)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SaxsScatteringPoint
    SasScatteringPoint *SasScatteringPoint::Clone() const
    {
      return new SasScatteringPoint( *this);
    }

    //! @brief Boolean function to return state of error value
    //! @return return true if the error is defined

    bool SasScatteringPoint::IsErrorDefined() const
    {
      return ( util::IsDefined( m_Error) ? true : false);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SasScatteringPoint::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SasScatteringPoint::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Qvalue, ISTREAM);
      io::Serialize::Read( m_Intensity, ISTREAM);
      io::Serialize::Read( m_Error, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SasScatteringPoint::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Qvalue, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Intensity, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Error, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_sas_transformation.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_histogram.h"
#include "restraint/bcl_restraint_sas_analysis.h"
#include "restraint/bcl_restraint_sas_debye_interface.h"
#include "restraint/bcl_restraint_sas_scattering_data.h"
#include "restraint/bcl_restraint_saxs_opti_result.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

    //! @brief ScoreFunction as string
    //! @param SCORE_FUNCTION the ScoreFunction
    //! @return the string for the ScoreFunction
    const std::string &SasTransformation::GetFunctionDescriptor( const TransformationType &TRANSFORMATION_TYPE)
    {
      static const std::string s_descriptors[] =
      {
        "Absolute",
        "Scale",
        "Normalize",
        "SetYMax",
        "Log10",
        "Derivative",
        "None",
        GetStaticClassName< TransformationType>()
      };

      return s_descriptors[ TRANSFORMATION_TYPE];
    }

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SasTransformation::s_Instance
    (
      util::Enumerated
      <
        util::FunctionInterfaceSerializable
        <
          SasExperimentalAndCalculatedData,
          SasExperimentalAndCalculatedData
        >
      >::AddInstance( new SasTransformation())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SasTransformation::SasTransformation() : m_OutputTransformations( false), m_UseErrors( false)
    {
      storage::Vector< TransformationTypeEnum> default_transforms;

      default_transforms.PushBack( SasTransformation::e_ScaleCalculatedProfile);
      default_transforms.PushBack( SasTransformation::e_SetYScale);
      default_transforms.PushBack( SasTransformation::e_Log10Profiles);
      default_transforms.PushBack( SasTransformation::e_DerivativeProfiles);

      m_Transformations = default_transforms;
      m_Ymax = 1.0;
    }

    //! @brief constructor taking member variables
    //! @param OUTPUT_TRANSFORMATIONS flag to print transformations to a file
    //! @param USE_ERRORS DATA_TYPE flag to use errors
    //! @param Y_MAX the max value for the y scale
    SasTransformation::SasTransformation
    (
      const storage::Vector< TransformationTypeEnum> &TRANSFORMATION_TYPES,
      const bool &PRINT_TRANSFORMATIONS,
      const bool &USE_ERRORS,
      const double &Y_MAX
    )
      :
      m_Transformations( TRANSFORMATION_TYPES),
      m_OutputTransformations( PRINT_TRANSFORMATIONS),
      m_UseErrors( USE_ERRORS),
      m_Ymax( Y_MAX)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AnalyzeAtomDistanceHeatmap
    SasTransformation *SasTransformation::Clone() const
    {
      return new SasTransformation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SasTransformation::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &SasTransformation::GetAlias() const
    {
      static const std::string s_Name( "Transform");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

    void SasTransformation::ScaleCalculatedProfile( SasExperimentalAndCalculatedData &DATA_OBJECT) const
    {
      const double scaling_factor
      (
        SasAnalysis::CalculateScalingWeight( DATA_OBJECT, m_UseErrors)
      );

      DATA_OBJECT.ScaleCalculatedData( scaling_factor);
    }

    void SasTransformation::NormalizeData( SasExperimentalAndCalculatedData &DATA_OBJECT) const
    {
      DATA_OBJECT.NormalizeData();
    }

    void SasTransformation::SetYScale( SasExperimentalAndCalculatedData &DATA_OBJECT) const
    {
      DATA_OBJECT.SetYScale( m_Ymax);
    }

    void SasTransformation::LogProfile( SasExperimentalAndCalculatedData &DATA_OBJECT) const
    {
      DATA_OBJECT.Log10();
    }

    void SasTransformation::DerivativeProfile( SasExperimentalAndCalculatedData &DATA_OBJECT) const
    {
      DATA_OBJECT.Derivative();
    }

    void SasTransformation::LogtoAbsolute( SasExperimentalAndCalculatedData &DATA_OBJECT) const
    {
      DATA_OBJECT.LogtoAbsolute();
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    SasExperimentalAndCalculatedData SasTransformation::operator()
     ( const SasExperimentalAndCalculatedData &ORIGINAL_DATA) const
    {
      SasExperimentalAndCalculatedData data_object( ORIGINAL_DATA);

      if( m_OutputTransformations)
      {
        data_object.WriteToGnuplotFileName( "raw.data");
      }

      // Iterate over the storage vector of Enums
      for
      (
          storage::Vector< TransformationTypeEnum>::const_iterator
            enum_itr( m_Transformations.Begin()),
            enum_itr_end( m_Transformations.End());
            enum_itr != enum_itr_end;
            ++enum_itr
      )
      {
        switch( *enum_itr)
        {
          case e_ScaleCalculatedProfile:
            ScaleCalculatedProfile( data_object);
            if( m_OutputTransformations)
            {
              data_object.WriteToGnuplotFileName( "scaled.data");
            }
            break;
          case e_SetYScale:
            SetYScale( data_object);
            if( m_OutputTransformations)
            {
              data_object.WriteToGnuplotFileName( "scaled_y_axis.data");
            }
            break;
          case e_NormalizeData:
            NormalizeData( data_object);
            if( m_OutputTransformations)
            {
              data_object.WriteToGnuplotFileName( "normalized.data");
            }
            break;
          case e_Log10Profiles:
            LogProfile( data_object);
            if( m_OutputTransformations)
            {
              data_object.WriteToGnuplotFileName( "logarithmic.data");
            }
            break;
          case e_DerivativeProfiles:
            DerivativeProfile( data_object);
            if( m_OutputTransformations)
            {
              data_object.WriteToGnuplotFileName( "derivative.data");
            }
            break;
          case e_Absolute:
            LogtoAbsolute( data_object);
            if( m_OutputTransformations)
            {
              data_object.WriteToGnuplotFileName( "absolute.data");
            }
            break;
          case e_None:
            // Do nothing
            break;
          default:
            BCL_Assert
            (
              false,
              "Unknown Transformation Procedure"
            );
            break;
        }
      }
      return data_object;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SasTransformation::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Transforms Computed SAXS Profile and Experimental Saxs Profile into different forms based on user input."
      );

      SasTransformation default_object;

      parameters.AddInitializer
      (
        "transformations",
        "computed scaling factor between experimental and calculated profiles and scale calculated profile",
        io::Serialization::GetAgent( &m_Transformations),
        io::Serialization::GetAgent( &default_object.m_Transformations)->GetLabel().ToString()
      );

      parameters.AddInitializer
      (
        "print_transformations",
        "Print transformations performed",
        io::Serialization::GetAgent( &m_OutputTransformations),
        "false"
      );

      parameters.AddInitializer
      (
        "use_errors",
        "use experimental errors in scaling saxs profiles",
        io::Serialization::GetAgent( &m_UseErrors)
      );

      parameters.AddInitializer
      (
        "y_max",
        "set max value on y scale for graphing",
        io::Serialization::GetAgent( &m_Ymax),
        "-9999"
      );

      return parameters;
    }
  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_saxs_data.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_mutate_tree.h"
#include "fold/bcl_fold_score_weight_set.h"
#include "restraint/bcl_restraint_sas_analysis.h"
#include "restraint/bcl_restraint_sas_debye.h"
#include "score/bcl_score_restraint_saxs.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    // initialize score
    fold::Score SaxsData::e_ScoreSaxsRestraint( fold::GetScores().e_Undefined);

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SaxsData::s_Instance
    (
      util::Enumerated< Interface>::AddInstance( new SaxsData())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from extension
    //! @param EXTENSION the extension used to identify files containing data for this restraint
    SaxsData::SaxsData() :
      m_Data(),
      m_Extension( GetDefaultExtension())
    {
    }

    //! @brief Clone function
    //! @return pointer to new SaxsData
    SaxsData *SaxsData::Clone() const
    {
      return new SaxsData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &SaxsData::GetAlias() const
    {
      static const std::string s_name( "SAXS");
      return s_name;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SaxsData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the default file extension
    //! @return the default file extension
    const std::string &SaxsData::GetDefaultExtension() const
    {
      static const std::string s_extension( ".saxs");
      return s_extension;
    }

    //! @brief gives reference to the specific type of data this restraint uses
    //! @return gives reference to the specific type of data this restraint uses
    const SasScatteringData &SaxsData::GetScatteringData() const
    {
      return *m_Data;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initialize the scores and add them to Scores enumerator
    void SaxsData::InitializeScores()
    {
      if( !e_ScoreSaxsRestraint.IsDefined())
      {
        // read the restraints from the file, if they aren't already defined
        if( !m_Data.IsDefined())
        {
          // reset the data
          m_Data = util::ShPtr< SasScatteringData>( new SasScatteringData());
          *m_Data = m_Data->ReadRestraintsFromFile();
        }

        bool approximate_loops( true);
        bool approximate_side_chains( true);
        double c1( 1.0);
        double c2( 0.0);
        bool cpu( false);
        bool sans( false);
        double deuterium_exchange( 0.0);

        // Setup Commandline Strings for either the opencl or non-opencl version of the code
        util::Implementation< SasDebyeInterface> saxs
        (
          SasAnalysis::SetDebyeImplementation
          (
            approximate_loops,
            approximate_side_chains,
            c1,
            c2,
            cpu,
            sans,
            deuterium_exchange
          )
        );

        //SasDebye saxs( true);
        //saxs.SetExperimentalData( m_Data);

        saxs->SetExperimentalData( m_Data);
        e_ScoreSaxsRestraint = fold::GetScores().AddScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::RestraintSaxs( saxs, score::SasType())
          )
        );
      }
    }

    //! @brief sets the weights of scores in a weight set
    //! @param SCORE_WEIGHT_SET the score weight set that will be modified
    void SaxsData::ModifyScoreWeightSet( fold::ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
      SCORE_WEIGHT_SET.SetWeight( e_ScoreSaxsRestraint, 500);
    }

    //! @brief modify the mutate tree used
    //! @param MUTATE_TREE MutateTree to be modified
    void SaxsData::ModifyMutateTree( fold::MutateTree &MUTATE_TREE) const
    {
      // set SSE remove to 0
      MUTATE_TREE.SetMutateTypeProbability( fold::MutateTree::e_Remove, 0);
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void SaxsData::MergeAndModifyMutateTree( fold::MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< fold::MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }
 
    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SaxsData::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        "Allows use of small angle x-ray scattering data restraints. "
      );
      serializer.AddInitializer
      (
        "extension",
        "restraints will be read in from {path/prefix given by -restraint_prefix}{extension}",
        io::Serialization::GetAgent( &m_Extension),
        GetDefaultExtension()
      );
      return serializer;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads restraints formatted for this restraint type from an istream
    //! @return istream restraints formatted for this restraint type were read from
    std::istream &SaxsData::ReadRestraints( std::istream &ISTREAM)
    {
      // reset the data
      m_Data = util::ShPtr< SasScatteringData>( new SasScatteringData( m_Extension));

      // read from the stream
      m_Data->ReadFromDataFile( ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SaxsData::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Data, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SaxsData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Data, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_saxs_data_reduction.h"

// includes from bcl - sorted alphabetically
#include "restraint/bcl_restraint_sas_analysis.h"
#include "restraint/bcl_restraint_sas_debye_interface.h"
#include "restraint/bcl_restraint_sas_transformation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SaxsDataReduction::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  /////////////////
  // operations  //
  /////////////////

    //! @brief compute the shannon numbers for saxs profile generation from experimental data
    //! @param DMAX - Maximum dimension of protein inferred from GNOM
    //! @param QMAX - Maximum momentum transfer value
    size_t SaxsDataReduction::ComputeShannonNumber( const double &DMAX, const double &QMAX)
    {
      return std::ceil( DMAX * QMAX / math::g_Pi);
    }

    //! @brief samples shannon bins n times ( n typically is 1000 - 3000) to produce noise free scattering profile
    //! @param PROTEIN_MODEL - protein model that contains the coordinates for current model
    //! @param NUMBER_OF_ITERATIONS - the number of iterations to compute ( n typically is 1000 - 3000)
    //! @param ORIGINAL_DATA - pointer to the preprocessed experimental data
    util::ShPtr< SasScatteringData> SaxsDataReduction::SasSignalRecovery
    (
      assemble::ProteinModel &PROTEIN_MODEL,
      size_t NUMBER_OF_ITERATIONS,
      const util::ShPtr< SasScatteringData> ORIGINAL_DATA,
      const double DMAX
    )
    {
      // Tasks
      // 1) partition a SAS dataset into ns equal bins for a given dmax
      // compute the number of bins on the saxs interval
      size_t number_of_bins
      (
        ComputeShannonNumber
        (
          DMAX,
          SasAnalysis::ComputeQmax( *ORIGINAL_DATA)
        )
      );

      double first_point( ORIGINAL_DATA->GetScatteringData().FirstElement().GetQvalue());
      double last_point( ORIGINAL_DATA->GetScatteringData().LastElement().GetQvalue());

      // compute the size of the bin on the interval based on the first and last point
      double binsize( ( last_point - first_point) / number_of_bins);

      double right_bin( binsize + first_point);
      int count( 0);

      storage::Vector< int> index_selected_points;

      // identify the index values that contain the desired ranges
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
         data_itr( ORIGINAL_DATA->Begin()),
         data_itr_end( ORIGINAL_DATA->End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        if( data_itr->GetQvalue() >= right_bin)
        {
          if( data_itr->GetQvalue() == right_bin)
          {
            index_selected_points.PushBack( count);
          }
          else
          {
            index_selected_points.PushBack( count - 1);
          }
          right_bin = right_bin + binsize;
        }
        ++count;
      }

      int left_index = 1;
      int right_index = index_selected_points( 0);

      // Set up containers to hold the data
      storage::Vector< SasScatteringPoint> subset;
      storage::Vector< storage::Vector< SasScatteringPoint> > reduced_data_container;
      storage::VectorND< 2, float> chi_score;
      storage::Vector< storage::VectorND< 2, float> > chi_container;

      std::string opencl_parameters
      (
        "OpenCLSaxsDebye(consider loops=0, analytic=0, excluded volume=1.0, hydration shell=0.0)"
      );
      std::string parameters
      (
        "SasDebye(consider loops=0, analytic=0, excluded volume=1.0, hydration shell=0.0)"
      );

      // try to create an opencl instance for this class, if that fails, use the non-opencl version
      std::stringstream err_stream;
      util::Implementation< SasDebyeInterface> saxs( opencl_parameters, err_stream);

      if( !saxs.IsDefined())
      {
        // Use the non-opencl version
        saxs = parameters;
      }
      saxs->SetExperimentalData( ORIGINAL_DATA);

      for( size_t iteration( 0); iteration < NUMBER_OF_ITERATIONS; ++iteration)
      {
        // Always select the first point
        subset.PushBack( ORIGINAL_DATA->GetScatteringLocation( 0));

        // randomly select points from each bin
        for( size_t i( 0); i < number_of_bins; ++i)
        {
          int index( SelectIndex( left_index, right_index));
          subset.PushBack( ORIGINAL_DATA->GetScatteringLocation( index));

          if( i != number_of_bins - 1)
          {
            left_index = right_index + 1;
            right_index = index_selected_points( i + 1);
          }
        } // end inner for loop

        reduced_data_container.PushBack( subset);
        subset.Reset();

        // 2) compute chi^2 for K iterations of shannon points
        storage::Vector< SasScatteringPoint> &data_point( reduced_data_container( iteration));

        // use the Clone to ShPtr to create a shared pointer to experimental data
        util::ShPtr< storage::Vector< SasScatteringPoint> > sp_reduced_data( util::CloneToShPtr( data_point));
        saxs->SetReducedExperimentalData( sp_reduced_data);

        // Compute the raw saxs profile
        SasExperimentalAndCalculatedData data_sets( saxs->operator()( PROTEIN_MODEL));

        // transform the profile into the desired form ( default of derivative score)
        SasExperimentalAndCalculatedData transformed_data( SasTransformation()( data_sets));

        // calculate the score
        float chi( score::SasType()( transformed_data));

        chi_score( 0) = chi;
        chi_score( 1) = iteration;

        chi_container.PushBack( chi_score);
        BCL_Message( util::Message::e_Standard, " iteration: " + util::Format()( iteration));

        // Reset indices for next iteration
        left_index = 0;
        right_index = index_selected_points( 0);
      } // end outer for loop

      // 3) chi^2 is taken as the median over k sampling rounds
      std::sort( chi_container.Begin(), chi_container.End());

      // compute the median and subtract 1 because of zero based indexing
      int median( ( NUMBER_OF_ITERATIONS / 2) - 1);

      int median_index( ( int)chi_container( median).Second());

      storage::Vector< SasScatteringPoint> &final_data( reduced_data_container( median_index));

      // Create a pointer to a new SasScatteringData Object with the reduced values

      util::ShPtr< SasScatteringData> sp_reduced_data
      (
        new SasScatteringData( final_data)
      );

      return sp_reduced_data;
    } // End Function

    //! @brief splits scatting profile into n (shannon) bins.  Selects datapoint to represent the curve from
    //!        each bin that has the least error associated with it.
    //! @param ORIGINAL_DATA - pointer to the preprocessed experimental data
    //! @param DENSITY_DATA - pointer to the transformed experimental data to P(r) domain
    util::ShPtr< SasScatteringData> SaxsDataReduction::SasSignalRecoveryEstimate
    (
      const util::ShPtr< SasScatteringData> ORIGINAL_DATA,
      const double DMAX
    )
    {
      // Tasks
      // 1) partition a SAS dataset into ns equal bins for a given dmax

      // compute the number of bins on the saxs interval
      size_t number_of_bins
      (
        ComputeShannonNumber
        (
          DMAX,
          SasAnalysis::ComputeQmax( *ORIGINAL_DATA)
        )
      );

      // compute the size of the bin on the interval

      double first_point( ORIGINAL_DATA->GetScatteringData().FirstElement().GetQvalue());
      double last_point( ORIGINAL_DATA->GetScatteringData().LastElement().GetQvalue());

      // compute the size of the bin on the interval based on the first and last point
      double binsize( ( last_point - first_point) / number_of_bins);
      double right_bin( binsize + first_point);
      int count( 0);

      storage::Vector< int> index_selected_points;

      // identify the index values that contain the desired ranges
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
        data_itr( ORIGINAL_DATA->Begin()),
        data_itr_end( ORIGINAL_DATA->End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        if( data_itr->GetQvalue() >= right_bin)
        {
          if( data_itr->GetQvalue() == right_bin)
          {
            index_selected_points.PushBack( count);
          }
          else
          {
            index_selected_points.PushBack( count - 1);
          }
          right_bin = right_bin + binsize;
        }
        ++count;
      }

      int left_index = 1;
      int right_index = index_selected_points( 0);

      // Set up containers to hold the data
      storage::Vector< SasScatteringPoint> estimated_data;

      // Always start with first point
      estimated_data.PushBack( ORIGINAL_DATA->GetScatteringLocation( 0));

      // select points with minimum error value from each bin
      for( size_t i( 0); i < index_selected_points.GetSize(); ++i)
      {
        // select points beginning with a left index of 1
        int index( SelectIndexMinError( ORIGINAL_DATA, left_index, right_index));
        estimated_data.PushBack( ORIGINAL_DATA->GetScatteringLocation( index));

        if( i != index_selected_points.GetSize() - 1)
        {
          left_index = right_index + 1;
          right_index = index_selected_points( i + 1);
        }
      } // end inner for loop

      util::ShPtr< SasScatteringData> sp_reduced_data
      (
        new SasScatteringData( estimated_data)
      );

      return sp_reduced_data;

    } // End Function

    //! @brief randomly select a point inside the given boundary
    //! @param LEFT - the left boundary of the interval
    //! @param RIGHT - the right boundary of the interval
    //! @return the randomly selected value inside the interval [ LEFT, RIGHT]
    int SaxsDataReduction::SelectIndex( int LEFT, int RIGHT)
    {
      return random::GetGlobalRandom().Random( LEFT, RIGHT);
    }

    //! @brief randomly select a point inside the given boundary
    //! @param ORIGINAL_DATA - pointer to the experimental data
    //! @param LEFT - the left boundary of the interval
    //! @param RIGHT - the right boundary of the interval
    //! @return the randomly selected value inside the interval [ LEFT, RIGHT]
    int SaxsDataReduction::SelectIndexMinError
    (
      const util::ShPtr< SasScatteringData> ORIGINAL_DATA,
      int LEFT,
      int RIGHT
    )
    {
      // Pick a point that minimizes the error in a bin
      double error( ORIGINAL_DATA->GetScatteringLocation( LEFT).GetError());

      int index( LEFT);

      for( int i( LEFT); i <= ( RIGHT); ++i)
      {
        if( ORIGINAL_DATA->GetScatteringLocation( i).GetError() < error)
        {
          error = ORIGINAL_DATA->GetScatteringLocation( i).GetError();
          index = i;
        }
        else
        {
          // Do nothing
        }
      }

      return index;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_saxs_opti_result.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SaxsOptiResult::s_Instance
    (
      GetObjectInstances().AddInstance( new SaxsOptiResult())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SaxsOptiResult::SaxsOptiResult() :
      m_C1( util::GetUndefined< float>()),
      m_C2( util::GetUndefined< float>()),
      m_Chi( util::GetUndefined< float>())
    {
    }

    //! @brief constructor from given input data
    //! @param QVALUE the scattering angle from SAXS curve ( x-axis)
    //! @param INTENSITY the wavelength intensity for a given scattering angle ( y-axis)
    //! @param ERROR - the variation in measurement for a given intensity
    SaxsOptiResult::SaxsOptiResult
    (
      const float C1,
      const float C2,
      const float CHI
    ) :
       m_C1( C1),
       m_C2( C2),
       m_Chi( CHI)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SaxsOptiResult
    SaxsOptiResult *SaxsOptiResult::Clone() const
    {
      return new SaxsOptiResult( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SaxsOptiResult::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SaxsOptiResult::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_C1, ISTREAM);
      io::Serialize::Read( m_C2, ISTREAM);
      io::Serialize::Read( m_Chi, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SaxsOptiResult::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_C1, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_C2, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Chi, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
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
#include "restraint/bcl_restraint_xlink_data.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_score_weight_set.h"
#include "score/bcl_score_restraint_xlink.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> XlinkData::s_Instance
    (
      util::Enumerated< Interface>::AddInstance( new XlinkData())
    );

    //! score for evaluating the agreement of a protein model with cross-linking restraints
    fold::Score XlinkData::e_ScoreXlinkRestraint( fold::GetScores().e_Undefined);

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    XlinkData::XlinkData() :
      m_Handler( GetDefaultExtension())
    {
    }

    //! @brief returns a pointer to a new XlinkData
    //! @return pointer to a new XlinkData
    XlinkData *XlinkData::Clone() const
    {
      return new XlinkData( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &XlinkData::GetAlias() const
    {
      static const std::string s_name( "Xlink");
      return s_name;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &XlinkData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the default file extension of files containing cross-link restraints
    //! @return the default file extension of files containing cross-link restraints
    const std::string &XlinkData::GetDefaultExtension() const
    {
      static const std::string s_extension( ".xlink_bcl");
      return s_extension;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initialize the scores and add them to Scores enumerator
    void XlinkData::InitializeScores()
    {
      // read in restraints
      if( !m_Restraints.IsDefined() && m_Handler.Exists())
      {
        m_Restraints = util::CloneToShPtr( m_Handler.ReadRestraintsFromFile());
      }
      if( !e_ScoreXlinkRestraint.IsDefined())
      {
        e_ScoreXlinkRestraint = fold::GetScores().AddScore
        (
          util::ShPtr< score::ProteinModel>( new score::RestraintXlink( m_Restraints))
        );
      }
    }

    //! @brief sets the weights of scores in a weight set
    //! @param SCORE_WEIGHT_SET the score weight set that will be modified
    void XlinkData::ModifyScoreWeightSet( fold::ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
      SCORE_WEIGHT_SET.SetWeight( e_ScoreXlinkRestraint, 10);
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void XlinkData::MergeAndModifyMutateTree( fold::MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< fold::MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }
 
    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer XlinkData::GetSerializer() const
    {
      io::Serializer serial;
      serial.SetClassDescription( "Cross linking restraints");
      serial.AddInitializer
      (
        "",
        "Handler for reading cross-link data",
        io::Serialization::GetAgent( &m_Handler),
        HandlerAtomDistanceAssigned( GetDefaultExtension()).GetString()
      );
      return serial;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads members from an input stream
    //! @param ISTREAM input stream to read members from
    //! @return returns the input stream
    std::istream &XlinkData::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Restraints, ISTREAM);

      return ISTREAM;
    }

    //! @brief writes members into an output stream
    //! @param OSTREAM output stream to write members into
    //! @INDENT number of indentations to use
    //! @return returns the output stream
    std::ostream &XlinkData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Restraints, OSTREAM, INDENT);

      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
