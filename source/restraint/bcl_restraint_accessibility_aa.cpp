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
