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
