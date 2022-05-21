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
#include "score/bcl_score_aa_assignment_property.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AAAssignmentProperty::s_Instance
    (
      GetObjectInstances().AddInstance( new AAAssignmentProperty())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor
    //! @param PROPERTY property that will be compared
    AAAssignmentProperty::AAAssignmentProperty( const biol::AATypeData::PropertyType PROPERTY) :
      m_Property( PROPERTY)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new AAAssignmentProperty copied from this one
    AAAssignmentProperty *AAAssignmentProperty::Clone() const
    {
      return new AAAssignmentProperty( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AAAssignmentProperty::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief operator that calculates the score between two assigned members
    //! @param MEMBER_A amino acid A that is compared
    //! @param MEMBER_B amino acid A that is compared
    //! @return difference of property values of the two aminoacids < 0.0
    double AAAssignmentProperty::operator()( const biol::AABase &MEMBER_A, const biol::AABase &MEMBER_B) const
    {
      return -math::Absolute( MEMBER_A.GetType()->GetAAProperty( m_Property) - MEMBER_B.GetType()->GetAAProperty( m_Property));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @param INDENT indentation
    //! @return returns the output stream
    std::ostream &AAAssignmentProperty::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Property, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &AAAssignmentProperty::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Property, ISTREAM);

      // end
      return ISTREAM;
    }

  } // namespace score
} // namespace bcl
