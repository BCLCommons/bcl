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
#include "assemble/bcl_assemble_sse_geometry_packing_criteria.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_geometry_packing.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSEGeometryPackingCriteriaCombine::s_Instance
    (
      GetObjectInstances().AddInstance( new SSEGeometryPackingCriteriaCombine())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SSEGeometryPackingCriteriaCombine::SSEGeometryPackingCriteriaCombine() :
      m_CriteriaVector()
    {
    }

    //! @brief constructor from a criteria vector
    //! @param CRITERIA_VECTOR vector of criteria
    SSEGeometryPackingCriteriaCombine::SSEGeometryPackingCriteriaCombine
    (
      const util::ShPtrVector< math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> > &CRITERIA_VECTOR
    ) :
      m_CriteriaVector( CRITERIA_VECTOR)
    {
    }

    //! @brief virtual copy constructor
    SSEGeometryPackingCriteriaCombine *SSEGeometryPackingCriteriaCombine::Clone() const
    {
      return new SSEGeometryPackingCriteriaCombine( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &SSEGeometryPackingCriteriaCombine::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief evaluate the given SSEGeometryPacking according to the group of criteria
    //! @param SSE_GEOMETRY_PACKING SSEGeometryPacking to be evaluated
    //! @return true if given SSEGeometryPacking follows all criteria
    bool SSEGeometryPackingCriteriaCombine::operator()
    (
      const SSEGeometryPacking &SSE_GEOMETRY_PACKING
    ) const
    {
      // iterate over criteria
      for
      (
        util::ShPtrVector< math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> >::const_iterator
          criteria_itr( m_CriteriaVector.Begin()), criteria_itr_end( m_CriteriaVector.End());
        criteria_itr != criteria_itr_end; ++criteria_itr
      )
      {
        // return false if the criteria is not correct
        if( !( *criteria_itr)->operator ()( SSE_GEOMETRY_PACKING))
        {
          return false;
        }
      }

      // return true if this point is reached, meaning all the criteria were followed
      return true;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEGeometryPackingCriteriaCombine::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_CriteriaVector, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSEGeometryPackingCriteriaCombine::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_CriteriaVector, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSEGeometryPackingCriteriaDistance::s_Instance
    (
      GetObjectInstances().AddInstance( new SSEGeometryPackingCriteriaDistance())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SSEGeometryPackingCriteriaDistance::SSEGeometryPackingCriteriaDistance() :
      m_DistanceCutoff(),
      m_Comparison()
    {
    }

    //! @brief constructor from a distance cutoff and a comparison
    //! @param DISTANCE_CUTOFF distance cutoff
    //! @param COMPARISON Comparison operator
    SSEGeometryPackingCriteriaDistance::SSEGeometryPackingCriteriaDistance
    (
      const double DISTANCE_CUTOFF,
      const math::Comparisons< double>::Comparison &COMPARISON
    ) :
      m_DistanceCutoff( DISTANCE_CUTOFF),
      m_Comparison( COMPARISON)
    {
    }

    //! @brief virtual copy constructor
    SSEGeometryPackingCriteriaDistance *SSEGeometryPackingCriteriaDistance::Clone() const
    {
      return new SSEGeometryPackingCriteriaDistance( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &SSEGeometryPackingCriteriaDistance::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief evaluate the given SSEGeometryPacking according to distance cutoff criteria
    //! @param SSE_GEOMETRY_PACKING SSEGeometryPacking to be evaluated
    //! @return true if SSEGeometryPacking follows distance cutoff criteria
    bool SSEGeometryPackingCriteriaDistance::operator()( const SSEGeometryPacking &SSE_GEOMETRY_PACKING) const
    {
      // evaluate the distance cutoff using the comparison operator
      return ( *m_Comparison)->operator ()( SSE_GEOMETRY_PACKING.GetDistance(), m_DistanceCutoff);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEGeometryPackingCriteriaDistance::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_DistanceCutoff, ISTREAM);
      io::Serialize::Read( m_Comparison, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSEGeometryPackingCriteriaDistance::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_DistanceCutoff, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Comparison, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSEGeometryPackingCriteriaContactType::s_Instance
    (
      GetObjectInstances().AddInstance( new SSEGeometryPackingCriteriaContactType())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SSEGeometryPackingCriteriaContactType::SSEGeometryPackingCriteriaContactType() :
      m_ContactTypes()
    {
    }

    //! @brief constructor from a contact type
    //! @param CONTACT_TYPE requested contact type
    SSEGeometryPackingCriteriaContactType::SSEGeometryPackingCriteriaContactType
    (
      const contact::Type &CONTACT_TYPE
    ) :
      m_ContactTypes( CONTACT_TYPE)
    {
    }

    //! @brief constructor from a set of contact types
    //! @param CONTACT_TYPES set of requested contact types
    SSEGeometryPackingCriteriaContactType::SSEGeometryPackingCriteriaContactType
    (
      const storage::Set< contact::Type> &CONTACT_TYPES
    ) :
      m_ContactTypes( CONTACT_TYPES)
    {
    }

    //! @brief virtual copy constructor
    SSEGeometryPackingCriteriaContactType *SSEGeometryPackingCriteriaContactType::Clone() const
    {
      return new SSEGeometryPackingCriteriaContactType( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &SSEGeometryPackingCriteriaContactType::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief evaluate the given SSEGeometryPacking according to contact type criteria
    //! @param SSE_GEOMETRY_PACKING SSEGeometryPacking to be evaluated
    //! @return true if SSEGeometryPacking follows contact type criteria
    bool SSEGeometryPackingCriteriaContactType::operator()
    (
      const SSEGeometryPacking &SSE_GEOMETRY_PACKING
    ) const
    {
      // search for the corresponding contact type
      return m_ContactTypes.Contains( SSE_GEOMETRY_PACKING.GetContactType());
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEGeometryPackingCriteriaContactType::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ContactTypes, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSEGeometryPackingCriteriaContactType::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ContactTypes, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSEGeometryPackingCriteriaInteractionWeight::s_Instance
    (
      GetObjectInstances().AddInstance( new SSEGeometryPackingCriteriaInteractionWeight())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SSEGeometryPackingCriteriaInteractionWeight::SSEGeometryPackingCriteriaInteractionWeight() :
      m_InteractionWeightCutoff(),
      m_Comparison()
    {
    }

    //! @brief constructor from a interaction weight cutoff and a comparison
    //! @param INTERACTION_WEIGHT_CUTOFF interaction weight cutoff
    //! @param COMPARISON Comparison operator
    SSEGeometryPackingCriteriaInteractionWeight::SSEGeometryPackingCriteriaInteractionWeight
    (
      const double INTERACTION_WEIGHT_CUTOFF,
      const math::Comparisons< double>::Comparison &COMPARISON
    ) :
      m_InteractionWeightCutoff( INTERACTION_WEIGHT_CUTOFF),
      m_Comparison( COMPARISON)
    {
    }

    //! @brief virtual copy constructor
    SSEGeometryPackingCriteriaInteractionWeight *SSEGeometryPackingCriteriaInteractionWeight::Clone() const
    {
      return new SSEGeometryPackingCriteriaInteractionWeight( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &SSEGeometryPackingCriteriaInteractionWeight::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief evaluate the given SSEGeometryPacking according to interaction weight cutoff criteria
    //! @param SSE_GEOMETRY_PACKING SSEGeometryPacking to be evaluated
    //! @return true if SSEGeometryPacking follows interaction weight cutoff criteria
    bool SSEGeometryPackingCriteriaInteractionWeight::operator()
    (
      const SSEGeometryPacking &SSE_GEOMETRY_PACKING
    ) const
    {
      // evaluate the interaction weight cutoff using the comparison operator
      return ( *m_Comparison)->operator ()( SSE_GEOMETRY_PACKING.GetInteractionWeight(), m_InteractionWeightCutoff);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEGeometryPackingCriteriaInteractionWeight::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_InteractionWeightCutoff, ISTREAM);
      io::Serialize::Read( m_Comparison, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSEGeometryPackingCriteriaInteractionWeight::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_InteractionWeightCutoff, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Comparison, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a strand weight cutoff and a comparison
    //! @param STRAND_WEIGHT_CUTOFF strand weight cutoff
    //! @param COMPARISON Comparison operator
    SSEGeometryPackingCriteriaStrandWeight::SSEGeometryPackingCriteriaStrandWeight
    (
      const double STRAND_WEIGHT_CUTOFF,
      const math::Comparisons< double>::Comparison &COMPARISON
    ) :
      m_StrandWeightCutoff( STRAND_WEIGHT_CUTOFF),
      m_Comparison( COMPARISON)
    {
    }

    //! @brief virtual copy constructor
    SSEGeometryPackingCriteriaStrandWeight *SSEGeometryPackingCriteriaStrandWeight::Clone() const
    {
      return new SSEGeometryPackingCriteriaStrandWeight( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &SSEGeometryPackingCriteriaStrandWeight::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief evaluate the given SSEGeometryPacking according to strand weight cutoff criteria
    //! @param SSE_GEOMETRY_PACKING SSEGeometryPacking to be evaluated
    //! @return true if SSEGeometryPacking follows strand weight cutoff criteria
    bool SSEGeometryPackingCriteriaStrandWeight::operator()
    (
      const SSEGeometryPacking &SSE_GEOMETRY_PACKING
    ) const
    {
      // evaluate the weight cutoff using the comparison operator
      return ( *m_Comparison)->operator ()( SSE_GEOMETRY_PACKING.GetStrandStrandPairingWeight(), m_StrandWeightCutoff);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEGeometryPackingCriteriaStrandWeight::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_StrandWeightCutoff, ISTREAM);
      io::Serialize::Read( m_Comparison, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSEGeometryPackingCriteriaStrandWeight::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_StrandWeightCutoff, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Comparison, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSEGeometryPackingCriteriaDistancePerType::s_Instance
    (
      GetObjectInstances().AddInstance( new SSEGeometryPackingCriteriaDistancePerType())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SSEGeometryPackingCriteriaDistancePerType::SSEGeometryPackingCriteriaDistancePerType()
    {
    }

    //! @brief virtual copy constructor
    SSEGeometryPackingCriteriaDistancePerType *SSEGeometryPackingCriteriaDistancePerType::Clone() const
    {
      return new SSEGeometryPackingCriteriaDistancePerType( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &SSEGeometryPackingCriteriaDistancePerType::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief evaluate the given SSEGeometryPacking according to the distance criteria
    //! @param SSE_GEOMETRY_PACKING SSEGeometryPacking to be evaluated
    //! @return true if SSEGeometryPacking follows the distance criteria
    bool SSEGeometryPackingCriteriaDistancePerType::operator()
    (
      const SSEGeometryPacking &SSE_GEOMETRY_PACKING
    ) const
    {
      // return true if the distance is within the range given by the contact type
      return SSE_GEOMETRY_PACKING.GetContactType()->GetDistanceRange().IsWithin( SSE_GEOMETRY_PACKING.GetDistance());
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEGeometryPackingCriteriaDistancePerType::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSEGeometryPackingCriteriaDistancePerType::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
