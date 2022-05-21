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
#include "score/bcl_score_body_assignment.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "io/bcl_io_serialization.h"
#include "restraint/bcl_restraint_assignment.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> BodyAssignment::s_Instance
    (
      GetObjectInstances().AddInstance( new BodyAssignment())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    BodyAssignment::BodyAssignment() :
      m_BodyAgreement()
    {
    }

    //! @brief construct from a ShPtr to FunctionInterface defining the measure of coord::GeometryInterface agreement
    //! @param BODY_AGREEMENT ShPtr to FunctionInterface defining the measure of coord::GeometryInterface agreement
    BodyAssignment::BodyAssignment
    (
      const util::ShPtr< math::BinaryFunctionInterfaceSerializable< assemble::SSEGeometryInterface, assemble::SSE, double> > &BODY_AGREEMENT
    ) :
      m_BodyAgreement( *BODY_AGREEMENT)
    {
    }

    //! @brief construct from a FunctionInterface defining the measure of coord::GeometryInterface agreement
    //! @param BODY_AGREEMENT FunctionInterface defining the measure of coord::GeometryInterface agreement
    BodyAssignment::BodyAssignment
    (
      const math::BinaryFunctionInterfaceSerializable< assemble::SSEGeometryInterface, assemble::SSE, double> &BODY_AGREEMENT
    ) :
      m_BodyAgreement( BODY_AGREEMENT.Clone())
    {
    }

    //! @brief Clone is the virtual copy constructor
    BodyAssignment *BodyAssignment::Clone() const
    {
      return new BodyAssignment( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &BodyAssignment::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get scheme
    //! @return scheme
    const std::string &BodyAssignment::GetScheme() const
    {
      return m_BodyAgreement->GetScheme();
    }

    //! @brief get scheme
    //! @return scheme
    const std::string &BodyAssignment::GetAlias() const
    {
      return m_BodyAgreement.GetAlias();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief operator() which takes an Assignment for calculating its score
    //! @param ASSIGNMENT the assignment which contains the coord::Bodies who will be scored
    //! @return return a double which is the score of the Assignment
    double BodyAssignment::operator()
    (
      const restraint::SSEAssignment &ASSIGNMENT
    ) const
    {
      // create double "score" for holding the score of "ASSIGNMENT"
      double score( 0);

      // iterate through the assigned Bodies in the Assignment
      for
      (
        restraint::GroupCollection< size_t, assemble::SSE>::const_iterator
          itr( ASSIGNMENT.GetGroupCollection().Begin()), itr_end( ASSIGNMENT.GetGroupCollection().End());
        itr != itr_end;
        ++itr
      )
      {
        // add the agreement of the two bodies indicated by "itr" to "score"
        const double agreement
        (
          m_BodyAgreement->operator()
          (
            *ASSIGNMENT.GetRestraint()->operator()( itr->first),
            *itr->second.FirstElement()
          )
        );
        BCL_MessageDbg
        (
          "score::BodyAssignment::operator() current assignment body agreement "
          + util::Format()( agreement)
        );
        score += agreement;
      }

      BCL_MessageDbg
      (
        "score::BodyAssignment::operator() assignment score "
        + util::Format()( score)
      );

      // return "score" which is the score of the assignment
      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer BodyAssignment::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        "computes body agreement for sses"
      );
      serializer.AddInitializer
      (
        "",
        "method used to compute body agreement",
        io::Serialization::GetAgent( &m_BodyAgreement)
      );
      return serializer;
    }

  } // namespace score
} // namespace bcl
