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

#ifndef BCL_SCORE_BODY_ASSIGNMENT_H_
#define BCL_SCORE_BODY_ASSIGNMENT_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "restraint/bcl_restraint.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_geometry_interface.h"
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BodyAssignment
    //! @brief is for scoring an restraint::Assignment which has coord::Bodies as restraints assigned
    //! with coord::Bodies
    //!
    //! @see @link example_score_body_assignment.cpp @endlink
    //! @author alexanns
    //! @date 07/05/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API BodyAssignment :
      public math::FunctionInterfaceSerializable< restraint::SSEAssignment, double>
    {
    private:

    //////////
    // data //
    //////////

      //! ShPtr to a FunctionInterface "m_BodyAgreement" defines how the agreement of two Bodies is calculated
      util::Implementation
      <
        math::BinaryFunctionInterfaceSerializable< assemble::SSEGeometryInterface, assemble::SSE, double>
      > m_BodyAgreement;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      BodyAssignment();

      //! @brief construct from a ShPtr to FunctionInterface defining the measure of coord::GeometryInterface agreement
      //! @param BODY_AGREEMENT ShPtr to FunctionInterface defining the measure of coord::GeometryInterface agreement
      BodyAssignment
      (
        const util::ShPtr< math::BinaryFunctionInterfaceSerializable< assemble::SSEGeometryInterface, assemble::SSE, double> > &BODY_AGREEMENT
      );

      //! @brief construct from a FunctionInterface defining the measure of coord::GeometryInterface agreement
      //! @param BODY_AGREEMENT FunctionInterface defining the measure of coord::GeometryInterface agreement
      BodyAssignment
      (
        const math::BinaryFunctionInterfaceSerializable< assemble::SSEGeometryInterface, assemble::SSE, double> &BODY_AGREEMENT
      );

      //! @brief Clone is the virtual copy constructor
      BodyAssignment *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get scheme
      //! @return scheme
      const std::string &GetScheme() const;

      //! @brief get scheme
      //! @return scheme
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief operator() which takes an Assignment for calculating its score
      //! @param ASSIGNMENT the assignment which contains the coord::Bodies who will be scored
      //! @return return a double which is the score of the Assignment
      double operator()
      (
        const restraint::SSEAssignment &ASSIGNMENT
      ) const;

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class BodyAssignment

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_BODY_ASSIGNMENT_H_
