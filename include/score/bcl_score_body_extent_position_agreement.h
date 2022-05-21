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

#ifndef BCL_SCORE_BODY_EXTENT_POSITION_AGREEMENT_H_
#define BCL_SCORE_BODY_EXTENT_POSITION_AGREEMENT_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "coord/bcl_coord.fwd.hh"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_binary_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BodyExtentPositionAgreement
    //! @brief scores the agreement of two coord::Bodies based on the deviation between
    //! the positions of their x, y, and z extents, where the position is given by (origin + extent).
    //!
    //! @see @link example_score_body_extent_position_agreement.cpp @endlink
    //! @author alexanns
    //! @date 07/05/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API BodyExtentPositionAgreement :
      public math::BinaryFunctionInterfaceSerializable< assemble::SSEGeometryInterface, assemble::SSE, double>
    {

    private:

    //////////
    // data //
    //////////

      double m_LowerTolerance;
      double m_LowerTransitionWidth;
      double m_UpperTolerance;
      double m_UpperTransitionWidth;
      double m_EnergyWellDepth;

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
      BodyExtentPositionAgreement();

      //! @brief Clone is the virtual copy constructor
      BodyExtentPositionAgreement *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the name of the object
      //! @return the name of the object
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief operator() which takes a density rod and an assigned SSE
      //! @param BODY the density rod
      //! @param SSE the secondary structure element assigned to that body
      //! @return return a double which is the agreement the body with the sse
      double operator()( const assemble::SSEGeometryInterface &BODY, const assemble::SSE &SSE) const;

      //! @brief GetAbsoluteExtentDeviation gives the absolute value of the extent position difference for an axis
      //! @param BODY_A is the first coord::GeometryInterface used in the calculation
      //! @param BODY_B is the second coord::GeometryInterface used in the calculation
      //! @param AXIS is the RotationAxis for which the calculation will be performed
      //! @return returns the absolute value of the difference in position of the extent of an axis
      double GetAbsoluteExtentDeviation
      (
        const coord::GeometryInterface &BODY_A, const coord::GeometryInterface &BODY_B, const coord::Axis &AXIS
      ) const;

      //! @brief GetExtentPosition calculates the position of an extent for a defined axis for a coord::GeometryInterface
      //! @param BODY the coord::GeometryInterface for which the position of its extent in a given axis will be calculated
      //! @param AXIS is the RotationAxis for which the extent position will be calculated
      //! @return double which is the value of the origin in the axis "XYZ" plus one half the extent of BODY in that
      //!         axis direction
      double GetExtentPosition( const coord::GeometryInterface &BODY, const coord::Axis &AXIS) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief write to ostream
      //! @param OSTREAM is the output stream
      //! @return returns the output stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! @brief read from istream
      //! @param ISTREAM is the input stream
      //! @return returns the input stream
      std::istream &Read( std::istream &ISTREAM);

    }; // class BodyExtentPositionAgreement

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_BODY_EXTENT_POSITION_AGREEMENT_H_ 
