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

#ifndef BCL_SCORE_BODY_EXTENT_AGREEMENT_H_
#define BCL_SCORE_BODY_EXTENT_AGREEMENT_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "coord/bcl_coord.fwd.hh"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_binary_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BodyExtentAgreement
    //! @brief scores the agreement of two coord::Bodies based on the deviation between
    //! the x, y, or z extents.
    //!
    //! @see @link example_score_body_extent_agreement.cpp @endlink
    //! @author alexanns
    //! @date 07/14/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API BodyExtentAgreement :
      public math::BinaryFunctionInterfaceSerializable< assemble::SSEGeometryInterface, assemble::SSE, double>
    {

    private:

    //////////
    // data //
    //////////

      double      m_LowerTolerance;
      double      m_LowerTransitionWidth;
      double      m_UpperTolerance;
      double      m_UpperTransitionWidth;
      double      m_EnergyWellDepth;
      coord::Axis m_Axis;

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
      BodyExtentAgreement();

      //! @brief constructor
      //! @param LOWER_TOLERANCE the amount an extent can be shorter than the restraint and still be considered perfect
      //! @param LOWER_TRANSITION_WIDTH width of range over which a too short extent goes from perfect to non-agreeing
      //! @param UPPER_TOLERANCE the amount an extent can be longer than the restraint and still be considered perfect
      //! @param UPPER_TRANSITION_WIDTH width of range over which a too long extent goes from perfect to non-agreeing
      //! @param ENERGY_WELL_DEPTH the amount that an agreeing restraint will give a bonus to the energy
      //! @param AXIS defines which extent (x, y, z) is being checked for restraint agreement
      BodyExtentAgreement
      (
        const double LOWER_TOLERANCE, const double LOWER_TRANSITION_WIDTH, const double UPPER_TOLERANCE,
        const double UPPER_TRANSITION_WIDTH, const double ENERGY_WELL_DEPTH, const coord::Axis &AXIS
      );

      //! @brief Clone is the virtual copy constructor
      BodyExtentAgreement *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the scheme
      //! @return string that describes that score
      const std::string &GetScheme() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
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

      //! @brief GetExtentDeviation gives the value of the extent position difference for an axis
      //! @param BODY_A is the first coord::GeometryInterface used in the calculation
      //! @param BODY_B is the second coord::GeometryInterface used in the calculation
      //! @return returns the value of the difference in position of the extent of an axis
      double GetExtentDeviation
      (
        const coord::GeometryInterface &BODY_A, const coord::GeometryInterface &BODY_B
      ) const;

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

    }; // class BodyExtentAgreement

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_BODY_EXTENT_AGREEMENT_H_ 
