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

#ifndef BCL_RESTRAINT_CONTAINS_BODY_ORIGIN_H_
#define BCL_RESTRAINT_CONTAINS_BODY_ORIGIN_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "linal/bcl_linal.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_binary_function_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ContainsBodyOrigin
    //! @brief ContainsBodyOrigin is used to determine if a coord::GeometryInterface contains the origin of another coord::GeometryInterface
    //!        This is for use to determine if a body restraint is occupied or not.
    //!
    //! @see @link example_restraint_contains_body_origin.cpp @endlink
    //! @author alexanns, linders
    //! @date 06/30/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ContainsBodyOrigin :
      public util::BinaryFunctionInterface< assemble::SSEGeometryInterface, assemble::SSE, bool>
    {

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
      ContainsBodyOrigin();

      //! @brief virtual copy constructor
      ContainsBodyOrigin *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator() determines if the BODY is occupied by any of the geometries of the SSE
      //! @details this is done by checking if any of the fragments in the body will contain the center of any of the
      //!          fragments in the sse
      //! @param BODY representing the density rod
      //! @param SSE represents the secondary structure element which contains geometries/fragments
      //! @return returns true if BODY contains any of the SSE's geometries
      bool operator()( const assemble::SSEGeometryInterface &BODY, const assemble::SSE &SSE) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read distance from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write distance to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief ContainsBodyOriginCheck determines if the coord::GeometryInterface envelops the origin
      //! @param BODY coord::GeometryInterface which is the body involved
      //! @param ORIGIN origin for which is checked whether it is within body
      //! @return returns true if Body envelops the origin
      bool ContainsBodyOriginCheck( const assemble::SSEGeometryInterface &BODY, const linal::Vector3D &ORIGIN) const;

    }; //end class ContainsBodyOrigin

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_CONTAINS_BODY_ORIGIN_H_
