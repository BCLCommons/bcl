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

#ifndef BCL_ASSEMBLE_SSE_GEOMETRY_COMPARE_H_
#define BCL_ASSEMBLE_SSE_GEOMETRY_COMPARE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEGeometryWithinSizeTolerance
    //! @brief returns whether an SSE and an SSEGeometryPhiPsi are within a specified size tolerance
    //! @details Compares SSE and SSEGeometryPhiPsi objects and returns whether they are within +/- one tolerance
    //!          length.  This is used for deciding if an SSE should be fit into a given geometry from a fold template.
    //!
    //! @see @link example_assemble_sse_geometry_compare.cpp @endlink
    //! @author weinerbe
    //! @date Nov 11, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEGeometryWithinSizeTolerance :
      public math::BinaryFunctionInterfaceSerializable< SSE, SSEGeometryPhiPsi, bool>
    {

    private:

    //////////
    // data //
    //////////

      //! size tolerance allowed for each sstype
      storage::Map< biol::SSType, size_t> m_SizeTolerances;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from helix and strand length tolerances
      //! @param HELIX_TOLERANCE helix tolerance
      //! @param STRAND_TOLERANCE strand tolerance
      SSEGeometryWithinSizeTolerance( const size_t HELIX_TOLERANCE = 2, const size_t STRAND_TOLERANCE = 2);

      //! @brief Clone function
      //! @return pointer to new SSEGeometryWithinSizeTolerance
      SSEGeometryWithinSizeTolerance *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief returns whether two SSE geometry interface objects are within a fragment length of each other
      //! @param SSE_TO_COMPARE SSE to be compared
      //! @param GEOMTRY_PHI_PSI SSEGeometryPhiPsi to be compared
      //! @return bool whether two SSE geometry interface objects are within a fragment length of each other
      bool operator()( const SSE &SSE_TO_COMPARE, const SSEGeometryPhiPsi &GEOMTRY_PHI_PSI) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class SSEGeometryWithinSizeTolerance

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_SSE_GEOMETRY_COMPARE_H_ 
