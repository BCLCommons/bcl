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

#ifndef BCL_ASSEMBLE_SSE_GEOMETRY_PACKING_COMPARE_H_
#define BCL_ASSEMBLE_SSE_GEOMETRY_PACKING_COMPARE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "math/bcl_math_binary_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEGeometryPackingCompareInteractionWeight
    //! @brief class compares interaction weight of SSEGeometryPacking pair, on equality the distance
    //! @details This binary function interface derived class has an operator that returns true if first one
    //! has a better interaction weight than the second one, if they are equal, it returns true if the first one has a
    //! smaller distance.
    //!
    //! @see @link example_assemble_sse_geometry_packing_compare.cpp @endlink
    //! @author karakam
    //! @date Jul 27, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEGeometryPackingCompareInteractionWeight :
      public math::BinaryFunctionInterfaceSerializable< SSEGeometryPacking, SSEGeometryPacking, bool>
    {
    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief virtual copy constructor
      SSEGeometryPackingCompareInteractionWeight *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief return true if SSE_GEOMETRY_PACKING_A has better interaction weight than SSE_GEOMETRY_PACKING_B, or is closer by distance
      //! @param SSE_GEOMETRY_PACKING_A first SSEGeometryPacking
      //! @param SSE_GEOMETRY_PACKING_B second SSEGeometryPacking
      //! @return true if SSE_GEOMETRY_PACKING_A has better interaction weight than SSE_GEOMETRY_PACKING_B, or smaller distance
      bool operator()
      (
        const SSEGeometryPacking &SSE_GEOMETRY_PACKING_A,
        const SSEGeometryPacking &SSE_GEOMETRY_PACKING_B
      ) const;

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

    }; // class SSEGeometryPackingCompareInteractionWeight

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEGeometryPackingCompareDistance
    //! @brief compares distance of SSEGeometryPacking pair
    //! @details This binary function interface derived class has an operator that returns true if first one
    //! has a smaller distance than the second one.
    //!
    //! @see @link example_assemble_sse_geometry_packing_compare.cpp @endlink
    //! @author weinerbe
    //! @date Aug 31, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEGeometryPackingCompareDistance :
      public math::BinaryFunctionInterfaceSerializable< SSEGeometryPacking, SSEGeometryPacking, bool>
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new SSEGeometryPackingCompareDistance
      SSEGeometryPackingCompareDistance *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief return true if SSE_GEOMETRY_PACKING_A has a shorter distance than SSE_GEOMETRY_PACKING_B
      //! @param SSE_GEOMETRY_PACKING_A first SSEGeometryPacking
      //! @param SSE_GEOMETRY_PACKING_B second SSEGeometryPacking
      //! @return true if SSE_GEOMETRY_PACKING_A has a shorter distance than SSE_GEOMETRY_PACKING_B
      bool operator()
      (
        const SSEGeometryPacking &SSE_GEOMETRY_PACKING_A,
        const SSEGeometryPacking &SSE_GEOMETRY_PACKING_B
      ) const;

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

    }; // class SSEGeometryPackingCompareDistance

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_SSE_GEOMETRY_PACKING_COMPARE_H_ 
