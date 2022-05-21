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
#include "assemble/bcl_assemble_sse_geometry_packing_compare.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_geometry_packing.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief virtual copy constructor
    SSEGeometryPackingCompareInteractionWeight *SSEGeometryPackingCompareInteractionWeight::Clone() const
    {
      return new SSEGeometryPackingCompareInteractionWeight( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &SSEGeometryPackingCompareInteractionWeight::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief return true if SSE_GEOMETRY_PACKING_A has better interaction weight than SSE_GEOMETRY_PACKING_B, or is closer by distance
    //! @param SSE_GEOMETRY_PACKING_A first SSEGeometryPacking
    //! @param SSE_GEOMETRY_PACKING_B second SSEGeometryPacking
    //! @return true if SSE_GEOMETRY_PACKING_A has better interaction weight than SSE_GEOMETRY_PACKING_B, or smaller distance
    bool SSEGeometryPackingCompareInteractionWeight::operator()
    (
      const SSEGeometryPacking &SSE_GEOMETRY_PACKING_A,
      const SSEGeometryPacking &SSE_GEOMETRY_PACKING_B
    ) const
    {
      // compare the interaction weights
      if( SSE_GEOMETRY_PACKING_A.GetInteractionWeight() > SSE_GEOMETRY_PACKING_B.GetInteractionWeight())
      {
        return true;
      }
      // for equal interaction weight, consider the closest by distance
      else if( SSE_GEOMETRY_PACKING_A.GetInteractionWeight() == SSE_GEOMETRY_PACKING_B.GetInteractionWeight())
      {
        return SSE_GEOMETRY_PACKING_A.GetDistance() < SSE_GEOMETRY_PACKING_B.GetDistance();
      }

      // end
      return false;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEGeometryPackingCompareInteractionWeight::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSEGeometryPackingCompareInteractionWeight::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief virtual copy constructor
    SSEGeometryPackingCompareDistance *SSEGeometryPackingCompareDistance::Clone() const
    {
      return new SSEGeometryPackingCompareDistance( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &SSEGeometryPackingCompareDistance::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief return true if SSE_GEOMETRY_PACKING_A has a shorter distance than SSE_GEOMETRY_PACKING_B
    //! @param SSE_GEOMETRY_PACKING_A first SSEGeometryPacking
    //! @param SSE_GEOMETRY_PACKING_B second SSEGeometryPacking
    //! @return true if SSE_GEOMETRY_PACKING_A has a shorter distance than SSE_GEOMETRY_PACKING_B
    bool SSEGeometryPackingCompareDistance::operator()
    (
      const SSEGeometryPacking &SSE_GEOMETRY_PACKING_A,
      const SSEGeometryPacking &SSE_GEOMETRY_PACKING_B
    ) const
    {
      // compare the interaction weights
      return SSE_GEOMETRY_PACKING_A.GetDistance() > SSE_GEOMETRY_PACKING_B.GetDistance()
             ? false
             : SSE_GEOMETRY_PACKING_A.GetDistance() < SSE_GEOMETRY_PACKING_B.GetDistance()
               ? true
               : SSE_GEOMETRY_PACKING_A.GetInteractionWeight() > SSE_GEOMETRY_PACKING_B.GetInteractionWeight();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEGeometryPackingCompareDistance::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSEGeometryPackingCompareDistance::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
