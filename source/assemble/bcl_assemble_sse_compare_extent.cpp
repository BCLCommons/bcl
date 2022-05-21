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
#include "assemble/bcl_assemble_sse_compare_extent.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> SSECompareExtent::s_Instance
    (
      GetObjectInstances().AddInstance( new SSECompareExtent())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SSECompareExtent::SSECompareExtent() :
      m_ExtentTolerance(),
      m_Axis()
    {
    }

    //! @brief constructor from tolerance and axis
    //! @param TOLERANCE extent tolerance to be used
    //! @param AXIS axis to be used (e.g. z axis)
    SSECompareExtent::SSECompareExtent( const double TOLERANCE, const coord::Axis &AXIS) :
      m_ExtentTolerance( TOLERANCE),
      m_Axis( AXIS)
    {
    }

    //! @brief virtual copy constructor
    SSECompareExtent *SSECompareExtent::Clone() const
    {
      return new SSECompareExtent( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that tests whether the two sses have the same extent within some tolerance
    //! @param SSE_A first sse to be tested
    //! @param SSE_B second sse to be tested
    //! @return bool whether the two sses agree in their extent within some tolerance (m_ExtentTolerance)
    bool SSECompareExtent::operator()( const SSE &SSE_A, const SSE &SSE_B) const
    {
      // if the the difference of the extent of the two sses along the specified axis is not greater than the tolerance
      // return true, otherwise false
      return ( math::Absolute( SSE_A.GetExtent( m_Axis) - SSE_B.GetExtent( m_Axis)) < m_ExtentTolerance);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSECompareExtent::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSECompareExtent::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
