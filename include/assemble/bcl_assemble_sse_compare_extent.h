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

#ifndef BCL_ASSEMBLE_SSE_COMPARE_EXTENT_H_
#define BCL_ASSEMBLE_SSE_COMPARE_EXTENT_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_assemble_sse.h"
#include "util/bcl_util_binary_function_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSECompareExtent
    //! @brief class to match two SSEs by extent
    //! @details The operator of this class returns whether two SSEs match with respect to their extent along a certain axis.
    //!
    //! @see @link example_assemble_sse_compare_extent.cpp @endlink
    //! @author linders
    //! @date Aug 17, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSECompareExtent :
      public util::BinaryFunctionInterface< SSE, SSE, bool>
    {

    private:

    //////////
    // data //
    //////////

      //! extent tolerance
      double m_ExtentTolerance;

      //! axis along which the extents are to be compared
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
      SSECompareExtent();

      //! @brief constructor from tolerance and axis
      //! @param TOLERANCE extent tolerance to be used
      //! @param AXIS axis to be used (e.g. z axis)
      SSECompareExtent( const double TOLERANCE, const coord::Axis &AXIS);

      //! @brief Clone function
      //! @return pointer to new SSEsMatchExtent
      virtual SSECompareExtent *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that tests whether the two sses have the same extent within some tolerance
      //! @param SSE_A first sse to be tested
      //! @param SSE_B second sse to be tested
      //! @return bool whether the two sses agree in their extent within some tolerance (m_ExtentTolerance)
      bool operator()( const SSE &SSE_A, const SSE &SSE_B) const;

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

    }; // class SSECompareExtent

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_SSE_COMPARE_EXTENT_H_ 
