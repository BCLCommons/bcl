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

#ifndef BCL_ASSEMBLE_SSE_PAIR_TEMPLATE_H_
#define BCL_ASSEMBLE_SSE_PAIR_TEMPLATE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_assemble_sse_geometry_interface.h"
#include "bcl_assemble_sse_geometry_packing.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEPairTemplate
    //! @brief This class stores bodies of SSE pairs to be used for placements in protein folding
    //! @details This class describes a template for an SSE pair.  A sequential SSE can be placed based on statistics
    //!          generated for similar instances of SSE types and loop lengths found in the PDB.
    //!
    //! @see @link example_assemble_sse_pair_template.cpp @endlink
    //! @author weinerbe
    //! @date Oct 1, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEPairTemplate :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! First SSE geometry
      util::ShPtr< SSEGeometryInterface> m_FirstSSEGeometry;

      //! Second SSE geometry
      util::ShPtr< SSEGeometryInterface> m_SecondSSEGeometry;

      //! loop length
      size_t m_LoopLength;

      //! packing object
      SSEGeometryPacking m_Packing;

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
      SSEPairTemplate();

      //! @brief construct from two SSE geometries
      //! @param SP_FIRST_SSE_GEOMETRY first SSE geometry
      //! @param SP_SECOND_SSE_GEOMETRY second SSE geometry
      //! @param LOOP_LENGTH length of the loop connecting the two geometries
      SSEPairTemplate
      (
        const util::ShPtr< SSEGeometryInterface> &SP_FIRST_SSE_GEOMETRY,
        const util::ShPtr< SSEGeometryInterface> &SP_SECOND_SSE_GEOMETRY,
        const size_t LOOP_LENGTH
      );

      //! @brief Clone function
      //! @return pointer to new SSEPairTemplate
      SSEPairTemplate *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief gets the first geometry
      //! @return the first geometry
      const util::ShPtr< SSEGeometryInterface> &GetFirstGeometry() const
      {
        return m_FirstSSEGeometry;
      }

      //! @brief gets the second geometry
      //! @return the second geometry
      const util::ShPtr< SSEGeometryInterface> &GetSecondGeometry() const
      {
        return m_SecondSSEGeometry;
      }

      //! @brief get the number of residues between SSE_A and SSE_B
      //! @return the number of residues between SSE_A and SSE_B
      size_t GetLoopLength() const
      {
        return m_LoopLength;
      }

      //! @brief get the packing defined by the two geometries
      //! @return the the packing defined by the two geometries
      const SSEGeometryPacking &GetPacking() const
      {
        return m_Packing;
      }

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

    }; // class SSEPairTemplate

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_SSE_PAIR_TEMPLATE_H_ 
