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

#ifndef BCL_ASSEMBLE_SSE_TRANSFORMER_H_
#define BCL_ASSEMBLE_SSE_TRANSFORMER_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_function_interface_serializable.h"
#include "math/bcl_math_transformation_matrix_3d.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSETransformer
    //! @brief Assigns a new sequence and transformation to an SSE
    //! @details () operator clones an SSE, assigns a new sequence to it (same residues but potentially different chain
    //!          id and coordinates), and applies a transformation.  This is used by ChainMultiplier to generate new
    //!          SSEs for each new chain.
    //!
    //! @see @link example_assemble_sse_transformer.cpp @endlink
    //! @author weinerbe
    //! @date Nov 12, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSETransformer :
      public math::FunctionInterfaceSerializable< SSE, util::ShPtr< SSE> >
    {

    private:

    //////////
    // data //
    //////////

      //! associated sequence
      util::ShPtr< biol::AASequence> m_Sequence;

      //! transformation to apply to sses
      util::ShPtr< math::TransformationMatrix3D> m_Transformation;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SSETransformer();

      //! @brief construct from a transformation matrix and chain IDs
      //! @param SP_SEQUENCE sequence used to generate new SSE
      //! @param SP_TRANSFORMATION the transformation to apply to sse
      SSETransformer
      (
        const util::ShPtr< biol::AASequence> &SP_SEQUENCE,
        const util::ShPtr< math::TransformationMatrix3D> &SP_TRANSFORMATION
      );

      //! @brief Clone function
      //! @return pointer to new SSETransformer
      SSETransformer *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief returns SSE after transformation
      //! @param INITIAL_SSE SSE to be transformed
      //! @return a cloned SSE with associated AAData to the member sequence and transformed
      util::ShPtr< SSE> operator ()( const SSE &INITIAL_SSE) const;

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

    }; // class SSETransformer

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_SSE_TRANSFORMER_H_ 
