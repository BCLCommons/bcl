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

#ifndef BCL_ASSEMBLE_CHAIN_MULTIPLIER_H_
#define BCL_ASSEMBLE_CHAIN_MULTIPLIER_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"
#include "math/bcl_math_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ChainMultiplier
    //! @brief Creates a new chain from the given sequence and passed existing chain
    //! @details () operator creates a new chain from the sequence and copies the SSEs from the passed chain using an
    //!          SSE transformer.  This class is used by ProteinMultiplier to generate protein multimers.
    //!
    //! @see @link example_assemble_chain_multiplier.cpp @endlink
    //! @author weinerbe
    //! @date Nov 12, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ChainMultiplier :
      public math::FunctionInterfaceSerializable< Chain, util::ShPtr< Chain> >
    {

    private:

    //////////
    // data //
    //////////

      //! SSE transformer to use
      util::ShPtr< util::FunctionInterface< SSE, util::ShPtr< SSE> > > m_Transformer;

      //! initial chain id
      char m_InitialChainID;

      //! transformation that is applied to SSEs by the transformer
      util::ShPtr< math::TransformationMatrix3D> m_Transformation;

      //! new sequence for new chain
      util::ShPtr< biol::AASequence> m_Sequence;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ChainMultiplier();

      //! @brief construct from an SSE transformer, chain ID, and sequence
      //! @param SSE_TRANSFORMER transformer to apply to each SSE
      //! @param INITIAL_CHAIN_ID initial chain id
      //! @param TRANSFORMATION transformation that is applied to each SSE
      //! @param SEQUENCE sequence used to generate new chain
      ChainMultiplier
      (
        const util::ShPtr< util::FunctionInterface< SSE, util::ShPtr< SSE> > > &SSE_TRANSFORMER,
        const char INITIAL_CHAIN_ID,
        const util::ShPtr< math::TransformationMatrix3D> &TRANSFORMATION,
        const util::ShPtr< biol::AASequence> &SEQUENCE
      );

      //! @brief Clone function
      //! @return pointer to new ChainMultiplier
      ChainMultiplier *Clone() const;

      //! @brief hardcopy this multiplier
      //! @return harcopied multiplier
      util::ShPtr< ChainMultiplier> HardCopy() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the initial chain ID
      //! @return the initial chain ID
      char GetInitialChainID() const
      {
        return m_InitialChainID;
      }

      //! @brief get the new chain ID
      //! @return the new chain ID
      char GetNewChainID() const
      {
        return m_Sequence->GetChainID();
      }

      //! @brief gets the transformation matrix
      //! @return a ShPtr to the transformation matrix
      const util::ShPtr< math::TransformationMatrix3D> &GetTransformationMatrix() const
      {
        return m_Transformation;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief returns chain after transforming all SSEs
      //! @param CHAIN chain to be replicated
      //! @return chain after transforming all SSEs
      util::ShPtr< Chain> operator ()( const Chain &CHAIN) const;

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

    }; // class ChainMultiplier

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ChainMultiplierLessThan
    //! @brief has operator for returning whether one chain multiplier is less than another
    //!
    //! @remarks example unnecessary
    //! @author weinerbe
    //! @date Nov 12, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct BCL_API ChainMultiplierLessThan
    {
      //! @brief return whether one chain multiplier is less than another
      //! @param CHAIN_MULTIPLIER_LHS first chain multiplier
      //! @param CHAIN_MULTIPLIER_RHS second chain multiplier
      //! @return whether one chain multiplier is less than another
      bool operator()
      (
        const ChainMultiplier &CHAIN_MULTIPLIER_LHS,
        const ChainMultiplier &CHAIN_MULTIPLIER_RHS
      ) const;

      //! @brief return whether one chain multiplier is less than another
      //! @param PTR_CHAIN_MULTIPLIER_LHS first chain multiplier
      //! @param PTR_CHAIN_MULTIPLIER_RHS second chain multiplier
      //! @return whether one chain multiplier is less than another
      bool operator()
      (
        const util::PtrInterface< ChainMultiplier> &PTR_CHAIN_MULTIPLIER_LHS,
        const util::PtrInterface< ChainMultiplier> &PTR_CHAIN_MULTIPLIER_RHS
      ) const;

    }; // struct ChainMultiplierLessThan

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_CHAIN_MULTIPLIER_H_ 
