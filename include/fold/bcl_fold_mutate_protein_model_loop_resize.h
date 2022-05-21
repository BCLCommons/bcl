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

#ifndef BCL_FOLD_MUTATE_PROTEIN_MODEL_LOOP_RESIZE_H_
#define BCL_FOLD_MUTATE_PROTEIN_MODEL_LOOP_RESIZE_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"
#include "find/bcl_find.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_fold_mutate_protein_model_sse_resize.h"
#include "biol/bcl_biol_aa_sequence_flexibility.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateProteinModelLoopResize
    //! @brief resizes an sse - if residues are added the phi psi angles are given by a phi psi generator object
    //! @details Extends or shrinks an sse with a given probability a desired number of times randomly chosen from a
    //!          range. A maximum and minimum size possible for the sse is enforced. This class differs from
    //!          MutateProteinModelSSEResize because that class is hard coded to use ideal geometry to extend an sse.
    //!
    //! @see @link example_fold_mutate_protein_model_loop_resize.cpp @endlink
    //! @author alexanns
    //! @date Sep 6, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateProteinModelLoopResize :
      public MutateProteinModelSSEResize
    {

    private:

    //////////
    // data //
    //////////

      //! maximum SSE sizes to be allowed when growing
      storage::Map< biol::SSType, size_t> m_MaxSSESizes;

      //! the method that will be used in order to generate phi and psi angles as the loop is grown
      util::ShPtr< math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> > > m_PhiPsiGenerator;

      //! boolean true if extension of sse into neighboring sse is not allowed - false otherwise (resplaces both sses)
      bool m_DisallowOverlap;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateProteinModelLoopResize();

      //! @brief constructor taking member variable parameters
      //! @param LOCATOR locator that decides which sse in the proteinmodel to mutate
      //! @param EXTEND_PROBABILITY probability for extending (1.0 - EXTEND_PROBABILITY for shrinking)
      //! @param LENGTH_CHANGE_RANGE range of number of residues to be added or removed in one mutate to one end [min,max]
      //! @param MIN_SSE_SIZES minimum SSE sizes to be allowed when shrinking
      //! @param MAX_SSE_SIZES minimum SSE sizes to be allowed when extending
      //! @param GROWING_DIRECTION the side of the sse that changes will occur on
      //! @param PHI_PSI_GENERATOR the method that will be used in order to generate phi and psi angles if extending
      //! @param GROWING_DIRECTION the side of the sse that changes will occur on
      //! @param DISALLOW_OVERLAP true if extension of sse into neighboring sse is not allowed-else (replaces both sses)
      //! @param SCHEME the scheme of this mutate
      MutateProteinModelLoopResize
      (
        const util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> > &LOCATOR,
        const double &EXTEND_PROBABILITY,
        const math::Range< size_t> &LENGTH_CHANGE_RANGE,
        const storage::Map< biol::SSType, size_t> &MIN_SSE_SIZES,
        const storage::Map< biol::SSType, size_t> &MAX_SSE_SIZES,
        const biol::AASequenceFlexibility::SequenceDirection &GROWING_DIRECTION,
        const util::ShPtr< math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> > > &PHI_PSI_GENERATOR,
        const bool DISALLOW_OVERLAP = true,
        const std::string &SCHEME = GetStaticClassName< MutateProteinModelLoopResize>()
      );

      //! @brief Clone function
      //! @return pointer to new MutateProteinModelLoopResize
      MutateProteinModelLoopResize *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator taking an ARGUMENT and returning a mutated object of t_ArgumentType
      //! @param PROTEIN_MODEL Argument of interest
      //! @return MutateResult that results from mutating to the argument
      math::MutateResult< assemble::ProteinModel> operator()
      (
        const assemble::ProteinModel &PROTEIN_MODEL
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

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief prepends an aa sequence to the nterminal side of an sse
      //! @param SSE the sse to which a sequence will be prepended
      //! @param SEQUENCE the sequence that will be prepended to the sse
      //! @return sse which has been created by prepending the provided sequence onto n-terminus of the sse
      util::ShPtr< assemble::SSE> PrependSequence( const assemble::SSE &SSE, const biol::AASequence &SEQUENCE) const;

      //! @brief appends an aa sequence to the cterminal side of an sse
      //! @param SSE the sse to which a sequence will be appended
      //! @param SEQUENCE the sequence that will be appended to the sse
      //! @return sse which has been created by appending the provided sequence onto c-terminus of the sse
      util::ShPtr< assemble::SSE> AppendSequence( const assemble::SSE &SSE, const biol::AASequence &SEQUENCE) const;

      static bool WillOverlap
      (
        const assemble::ProteinModel &MODEL, const size_t EXTENSION_AMOUNT, const assemble::SSE &CURRENT_SSE,
        const biol::AASequenceFlexibility::SequenceDirection &EXTENSION_DIRECTION
      );

    }; // class MutateProteinModelLoopResize

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_PROTEIN_MODEL_LOOP_RESIZE_H_ 
