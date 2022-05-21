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

#ifndef BCL_FOLD_MUTATE_PROTEIN_MODEL_MOVE_AA_H_
#define BCL_FOLD_MUTATE_PROTEIN_MODEL_MOVE_AA_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence_flexibility.h"
#include "find/bcl_find_locator_interface.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_range.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateProteinModelMoveAA
    //! @brief move an amino acid between two sequence-adjacent sses
    //! @details Initially, an SSE is located. Afterwards, adjacent SSEs are found in the protein model. A number of
    //!          amino acids is removed from one SSE and added to the second one.
    //!
    //! @see @link example_fold_mutate_protein_model_move_aa.cpp @endlink
    //! @author woetzen
    //! @date Nov 3, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateProteinModelMoveAA :
      public math::MutateInterface< assemble::ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      //! locator that decides which SSE in the protein model to move aas from or to
      util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> > m_SSELocator;

      //! max number of residues to move
      math::Range< size_t> m_ResdiuesToMoveRange;

      //! side on located sse where adjacent sse is located and aas are moved
      biol::AASequenceFlexibility::DirectionEnum m_SSESide;

      //! min sse size - after moving aas, that sses should not be below the given length
      storage::Map< biol::SSType, size_t> m_MinSSESize;

      //! mutate to be applied to the previously located sse
      util::ShPtr< math::MutateInterface< assemble::SSE> > m_MutateLocatedSSE;

      //! scheme
      std::string m_Scheme;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from max number of residues to move
      //! @param SSE_LOCATOR locator that decides which SSE in the protein model to find an adjacent one to
      //! @param RESIDUES_TO_MOVE_RANGE range of number of residues to move
      //! @param SEQUENCE_SIDE side on sequence to move aas
      //! @param MIN_SSE_SIZE minimum SSE sizes to be allowed as result after moving aas
      //! @param SP_MUTATE_LOCATED_SSE optional mutate that is applied to the located and modified sse before it is added to the model
      //! @param SCHEME the scheme
      MutateProteinModelMoveAA
      (
        const util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> > &SSE_LOCATOR,
        const math::Range< size_t> &RESIDUES_TO_MOVE_RANGE,
        const biol::AASequenceFlexibility::SequenceDirection &SEQUENCE_SIDE,
        const storage::Map< biol::SSType, size_t> &MIN_SSE_SIZE,
        const util::ShPtr< math::MutateInterface< assemble::SSE> > &SP_MUTATE_LOCATED_SSE,
        const std::string &SCHEME
      );

      //! @brief Clone function
      //! @return pointer to new MutateProteinModelMoveAA
      MutateProteinModelMoveAA *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this mutate
      //! @return the scheme for this mutate
      const std::string &GetScheme() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that mutates a protein model by shifting AAs between SSEs
      //! @param PROTEIN_MODEL
      //! @return MutateResult that results from mutating to the PROTEIN_MODEL
      math::MutateResult< assemble::ProteinModel> operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

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

    public:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief move a given number of amino acids in a given direction between two adjacent sses
      //! @param SSE_LEFT n-terminal SSE
      //! @param SSE_RIGHT c-terminal SSE
      //! @param NUMBER_AAS the number of amino acids to move from one sse to another
      //! @param DIRECTION the sequence direstion - c-terminal move from left to right, n-terminal right to left
      //! @return 2 new SSE with the same type as the given ones - if the number of AAs was larger than the size, the
      //!         ShPtr will be undefined
      static storage::VectorND< 2, util::ShPtr< assemble::SSE> > MoveAAs
      (
        const assemble::SSE &SSE_LEFT,
        const assemble::SSE &SSE_RIGHT,
        const size_t NUMBER_AAS,
        const biol::AASequenceFlexibility::SequenceDirection &DIRECTION
      );

    }; // class MutateProteinModelMoveAA

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_PROTEIN_MODEL_MOVE_AA_H_
