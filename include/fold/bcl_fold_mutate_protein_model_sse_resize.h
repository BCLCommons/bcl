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

#ifndef BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_RESIZE_H_
#define BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_RESIZE_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence_flexibility.h"
#include "find/bcl_find_locator_interface.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_range.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateProteinModelSSEResize
    //! @brief resizes SSEs in the given Protein Model
    //! @details This class allows extending and shrinking the lengths of SSEs in the model, thus allowing to cover a larger
    //! space with a smaller pool
    //!
    //! @see @link example_fold_mutate_protein_model_sse_resize.cpp @endlink
    //! @author karakam
    //! @date Oct 16, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API MutateProteinModelSSEResize :
      public math::MutateInterface< assemble::ProteinModel>
    {

    protected:

    //////////
    // data //
    //////////

      //! locator that decides which sse in the proteinmodel to mutate
      util::Implementation< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> > m_SSELocator;

      //! probability for extending (1.0-m_ExtendProbability is the shrink probability
      double m_ExtendProbability;

      //! range of number of residues that are to be added or removed in one mutate to one end [min,max]
      math::Range< size_t> m_LengthChangeRange;

      //! side of sequence to be modified
      biol::AASequenceFlexibility::DirectionEnum m_Side;

      //! boolean to whether recenter the sse after resize to its original center
      bool m_RecenterAfterResize;

      //! minimum SSE sizes to be allowed when shrinking
      storage::Map< biol::SSType, size_t> m_MinSSESizes;

      //! scheme
      std::string m_Scheme;

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
      MutateProteinModelSSEResize();

      //! @brief constructor from a locator, an extend/shrink probability, length increment and a boolean flag
      //! @param SSE_LOCATOR locator that decides which SSE in the protein model to mutate
      //! @param EXTEND_PROBABILITY probability for extending (1.0 - EXTEND_PROBABILITY for shrinking)
      //! @param LENGTH_CHANGE_RANGE range of number of residues that are to be added or removed in one mutate to one end
      //! @param SEQUENCE_SIDE side of sequence to modify
      //! @param RECENTER_AFTER_RESIZE boolean to whether recenter the SSE after resize to its original center
      //! @param MIN_SSE_SIZES map of minimum SSE sizes to be allowed when shrinking
      //! @param SCHEME Scheme to be used
      MutateProteinModelSSEResize
      (
        const util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> > &SSE_LOCATOR,
        const double &EXTEND_PROBABILITY,
        const math::Range< size_t> &LENGTH_CHANGE_RANGE,
        const biol::AASequenceFlexibility::SequenceDirection &SEQUENCE_SIDE,
        const bool RECENTER_AFTER_RESIZE,
        const storage::Map< biol::SSType, size_t> &MIN_SSE_SIZES,
        const std::string &SCHEME = GetStaticClassName< MutateProteinModelSSEResize>()
      );

      //! @brief Clone function
      //! @return pointer to new MutateProteinModelSSEResize
      MutateProteinModelSSEResize *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return extend probability
      //! @return extend probability
      const double &GetExtendProbability() const;

      //! @brief return length change range
      //! @return length change range
      const math::Range< size_t> &GetLengthChangeRange() const;

      //! @brief return which ends are changed
      //! @return sequence direction
      biol::AASequenceFlexibility::SequenceDirection GetSide() const;

      //! @brief returns min sse sizes
      //! @return min sse sizes
      const storage::Map< biol::SSType, size_t> &GetMinSSESizes() const;

      //! @brief gets the scheme for this mutate
      //! @return the scheme for this mutate
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
      //! @param PROTEIN_MODEL protein model interest
      //! @return MutateResult with ProteinModel after the mutate
      math::MutateResult< assemble::ProteinModel> operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // class MutateProteinModelSSEResize

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_RESIZE_H_
