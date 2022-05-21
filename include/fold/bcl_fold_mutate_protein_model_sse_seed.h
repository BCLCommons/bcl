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

#ifndef BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_SEED_H_
#define BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_SEED_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence_flexibility.h"
#include "find/bcl_find_locator_interface.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_range.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateProteinModelSSESeed
    //! @brief Add an adjacent sse next to a located sse of a given type, e.g. as a seed for growing loops
    //! @details
    //!
    //! @see @link example_fold_mutate_protein_model_sse_seed.cpp @endlink
    //! @author woetzen, alexanns
    //! @date Nov 5, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateProteinModelSSESeed :
      public math::MutateInterface< assemble::ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      //! locator that decides which sse in the proteinmodel to mutate
      util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> > m_SSELocator;

      //! range of number of residues that are to be added in one mutate
      math::Range< size_t> m_SeedLengthRange;

      //! side of located sse to add seed to
      biol::AASequenceFlexibility::DirectionEnum m_Direction;

      //! mutate to be applied to the added seed coil
      util::ShPtr< math::MutateInterface< assemble::SSE> > m_MutateSeed;

      //! range of amino acids to cut into the located sse
      math::Range< size_t> m_CutInRange;

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
      MutateProteinModelSSESeed();

      //! @brief constructor from a locator and other things
      //! @param SSE_LOCATOR locator that decides to which sse the seed sse is added
      //! @param SEED_LENGTH_RANGE number of residues to construct in seed sse, in addition to the cut
      //! @param DIRECTION the side of the located sse to which the seed sse is attached to
      //! @param SP_MUTATE_SEED optional mutate that is applied to the generated seed before it is added to the model
      //! @param CUT_IN_RANGE number of amino acids to cut into the located sse
      //! @param SCHEME Scheme to be used
      MutateProteinModelSSESeed
      (
        const util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> > &SSE_LOCATOR,
        const math::Range< size_t> &SEED_LENGTH_RANGE,
        const biol::AASequenceFlexibility::SequenceDirection &DIRECTION,
        const util::ShPtr< math::MutateInterface< assemble::SSE> > &SP_MUTATE_SEED,
        const math::Range< size_t> &CUT_IN_RANGE,
        const std::string &SCHEME = GetStaticClassName< MutateProteinModelSSESeed>()
      );

      //! @brief Clone function
      //! @return pointer to new MutateProteinModelSSESeed
      MutateProteinModelSSESeed *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this mutate
      //! @return the scheme for this mutate
      const std::string &GetScheme() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator taking an PROTEIN_MODEL and returning a protein model with an additional sse
      //! @param PROTEIN_MODEL protein model interest
      //! @return MutateResult with ProteinModel after the mutate
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

    }; // class MutateProteinModelSSESeed

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_SEED_H_ 
