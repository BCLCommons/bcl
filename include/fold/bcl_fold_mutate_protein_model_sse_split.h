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

#ifndef BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_SPLIT_H_
#define BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_SPLIT_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_ss_types.h"
#include "find/bcl_find_locator_interface.h"
#include "math/bcl_math_mutate_interface.h"
#include "sspred/bcl_sspred_methods.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateProteinModelSSESplit
    //! @brief mutates a pool by picking an sse and splitting it into two SSEs at the largest drop in the ss prediction
    //!
    //! @see @link example_fold_mutate_protein_model_sse_split.cpp @endlink
    //! @author karakam, woetzen
    //! @date Aug 3, 2011
    //!
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////d///

    class BCL_API MutateProteinModelSSESplit :
      public math::MutateInterface< assemble::ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      //! range for how many residues from the middle should be split
      sspred::Method m_Method;

      //! locator of an sse in a domain
      util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> > m_SSELocator;

      //! min sse sizes for the resultant SSEs
      storage::Map< biol::SSType, size_t> m_MinSSESizes;

      //! length of the coil to be inserted after split
      math::Range< size_t> m_SplitCoilLengthRange;

      //! scheme for this mutate
      std::string m_Scheme;

    public:

      //! @brief the default scheme for this class
      static const std::string &GetDefaultScheme();

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from locator and mutate
      //! @param SS_METHOD the method to use to locate the largest drop in the ss prediction for the located sse
      //! @param SP_LOCATE_SSE locator of sse from domain
      //! @param MIN_SSE_SIZES min sse sizes for the resultant SSEs
      //! @param SPLIT_COIL_LENGTH_RANGE range of length of the coil to be inserted after the split
      //! @param SCHEME the scheme
      MutateProteinModelSSESplit
      (
        const sspred::Method &SS_METHOD,
        const util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> > &SP_LOCATE_SSE,
        const storage::Map< biol::SSType, size_t> &MIN_SSE_SIZES,
        const math::Range< size_t> &SPLIT_COIL_LENGTH_RANGE,
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief Clone function
      //! @return pointer to new MutateProteinModelSSESplit
      MutateProteinModelSSESplit *Clone() const;

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

      //! @brief operator that mutates a single SSE in a given ProteinModel by splitting a single sse
      //! @param PROTEIN_MODEL ProteinModel from which a single SSE will be splitted
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

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief find the amino acid right to the highest difference for the given sspred::Method and SSType and report the difference
      //! @param ELEMENT the secondary structure element
      //! @param SSMETHOD the ssprediction method
      //! @param NR_RESIDUES_TO_EXCLUDE number of residues to exclude from both sides
      //! @return Pair of difference and the amino acid to the left
      static storage::Pair< util::SiPtr< const biol::AABase>, double> FindLargestSSPredDifferencePosition
      (
        const assemble::SSE &ELEMENT,
        const sspred::Method &SSMETHOD,
        const size_t NR_RESIDUES_TO_EXCLUDE
      );

      //! @brief find the amino acid right to the highest difference for the given sspred::Method and SSType and report the difference
      //! @param ELEMENT the secondary structure element
      //! @param SSMETHOD the ssprediction method
      //! @param NR_RESIDUES_TO_EXCLUDE_LEFT number of residues to exclude from the left side
      //! @param NR_RESIDUES_TO_EXCLUDE_RIGHT number of residues to exclude from the right side
      //! @return Pair of difference and the amino acid to the left
      static util::SiPtr< const biol::AABase> FindSplitPosition
      (
        const assemble::SSE &ELEMENT,
        const sspred::Method &SSMETHOD,
        const size_t NR_RESIDUES_TO_EXCLUDE_LEFT,
        const size_t NR_RESIDUES_TO_EXCLUDE_RIGHT
      );

    }; // class MutateProteinModelSSESplit

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_SPLIT_H_ 
