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

#ifndef BCL_FOLD_MUTATE_PROTEIN_MODEL_THREAD_SEQUENCE_H_
#define BCL_FOLD_MUTATE_PROTEIN_MODEL_THREAD_SEQUENCE_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "align/bcl_align.fwd.hh"
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_mutate_interface.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateProteinModelThreadSequence
    //! @brief uses an alignment to assign coordinates from a known structure onto a sequence of unknown structure
    //! @details A sequence alignment between a sequence with known coordinates and a sequence with unknown
    //!          coordinates is used to assign the coordinates from the given protein model to the second sequence in
    //!          the alignment. The protein model must have the same sequence as the first sequence in the alignment.
    //!
    //! @see @link example_fold_mutate_protein_model_thread_sequence.cpp @endlink
    //! @author alexanns
    //! @date Sep 14, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateProteinModelThreadSequence :
      public math::MutateInterface< assemble::ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      //! map to keep track of alignments for chains
      storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > m_Alignment;

      bool m_PrintStructuralProblemResidues;

      //! scheme
      std::string m_Scheme;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateProteinModelThreadSequence();

      //! @brief constructor taking member variable parameters
      //! @param CHAIN_ALIGNMENT for each chain, the alignment that should be used
      //! @param SCHEME the scheme of the mutate
      MutateProteinModelThreadSequence
      (
        const storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > &CHAIN_ALIGNMENT,
        const std::string &SCHEME = GetStaticClassName< MutateProteinModelThreadSequence>(),
        const bool PRINT_STRUCTURAL_PROBLEM_AAS = true
      );

      //! @brief Clone function
      //! @return pointer to new MutateProteinModelThreadSequence
      MutateProteinModelThreadSequence *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this mutate
      //! @return the scheme for this mutate
      const std::string &GetScheme() const;

    ////////////////
    // operations //
    ////////////////

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

    public:

      //! @brief uses a sequence alignment to assign coordinates from a chain onto a sequence with unknown structure
      //! @param CHAIN the coordinates that will be assigned to sequence of unknown structure
      //! @param ALIGNMENT the sequence alignment indicating which residues to assign which coordinates to
      //! @return shptr to chain which has the second sequence from the alignment but coordinates from the CHAIN param
      static util::ShPtr< assemble::Chain> ThreadSequenceOntoChain
      (
        const assemble::Chain &CHAIN, const align::AlignmentInterface< biol::AABase> &ALIGNMENT
      );

      //! @brief replaces a chain in a model
      //! @param MODEL the model in which the chain will be replaced
      //! @param CHAIN the chain which will replace existing chain in the provided model
      //! @return shptr to a new protein model with the desired chain replaced
      static util::ShPtr< assemble::ProteinModel>
      UpdateChainInModel
      (
        const assemble::ProteinModel &MODEL, const util::ShPtr< assemble::Chain> &CHAIN
      );

      //! @brief uses a sequence alignment to and chain to determine which stretches of unknown structure need rebuilt
      //!        This could be due e.g. to not having coordinates from the template structure, or if there are gaps in
      //!        the sequence of unknown structure (since this will likely lead to a chain break).
      //! @param CHAIN the coordinates that will be assigned to sequence of unknown structure
      //! @param ALIGNMENT the sequence alignment indicating which residues to assign which coordinates to
      //! @return string which has the regions of the sequence of unknown structure which most likely need rebuilt
      static std::string OutputStructurallyProbematicSequenceRegions
      (
        const assemble::Chain &CHAIN, const align::AlignmentInterface< biol::AABase> &ALIGNMENT
      );

      //! @brief determines if a residue has any defined coordinates
      //! @param RESI the residue that will be checked to see if it has any defined coordinates
      //! @return boolean true if there are any defined atom coordinates in RESI - false otherwise
      static bool HasAnyDefinedCoordinates( const biol::AABase &RESI);

    }; // class MutateProteinModelThreadSequence

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_PROTEIN_MODEL_THREAD_SEQUENCE_H_ 
