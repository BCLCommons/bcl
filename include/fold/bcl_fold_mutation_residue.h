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

#ifndef BCL_FOLD_MUTATION_RESIDUE_H_
#define BCL_FOLD_MUTATION_RESIDUE_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutationResidue
    //! @brief is for holding information relevant for changing the dihedral angles of a residue of
    //! interest who might be mutated in some way, especially its phi and psi angles.
    //!
    //! @see @link example_fold_mutation_residue.cpp @endlink
    //! @author alexanns
    //! @date Sep 3, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutationResidue :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! the residue of interest
      util::ShPtr< biol::AABase> m_MutationResidue;

      //! the residue preceding the residue of interest in sequence
      util::ShPtr< biol::AABase> m_PreviousResidue;

      //! the residue proceeding the residue of interest in sequence
      util::ShPtr< biol::AABase> m_FollowingResidue;

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
      MutationResidue();

      //! @brief constructor taking member variables
      //! @param RESIDUE_TO_MUTATE the residue of interest
      //! @param RESIDUE_A residue to either side of the residue to mutate
      //! @param RESIDUE_B residue to either side of residue to mutate
      MutationResidue
      (
        const util::ShPtr< biol::AABase> &RESIDUE_TO_MUTATE,
        const util::ShPtr< biol::AABase> &RESIDUE_A,
        const util::ShPtr< biol::AABase> &RESIDUE_B
      );

      //! @brief constructor from iterator to residue to mutate and the map of all the residues
      //! @param RESIDUE_TO_MUTATE residue indicated to mutate
      //! @param ALL_RESIDUES all residues that could be mutated
      MutationResidue
      (
        const storage::List< storage::Pair< util::ShPtr< biol::AABase>, bool> >::const_iterator RESIDUE_TO_MUTATE,
        const storage::List< storage::Pair< util::ShPtr< biol::AABase>, bool> > &ALL_RESIDUES
      );

      //! @brief Clone function
      //! @return pointer to new MutationResidue
      MutationResidue *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief GetPreviousResidue gives the residue preceding the residue of interest in sequence
      //! @return the residue preceding the residue of interest in sequence
      const util::ShPtr< biol::AABase> &GetPreviousResidue() const;

      //! @brief GetFollowingResidue gives the residue proceeding the residue of interest in sequence
      //! @return the residue proceeding the residue of interest in sequence
      const util::ShPtr< biol::AABase> &GetFollowingResidue() const;

      //! @brief GetMutationResidue gives the residue of interest
      //! @return the residue of interest
      const util::ShPtr< biol::AABase> &GetMutationResidue() const;

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

    protected:

      //! @brief assigns the member variables to the proper residue
      //! @param RESIDUE_A first residue
      //! @param RESIDUE_B second residue
      void SortResidues
      (
        const util::ShPtr< biol::AABase> &RESIDUE_A, const util::ShPtr< biol::AABase> &RESIDUE_B
      );

    }; // class MutationResidue

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATION_RESIDUE_H_ 
