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

#ifndef BCL_CHEMISTRY_AA_FRAGMENT_COMPLETE_H_
#define BCL_CHEMISTRY_AA_FRAGMENT_COMPLETE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_fragment_complete.h"
#include "biol/bcl_biol_aa_sequence.h"
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AAFragmentComplete
    //! @brief Class that represents a fragment that is an amino acid or amino acid sequence
    //!
    //! @see @link example_chemistry_aa_fragment_complete.cpp @endlink
    //! @author mendenjl
    //! @date Sep 12, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AAFragmentComplete :
      public FragmentComplete
    {

    private:
    //////////
    // data //
    //////////

      storage::Vector< biol::AtomType>       m_BiolAtomType;    //!< Biol atom type, one per Atom
      util::SiPtrVector< const biol::AABase> m_Sequence;        //!< Sequence that is represented; one pointer per residue
      util::SiPtrVector< const biol::AABase> m_AtomAASequence;  //!< Sequence that is represented; one pointer per atom
      storage::Vector< size_t>               m_AAIndices;       //!< Indices of the AA for each atom

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      AAFragmentComplete();

      //! @brief construct from an AA sequence.  Normally the contained AA class should be AA complete
      //! @param ALLOW_UNDEFINED_POSITIONS true to allow undefined atom positions
      //! Default is to assert if heavy atom position is absent
      AAFragmentComplete
      (
        const biol::AASequence &SEQUENCE,
        const bool &ALLOW_UNDEFINED_POSITIONS = false
      );

      //! @brief construct from an AA sequence.  Normally the contained AA type should be AA complete
      //! @param ALLOW_UNDEFINED_POSITIONS true to allow undefined atom positions
      //! Default is to assert if heavy atom position is absent
      AAFragmentComplete
      (
        const util::SiPtrVector< const biol::AABase> &AA_SEQUENCE,
        const bool &ALLOW_UNDEFINED_POSITIONS = false
      );

      //! @brief Clone function
      //! @return pointer to new MolecularConformationShared
      AAFragmentComplete *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get the number of valences
      //! @return the number of valences (=0)
      size_t GetNumberValences() const
      {
        return 0;
      }

      //! @brief get the number of valences
      //! @return the number of valences (=0)
      void RemoveAtomsUndefinedPos();

      //! @brief get the underlying sequence
      //! @return the underlying sequence of AAs
      const util::SiPtrVector< const biol::AABase> &GetResidueSequence() const
      {
        return m_Sequence;
      }

      //! @brief get the corresponding residue given an atom
      //! @param ATOM_ID id of the desired atom
      const biol::AABase &GetAtomsResidue( const size_t &ATOM_ID) const
      {
        return *m_AtomAASequence( ATOM_ID);
      }

      //! @brief get the corresponding residue given an atom
      //! @param ATOM_ID id of the desired atom
      const biol::AtomType &GetAtomsType( const size_t &ATOM_ID) const
      {
        return m_BiolAtomType( ATOM_ID);
      }

      //! @brief get the corresponding residue indices (relative to the original sequence)
      const storage::Vector< size_t> &GetAtomsResidueIndices() const
      {
        return m_AAIndices;
      }

      //! @brief return each residue as a separate fragment
      //! @param CONSIDER_BACK_BONE whether to include the back bone atoms
      //! @param CONSIDER_SIDE_CHAIN whether to include the side chain atoms
      FragmentEnsemble GetResiduesAsFragments( const bool &CONSIDER_BACK_BONE, const bool &CONSIDER_SIDE_CHAIN) const;

      //! @brief reconstruct protein model
      //! @param MODEL the original protein model; needed to get the seqres
      //! @note the SSEs present in the model will be ignored entirely
      assemble::ProteinModel ReconstructProteinModel( const assemble::ProteinModel &MODEL);

    private:

      //! @brief construct this object from a siptr vector of AAs
      //! @param SEQUENCE the sequence of AAs of interest
      //! @param ALLOW_UNDEFINED_POSITIONS true to allow undefined atom positions
      //! Default is to assert if heavy atom position is absent
      void Construct
      (
        const util::SiPtrVector< const biol::AABase> &AA_SEQUENCE,
        const bool &ALLOW_UNDEFINED_POSITIONS
      );

    }; // class AAFragmentComplete

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_AA_FRAGMENT_COMPLETE_H_
