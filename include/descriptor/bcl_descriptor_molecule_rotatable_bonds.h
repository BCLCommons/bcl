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

#ifndef BCL_DESCRIPTOR_MOLECULE_ROTATABLE_BONDS_H_
#define BCL_DESCRIPTOR_MOLECULE_ROTATABLE_BONDS_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "chemistry/bcl_chemistry.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_sequence.h"
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeRotatableBonds
    //! @brief calculates the # of rotatable bonds in the molecule
    //! @details Rotatable bond definition is from:
    //!   Veber, et. al. 2002, Molecular Properties That Influence the Oral Bioavailability of Drug Candidates
    //!
    //!   Single
    //!   Not part of a ring
    //!   Not terminal (e.g. C-H)
    //!   Not connecting to an atom that has no other non-hydrogen bonds (e.g. C-CH3)
    //!   Not a C-N bond in an amide-like group (C float  bonded to a terminal group, single bonded to a nitrogen)
    //!
    //! @see @link example_descriptor_molecule_rotatable_bonds.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Feb 13, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeRotatableBonds :
      public BaseSequence< chemistry::AtomConformationalInterface, float>
    {

      bool m_SymmetryAware; //!< Whether to consider symmetry

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MoleculeRotatableBonds( const bool &SYMMETRY_AWARE = false) :
        m_SymmetryAware( SYMMETRY_AWARE)
      {
      }

      //! @brief Clone function
      //! @return pointer to new MoleculeRotatableBonds
      MoleculeRotatableBonds *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get name of the current class
      //! @return name of the class
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature siz  e under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const
      {
        return 1;
      }

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief calculate the descriptors
      //! @param STORAGE storage for the descriptor
      void Calculate( linal::VectorReference< float> &STORAGE);

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const
      {
        return e_PreferCache;
      }

    }; // class MoleculeRotatableBonds

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MOLECULE_ROTATABLE_BONDS_H_
