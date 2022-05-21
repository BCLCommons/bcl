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

#ifndef BCL_CHEMISTRY_SMALL_MOLECULE_FRAGMENT_MAPPING_H_
#define BCL_CHEMISTRY_SMALL_MOLECULE_FRAGMENT_MAPPING_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_rotamer_dihedral_bond_data.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SmallMoleculeFragmentMapping
    //! @brief helper class that creates RotamerDihedralBondData object for fragments of interest.
    //! @details creates Rotamer information about a fragment from the properties stored on the fragment
    //!
    //! @see @link example_chemistry_small_molecule_fragment_mapping.cpp @endlink
    //! @author kothiwsk
    //! @date Jul 17, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SmallMoleculeFragmentMapping :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new SmallMoleculeFragmentMapping
      SmallMoleculeFragmentMapping *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operators //
    ////////////////

      //! @brief get all possible mappings between fragments identified and the molecule of interest
      //! @param FRAGMENTS fragment isomorphisms of fragments contained in a molecule of interest
      //! @return map of number of rotatable bond of fragment to its RotateDihedralBondData
      util::ShPtrVector< RotamerDihedralBondData> MapFragmentIsomorphisms
      (
        const FragmentComplete &MOLECULE,
        const util::ShPtrVector< SmallMoleculeFragmentIsomorphism> &FRAGMENTS,
        const storage::Set< size_t> &SAMPLE_PARTS = storage::Set< size_t>()
      );

      //! @brief Filter out RDBD, choosing the best one to represent each case
      //! @param FRAGMENTS fragment isomorphisms of fragments contained in a molecule of interest
      //! @return map of number of rotatable bond of fragment to its RotateDihedralBondData
      util::ShPtrVector< RotamerDihedralBondData> Filter
      (
        const util::ShPtrVector< RotamerDihedralBondData> &ROTAMERS
      );

    ////////////////
    // operations //
    ////////////////

      //! @brief get all possible rings contained in the molecule of  interest
      //! @param MOLECULE molecule whose rings have to be found
      //! @return list of aotm vertices that are contained in rings
      static const storage::List< storage::Vector< size_t> > GetAllRings( const ConformationInterface &MOLECULE);

      //! @brief get all possible rings contained in the molecule of  interest
      //! @param MOLECULE molecule whose rings have to be found
      //! @return list of aotm vertices that are contained in rings
      static const storage::List< storage::Vector< size_t> > GetRingVertices( const ConformationInterface &MOLECULE);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class SmallMoleculeFragmentMapping

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_SMALL_MOLECULE_FRAGMENT_MAPPING_H_
