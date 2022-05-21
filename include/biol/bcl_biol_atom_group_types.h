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

#ifndef BCL_BIOL_ATOM_GROUP_TYPES_H_
#define BCL_BIOL_ATOM_GROUP_TYPES_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_biol_atom_group_type_data.h"
#include "util/bcl_util_enumerate.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomGroupTypes
    //! @brief Enumerator class to be used for accessing Atom Group type information in biology related classes.
    //! @details This enumerator class has an individual enum for each Atom Group type that are relevant to proteins
    //! thus it is more limited compared to chemistry::AtomTypes.
    //!
    //! @see @link example_biol_atom_group_types.cpp @endlink
    //! @author putnamdk
    //! @date 11/21/2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AtomGroupTypes :
      public util::Enumerate< AtomGroupTypeData, AtomGroupTypes>
    {
      friend class util::Enumerate< AtomGroupTypeData, AtomGroupTypes>;
    public:

    //////////
    // data //
    //////////

      // declare all atom types
      const AtomGroupType H;    //!< Hydrogen Atomic group
      const AtomGroupType C;    //!< Carbon Atomic group
      const AtomGroupType CH;   //!< Carbon with 1 Hydrogen
      const AtomGroupType CH2;  //!< Carbon with 2 Hydrogens
      const AtomGroupType CH3;  //!< Carbon with 3 Hydrogens
      const AtomGroupType N;    //!< Nitrogen Atomic group
      const AtomGroupType NH;   //!< Nitrogen with 1 Hydrogen
      const AtomGroupType NH2;  //!< Nitrogen with 2 Hydrogens
      const AtomGroupType NH3;  //!< Nitrogen with 3 Hydrogens
      const AtomGroupType O;    //!< Oxygen Atomic group
      const AtomGroupType OH;   //!< Oxygen with 1 Hydrogen
      const AtomGroupType S;    //!< Sulfur Atomic group
      const AtomGroupType SH;   //!< Sulfur with 1 Hydrogen
      const AtomGroupType SE;   //!< Selenium Atomic group

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct all AtomTypes
      AtomGroupTypes();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief map to hold Solvent Volume and Radius values for each heavy atom in each amino acid
      //! @param AA The amino acid types
      //! @param ATOMTYPE The atomtype
      //! @return Group Types
      //! The values are taken from Crysol - a program to evaluate x-ray solution scattering of biological
      //! macromolecules from Atomic Coordinates.  D. Svergun et all 1995
      static const AtomGroupType &GetType( const AAType &AA, const AtomType &ATOMTYPE);

      //! @brief map to hold Solvent Volume and Radius values for each heavy atom in each amino acid
      //! @param AA The amino acid types
      //! @param ATOMTYPE The atomtype
      //! @return string of group types
      //! The values are taken from Crysol - a program to evaluate x-ray solution scattering of biological
      //! macromolecules from Atomic Coordinates.  D. Svergun et all 1995
      static const std::string &GetTypeString( const AAType &AA, const AtomType &ATOMTYPE);

    }; // class AtomGroupTypes

    //! @brief access to the only instance of AtomGroupTypes
    //! @return reference to only instance of AtomGroupTypes
    BCL_API
    const AtomGroupTypes &GetAtomGroupTypes();

  } // namespace biol

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< biol::AtomGroupTypeData, biol::AtomGroupTypes>;

  } // namespace util
} // namespace bcl

#endif //BCL_BIOL_ATOM_GROUP_TYPES_H_
