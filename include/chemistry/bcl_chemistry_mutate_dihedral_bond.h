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

#ifndef BCL_CHEMISTRY_MUTATE_DIHEDRAL_BOND_H_
#define BCL_CHEMISTRY_MUTATE_DIHEDRAL_BOND_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_interface.h"
#include "bcl_chemistry_mutate_dihedrals_interface.h"
#include "bcl_chemistry_ring_fragment_map.h"
#include "bcl_chemistry_rotamer_dihedral_bond_data.h"
#include "coord/bcl_coord_line_segment_3d.h"
#include "graph/bcl_graph_const_graph.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"
#include "random/bcl_random_histogram_1d_distribution.h"
#include "storage/bcl_storage_vector.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateDihedralBond
    //! @brief This class is for randomly using different rotamers of a fragment to change conformation of a molecule of interest.
    //! @details
    //!
    //! @see @link example_chemistry_mutate_dihedral_bond.cpp @endlink
    //! @author mendenjl
    //! @date Oct 05, 2018
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateDihedralBond :
      public MutateDihedralsInterface
    {

    private:

    //////////
    // data //
    //////////

        //! Atom indices involved in the dihedral. Atoms A and D are chosen arbitrarily
        size_t                                                                               m_AtomIndexA;
        size_t                                                                               m_AtomIndexB;
        size_t                                                                               m_AtomIndexC;
        size_t                                                                               m_AtomIndexD;

        //! if true, only flip the bond (true for amides or double bonds (only if change chirality is given))
        bool                                                                                 m_FlipOnly;

        //! if true, only wiggle the dihedral bond (stay approximately within the same dihedral bin but just wiggle according
        //! to the dihedral bonds stability/rigidity, e.g. amide bonds rotate much less than H-O-C-C dihedrals
        bool                                                                                 m_WiggleOnly;

        //! atoms to transform if a dihedral bond is mutated. The keys are the various non-ring bonds in fragment
        //! and value is a vector of atom vertices whose position is changed when the bond (key) is mutated
        storage::Vector< size_t>                                                             m_ConnectedAtoms;

        //! molecule dihedral bond index
        storage::Vector< size_t>                                                             m_MoleculeDihedralBondIndex;

        //! original bin
        double m_OriginalAngle;

        //! whether an amide is present on either B or C
        bool m_HasAmide;

        //! whether an amide is present on either B or C
        bool m_IsAmide;

        //! whether this is a double bond (+/- 5 degree deflection is max allowed)
        bool m_IsDoubleBond;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor taking the member variable parameters
      //! @param ATOM_B, ATOM_C atoms involved in the bond
      //! @param GRAPH constitution graph of the molecule of interest
      //! @param MOLECULE molecule of interest
      //! @param FLIP_ONLY -- if true, only rotations of 180 degrees are allowed
      //! @param WIGGLE_ONLY -- if true, do not change dihedral bin, just wiggle the dihedral within its bin
      MutateDihedralBond
      (
        const AtomConformationalInterface &ATOM_B,
        const AtomConformationalInterface &ATOM_C,
        const graph::ConstGraph< size_t, size_t> &GRAPH,
        const FragmentComplete &MOLECULE,
        const bool &FLIP_ONLY,
        const bool &WIGGLE_ONLY
      );

      //! @brief Clone function
      //! @return pointer to new MutateDihedralBond
      MutateDihedralBond *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns atom indices that are connected to different bonds
      //! @return atom indices that are connected to different rotatable bonds
      const storage::Vector< size_t> &GetConnectedAtoms() const;

      //! @brief get the dihedral bond indices referenced by this mutate
      const storage::Vector< size_t> &GetMoleculeDihedralBondIndices() const;

      //! @brief get atom index A of A->B->C->D
      const size_t &GetAtomIndexA() const
      {
        return m_AtomIndexA;
      }

      //! @brief get atom index B of A->B->C->D
      const size_t &GetAtomIndexB() const
      {
        return m_AtomIndexB;
      }

      //! @brief get atom index C of A->B->C->D
      const size_t &GetAtomIndexC() const
      {
        return m_AtomIndexC;
      }

      //! @brief get atom index D of A->B->C->D
      const size_t &GetAtomIndexD() const
      {
        return m_AtomIndexD;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator taking a conformation and returns a mutated conformation
      //! @param MOLECULE conformation of interest
      //! @return MutateResult that results from mutating to the argument
      math::MutateResult< FragmentComplete> Mutate
      (
        const FragmentComplete &MOLECULE
      ) const;

      //! @brief remove the clash of two particular atoms by the most conservative movement of the given atoms if possible
      math::MutateResult< FragmentComplete> RemoveClash
      (
        const FragmentComplete &MOLECULE,
        const size_t &ATOM_A,
        const size_t &ATOM_B,
        const double &MIN_DISTANCE,
        const double &LOWER_BOUND=0.0,
        const double &UPPER_BOUND=400.0
      ) const;

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

    }; // class MutateDihedralBond

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_MUTATE_DIHEDRAL_BOND_H_
