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

#ifndef BCL_CHEMISTRY_PRIORITY_DIHEDRAL_ANGLES_H_
#define BCL_CHEMISTRY_PRIORITY_DIHEDRAL_ANGLES_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_bond_conformational.h"
#include "bcl_chemistry_conformation_interface.h"
#include "graph/bcl_graph_undirected_edge.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PriorityDihedralAngles
    //! @brief Class for determining priority angles of a given molecule.
    //!
    //! @see @link example_chemistry_priority_dihedral_angles.cpp @endlink
    //! @author kothiwsk
    //! @date Jul 05, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PriorityDihedralAngles :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      //! @brief get the wrap around angle; used for determining when to wrap angles in CalculateMinimumDihedralAngle
      //! @return the wrapping angle
      static double &GetChangeableWrappingAngle();

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new PriorityDihedralAngles
      PriorityDihedralAngles *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get the wrap around angle; used for determining when to wrap angles in CalculateMinimumDihedralAngle
      //! @return the wrapping angle
      static const double &GetWrappingAngle();

      //! @brief set the wrap around angle; used for determining when to wrap angles in CalculateMinimumDihedralAngle
      //! @param WRAPPING_ANGLE the wrapping angle
      static void SetWrappingAngle( const double &WRAPPING_ANGLE);

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief calculate prioritized dihedral angle and indices of atoms that make the angle
      //! @param MOLECULE the molecule of interest
      //! @return a pair of vector containing wrapped around dihedral angles in degrees, in ascending order of atom positions
      storage::Pair< storage::Vector< double>, storage::Vector< storage::VectorND< 4, size_t> > >
        operator()( const ConformationInterface &MOLECULE) const;

      //! @brief calculate maximal dihedral angle differences
      //! @param MOLECULE_A, MOLECULE_B the molecules to compare
      //! @return maximal dihedral angle differences
      storage::Vector< double> operator()
      (
        const ConformationInterface &MOLECULE_A,
        const ConformationInterface &MOLECULE_B
      ) const;

      //! @brief calculate the center atoms of dihedral edges
      //! @param MOLECULE the molecule of interest
      //! @return a vector of edges that represent center bonds of dihedral bonds
      static storage::Vector< graph::UndirectedEdge< ConfigurationalBondType> > GetDihedralEdges( const ConformationInterface &MOLECULE);

      //! @brief get the multiplicity of the priorities, to identify symmetric positions
      //! @param MOLECULE the molecule of interest
      static storage::Vector< size_t> CalculateMultiplicity( const ConformationInterface &MOLECULE);

      //! @brief get the priorities, to identify symmetric positions
      //! @param MOLECULE the molecule of interest
      static storage::Vector< size_t> GetPriority( const ConformationInterface &MOLECULE);

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

    ////////////////////
    // helper methods //
    ////////////////////

    private:

      //! @brief calculate priority of all atoms in the molecule and populate the member variable that stores priority
      //! @param MOLECULE the molecule of interest
      static void RecalculatePriority( const ConformationInterface &MOLECULE);

      //! @brief return priority atoms among all atoms connected to a single atom
      //! @param MOLECULE the molecule of interest
      //! @param BONDS all the bonds of the atom for which substituent priority has to be determined
      //! @param DIHEDRAL_BONDED_ATOM atom which needs to be excluded from list of priorities
      //! @return priority atoms among all atoms connected to a single atom
      static const storage::Vector< size_t> DetermineSubstituentPriority
      (
        const ConformationInterface &MOLECULE,
        const linal::Vector< float> &ATOM_PRIORITY,
        const storage::Vector< BondConformational> &BONDS,
        const AtomConformationalInterface &DIHEDRAL_BONDED_ATOM,
        const bool &IS_RING_BOND
      );

      //! @brief calculate prioritized dihedral angle
      //! @param MOLECULE the molecule of interest
      //! @param ATOM_A_PRIORITIES vector of indices of priority atoms connected to atom B
      //! @param ATOM_B atom involved in the central bond of the dihedral bond
      //! @param ATOM_C atom involved in the central bond of the dihedral bond
      //! @param ATOM_D_PRIORITIES vector of indices of priority atoms connected to atom C
      //! @param IS_RING_BOND true if B-C is in a ring
      storage::Pair< double, storage::VectorND< 4, size_t> > CalculateMinimumDihedralAngle
      (
        const ConformationInterface &MOLECULE,
        const storage::Vector< size_t> &ATOM_A_PRIORITIES,
        const AtomConformationalInterface &ATOM_B,
        const AtomConformationalInterface &ATOM_C,
        const storage::Vector< size_t> &ATOM_D_PRIORITIES,
        const bool &IS_RING_BOND
      ) const;

      //! @brief calculate prioritized dihedral angle
      //! @param MOLECULE the molecule of interest
      //! @param ATOM_A_PRIORITIES vector of indices of priority atoms connected to atom B
      //! @param ATOM_B atom involved in the central bond of the dihedral bond
      //! @param ATOM_C atom involved in the central bond of the dihedral bond
      //! @param ATOM_D_PRIORITIES vector of indices of priority atoms connected to atom C
      //! @param IS_RING_BOND true if B-C is in a ring
      storage::Vector< double> CalculateOrderedDihedralAngles
      (
        const ConformationInterface &MOLECULE,
        const storage::Vector< size_t> &ATOM_A_PRIORITIES,
        const AtomConformationalInterface &ATOM_B,
        const AtomConformationalInterface &ATOM_C,
        const storage::Vector< size_t> &ATOM_D_PRIORITIES,
        const bool &IS_RING_BOND
      ) const;

      //! @brief calculate maximum dihedral angle difference, given ordered equivalent sets of dihedral angles
      //! @param DIHEDRALS_A dihedrals about a single bond, ordered by priority, for the first molecule
      //! @param DIHEDRALS_B dihedrals about a single bond, ordered by priority, for the second molecule
      double CalculateMaximumDihedralAngleDifference
      (
        const storage::Vector< double> &DIHEDRALS_A,
        const storage::Vector< double> &DIHEDRALS_B
      ) const;

    }; // class PriorityDihedralAngles
  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_PRIORITY_DIHEDRAL_ANGLES_H_
