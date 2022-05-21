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

#ifndef BCL_CHEMISTRY_BOND_ANGLE_ASSIGNMENT_H_
#define BCL_CHEMISTRY_BOND_ANGLE_ASSIGNMENT_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_fragment_complete.h"
#include "bcl_chemistry_priority_dihedral_angles.h"
#include "bcl_chemistry_rotamer_library_interface.h"
#include "coord/bcl_coord_point_cloud.h"
#include "graph/bcl_graph_subgraph_isomorphism.h"
#include "math/bcl_math_running_average_sd.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_si_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BondAngleAssignment
    //! @brief data class that stores isomorphism object between Scaffold(fragment) and a Molecule of interest.
    //! @details stores isomorphism object between a Scaffold and the Molecule, and scaffold
    //!
    //!
    //! @see @link example_chemistry_bond_angle_assignment.cpp @endlink
    //! @author mendenjl
    //! @date Jun 19, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API BondAngleAssignment :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      //!< Bond angles with counts
      util::SiPtr< const RotamerLibraryInterface::t_BondAnglesWithCounts> m_BondAngleCounts;

      //! storage for m_BondAngleCounts if not considering alternate chiralities and a chirality is present
      RotamerLibraryInterface::t_BondAnglesWithCounts m_BondAngleCountsStorage;

      //!< atom index in the corresponding molecule
      size_t                   m_CentralAtomIndex;

      //! use m_AttachedAtomIndices to go from bond-of-atom index to row-of-matrix index
      storage::Vector< size_t> m_AttachedAtomIndices;

      //! use m_AttachedAtomInverseIndices to go from row-of-bond angle matrix to bond-of-atom index
      storage::Vector< size_t> m_AttachedAtomInverseIndices;

      //! Same as m_AttachedAtomInverseIndices except has raw atom indices rather than bond indices
      storage::Vector< size_t> m_AttachedAtomInverseRawAtomIndices;

      //! Free energy for each rotamer
      linal::Vector< float>    m_BondAngleFreeEnergies;

      //! free energy of unseen rotamer
      float                    m_FreeEnergyUnseenRotamer;

      //! expectation value
      math::RunningAverageSD< float> m_ExpectedBondAngleFreeEnergies;

      //! whether alternate chiralities are to be considered
      bool m_ConsiderAlternateChiralities;

      //! Number of ring bonds
      size_t m_NumberRingBonds;

      //! Whether the atoms nominal chirality can be flipped (primarily for nitrogen or tetrahedral boron/aluminum)
      bool m_CanFlip;

      //! whether the atom is connected by an amide bond to any other atom
      bool m_HasAmide;

      //! whether to trust the bond lengths (false if search was not specific for atom type)
      bool m_TrustBondLengths;

    public:

      //! @brief get the pseudocount that is used for all possible rotamers (default 10)
      //! @return the pseudocount that is used for all possible rotamers
      static double &GetPseudocount();

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      BondAngleAssignment();

      //! @brief Constructor from data members
      //! @param FRAGMENT fragment whose rotamer information this class stores
      //! @param FRAGMENTCSI isomorphism object of fragment and molecule
      BondAngleAssignment
      (
        const typename RotamerLibraryInterface::t_BondAnglesWithCounts &BOND_ANGLE_COUNTS,
        const typename RotamerLibraryInterface::t_BondAngleMapKey &BOND_ANGLE_MAP_KEY,
        const ConformationInterface &MOLECULE,
        const size_t &CENTRAL_ATOM_INDICES,
        const bool &CONSIDER_ALTERNATE_CHIRALITIES = false,
        const ConformationGraphConverter::AtomComparisonType &ATOM_LEVEL = ConformationGraphConverter::e_AtomType
      );

      //! @brief Clone function
      //! @return pointer to new BondAngleAssignment
      BondAngleAssignment *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the fragment
      //! @return fragment
      const size_t &GetCentralAtomIndex() const
      {
        return m_CentralAtomIndex;
      }

      //! @brief returns the isomorphism object of molecule and scaffold
      //! @return isomorphism object of molecule and scaffold
      const linal::Vector< float> &GetBondAngleFreeEnergies() const
      {
        return m_BondAngleFreeEnergies;
      }

      //! @brief get whether we can use the bond lengths / whether they are reliable
      const bool &GetCanUseBondLengths() const
      {
        return m_TrustBondLengths;
      }

      //! @brief returns the isomorphism object of molecule and scaffold
      //! @return isomorphism object of molecule and scaffold
      const RotamerLibraryInterface::t_BondAnglesWithCounts &GetBondAngleCounts() const
      {
        return *m_BondAngleCounts;
      }

      //! @brief returns the isomorphism object of molecule and scaffold
      //! @return isomorphism object of molecule and scaffold
      const storage::Vector< size_t> &GetAttachedAtomIndices() const
      {
        return m_AttachedAtomIndices;
      }

      //! @brief returns the isomorphism object of molecule and scaffold
      //! @return isomorphism object of molecule and scaffold
      const storage::Vector< size_t> &GetAttachedAtomIndicesInverse() const
      {
        return m_AttachedAtomInverseIndices;
      }

      //! @brief returns the isomorphism object of molecule and scaffold
      //! @return isomorphism object of molecule and scaffold
      const size_t &GetNumberRingBonds() const
      {
        return m_NumberRingBonds;
      }

      //! @brief return bin size used for binning of dihedral angles
      //! @return the bin size used for binning of dihedral angles
      const bool IsDefined() const
      {
        return m_BondAngleCounts.IsDefined();
      }

      //! @brief get the theoretical maximum # of rotamers
      double GetBestFreeEnergy() const
      {
        return m_BondAngleFreeEnergies.Min();
      }

      //! @brief Get the expected free energy
      const math::RunningAverageSD< float> &GetExpectedFreeEnergy() const
      {
        return m_ExpectedBondAngleFreeEnergies;
      }

      //! @brief get whether the atom can flip
      const bool &CanFlip() const
      {
        return m_CanFlip;
      }

      //! @brief get whether the atom has an amide bond
      const bool &HasAnAmideBond() const
      {
        return m_HasAmide;
      }

      //! @brief get the free energy of a particular rotamer
      //! @param ROTAMER rotamer signature of interest
      //! @param CONSIDER_ISOMETRY_CHANGES true - consider changes in isometry
      double GetFreeEnergy( const FragmentComplete &FRAGMENT) const;

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get a standardized coordinate matrix for a given atom
      //! @param ATOM the atom of interest
      static linal::Matrix< double> GetStandardizedCoordinateMatrix
      (
        const linal::Vector3D &ATOM_POS,
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const size_t &FULL_N_BONDS,
        const size_t &N_RING_BONDS
      );

      //! @brief get a standardized coordinate matrix for a given atom
      //! @param ATOM the atom of interest
      linal::Matrix< double> GetStandardizedCoordinateMatrix
      (
        const ConformationInterface &MOLECULE
      ) const;

      //! @brief Transfer standardized coordinate matrix (usually taken from a different molecule) to another
      //! @param MOLECULE the molecule of interest
      //! @param STANDARDIZED_COORDS the coordinates of interest
      void TransferStandardizedCoordinateMatrixForRingSystem
      (
        const storage::Vector< sdf::AtomInfo> &MOLECULE,
        linal::Matrix< double> &STANDARDIZED_COORDS
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
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class BondAngleAssignment

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_BOND_ANGLE_ASSIGNMENT_H_

