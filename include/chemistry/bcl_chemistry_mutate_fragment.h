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

#ifndef BCL_CHEMISTRY_MUTATE_FRAGMENT_H_
#define BCL_CHEMISTRY_MUTATE_FRAGMENT_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_bond_angle_assignment.h"
#include "bcl_chemistry_conformation_interface.h"
#include "bcl_chemistry_mutate_bond_angles.h"
#include "bcl_chemistry_mutate_chirality.h"
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
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateFragment
    //! @brief This class is for randomly using different rotamers of a fragment to change conformation of a molecule of interest.
    //! @details
    //!
    //! @see @link example_chemistry_mutate_fragment.cpp @endlink
    //! @author kothiwsk
    //! @date Sep 04, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateFragment :
      public MutateDihedralsInterface
    {

    private:

    //////////
    // data //
    //////////

      //! rotamer data for fragment which this class mutates
      util::ShPtr< RotamerDihedralBondData>                                               m_FragmentData;

      //! probability distribution for different rotamers of fragemnt
      mutable util::ShPtr< random::Histogram1DDistribution>                                m_ProbabilityDistribution;

      //! atoms to transform if a dihedral bond is mutated. The keys are the various non-ring bonds in fragment
      //! and value is a vector of atom vertices whose position is changed when the bond (key) is mutated
      storage::Map< size_t, storage::Vector< size_t> >                                      m_ConnectedAtoms;

      //! graph complement of fragment rings that exists in multiple conformations
      mutable storage::Vector< storage::Map< size_t, storage::Vector< RingFragmentMap> > > m_FragmentMap;

      //! if fragment ring conformations have to be applied to the molecule of interest
      //! Mutable for when chirality or isometry are constrained, since it is possible that all of the ring conformations
      //! yield incorrect chirality, in which case we shouldn't update the rings
      mutable bool                                                                          m_UpdateRings;

      size_t                                                                                m_MaxRotamer;

      bool                                                                                  m_ChangeChirality;

      bool                                                                                  m_ChangeIsometry;

      //! Whether to perform unbiased sampling
      bool                                                                                  m_Unbiased;

      //! Whether to sample the exact dihedral angle from a distribution
      bool                                                                                  m_SampleDihedralAngles;

      //! original chirality string
      std::string                                                                           m_ChiralityString;
      std::string                                                                           m_IsometryString;

      //! disallowed pairings of chosen isometry and fragment
      mutable storage::Set< storage::Pair< size_t, size_t> >                                m_DisallowedPairs;

      //! store the last chosen isometry and fragment
      mutable storage::Pair< size_t, size_t>                                                m_LastChosen;

      // Mutate chirality. Only saved if change chirality or isometry is set and molecule is chiral
      util::SiPtr< const MutateChirality>                                                   m_MutateChirality;

      //
      storage::Vector< size_t>                                                              m_DihedralIndicesMutated;

      // whether to check amide bond planarity
      bool                                                                                  m_CheckAmideBondPlanarity;

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
      MutateFragment();

      //! @brief constructor taking the member variable parameters
      //! @param FRAGMENT fragment which the mutate mutates
      //! @param GRAPH constitution graph of the molecule of interest
      //! @param MOLECULE molecule of interest
      MutateFragment
      (
        const util::ShPtr< RotamerDihedralBondData> &FRAGMENT,
        const graph::ConstGraph< size_t, size_t> &GRAPH,
        const FragmentComplete &MOLECULE,
        const MutateChirality &MUTATE_CHIRALITY,
        bool CHANGE_CHIRALITY = false,
        bool UNBIASED = false,
        bool SAMPLE_DIHEDRALS = false,
        bool CHECK_AMIDE_PLANARITY = true
      );

      //! @brief Clone function
      //! @return pointer to new MutateFragment
      MutateFragment *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetScheme() const;

      //! @brief returns the fragment that the mutate uses
      //! @return the fragment that the mutate uses
      const util::ShPtr< RotamerDihedralBondData> &GetFragmentData() const;

      //! @brief returns atom indices that are connected to different bonds
      //! @return atom indices that are connected to different rotatable bonds
      const storage::Map< size_t, storage::Vector< size_t> > &GetConnectedAtoms() const;

      //! @brief returns probability distribution for rotamers of this particular fragments
      //! @return probability distribution for rotamers of this particular fragments
      const random::Histogram1DDistribution &GetProbability() const;

      //! @brief get the dihedral bond indices referenced by this mutate
      const storage::Vector< size_t> &GetMoleculeDihedralBondIndices() const;

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

      //! @brief rotate the molecule of interest at a particular bond at a particular angle
      //! @param ATOM_INFO the atom info that needs to be updated with coordinates after rotating about a particular bond
      //! @param CONNECTED_ATOMS the connected atoms associated with bond of interest
      //! @param NEW_ANGLE the new angle for the bond of interest
      //! @param EXISTING_ANGLE the angle for the bond of interest in the argument that is passed in
      static void RotateBond
      (
        storage::Vector< sdf::AtomInfo> &ATOM_INFO,
        const storage::Vector< size_t> &CONNECTED_ATOMS,
        double NEW_ANGLE,
        double EXISTING_ANGLE
      );

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
    /////////////////////

    public:

      //! @brief return a vector of atoms that are connected to a bond of interest
      //! @param ATOM_INDICES atom indices of the center bond around which connected bonds have to be found
      //! @param GRAPH constitution graph of the molecule of interest
      //! @param return a vector of atoms that are connected to a bond of interest
      static storage::Vector< size_t> GetConnectedAtoms
      (
        const graph::UndirectedEdge< ConfigurationalBondType> &ATOM_INDICES,
        const graph::ConstGraph< size_t, size_t> &GRAPH
      );

      //! @brief get a vector containing angle for each rotatable bond in the molecule
      //! @param ATOM_INFO the atom info that needs to be updated with coordinates after rotatating about a particular bond
      //! @param DIHEDRAL_BONDS atom indices making up the different rotatable bonds in the molecule
      //! @return vector containing angle for each rotatable bond in the molecule
      static storage::Vector< double> GetDihedralAngles
      (
        const storage::Vector< sdf::AtomInfo> &ATOM_INFO,
        const storage::Vector< storage::VectorND< 4, size_t> > &DIHEDRAL_BONDS
      );

      //! @brief get a vector containing angle for each rotatable bond in the molecule
      //! @param ATOM_INFO the atom info of molecule of interest
      //! @param DIHEDRAL_BONDS atom indices making up the different rotatable bonds in the molecule
      //! @return vector containing angle for each rotatable bond in the molecule
      static double GetDihedralAngle
      (
        const storage::Vector< sdf::AtomInfo> &ATOM_INFO,
        const storage::VectorND< 4, size_t> &DIHEDRAL_BOND
      );

    private:

      //! @brief get complement of fragment in the molecule and information about bonds that can connect fragment to other parts of the molecule
      //! @param GRAPH the constitution graph of the molecule on which this mutates operates
      //! @return
      storage::Vector< storage::Map< size_t, storage::Vector< RingFragmentMap> > > GetFragmentLinkWithComplement
      (
        const graph::ConstGraph< size_t, size_t> &GRAPH
      ) const;

      //! @brief get information about parts of molecule that are complementary to the fragment of interest
      //! @param ATOM_INFO atom info of the new conformation
      //! @param ISOMORPHISM isomorphism between fragment and molecule of interest
      //! @param MOLECULE_ATOM_INFO atom info of conformation that is given
      //! @param FRAGMENT_ATOM_INFO atom info of fragment that this mutate uses
      //! @param MAPPING_INFORMATION connectivity between fragment and rest of the molecule
      //! @return return informatin about parts of molecule that are complementary to the fragment of interest
      void MappingFragment
      (
        storage::Vector< sdf::AtomInfo> &ATOM_INFO,
        const storage::Vector< size_t> &ISOMORPHISM,
        const FragmentComplete &MOLECULE,
        const AtomVector< AtomComplete> &FRAGMENT,
        const std::pair< size_t, storage::Vector< RingFragmentMap> > &MAPPING_INFORMATION
      ) const;

      //! @brief get information about parts of molecule that are complementary to the fragment of interest
      //! @param ATOM_INFO atom info of the new conformation
      //! @param ISOMORPHISM isomorphism between fragment and molecule of interest
      //! @param MOLECULE_ATOM_INFO atom info of conformation that is given
      //! @param FRAGMENT_ATOM_INFO atom info of fragment that this mutate uses
      //! @param MAPPING_INFORMATION connectivity between fragment and rest of the molecule
      //! @return return informatin about parts of molecule that are complementary to the fragment of interest
      void FixDihedral
      (
        storage::Vector< sdf::AtomInfo> &ATOM_INFO,
        const FragmentComplete &MOLECULE,
        const std::pair< size_t, storage::Vector< RingFragmentMap> > &MAPPING_INFORMATION
      ) const;

    }; // class MutateFragment

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_MUTATE_FRAGMENT_H_
