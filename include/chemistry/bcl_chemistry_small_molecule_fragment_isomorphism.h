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

#ifndef BCL_CHEMISTRY_SMALL_MOLECULE_FRAGMENT_ISOMORPHISM_H_
#define BCL_CHEMISTRY_SMALL_MOLECULE_FRAGMENT_ISOMORPHISM_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_fragment_complete.h"
#include "bcl_chemistry_priority_dihedral_angles.h"
#include "coord/bcl_coord_point_cloud.h"
#include "graph/bcl_graph_subgraph_isomorphism.h"
#include "linal/bcl_linal_vector_operations.h"
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
    //! @class SmallMoleculeFragmentIsomorphism
    //! @brief data class that stores isomorphism object between Scaffold(fragment) and a Molecule of interest.
    //! @details stores isomorphism object between a Scaffold and the Molecule, and scaffold
    //!
    //!
    //! @see @link example_chemistry_small_molecule_fragment_isomorphism.cpp @endlink
    //! @author kothiwsk
    //! @date Aug 26, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SmallMoleculeFragmentIsomorphism :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      //!< pointer to the molecule for which this object stores fragment information
      util::ShPtr< FragmentComplete>                                    m_Molecule;

      //!< fragment of interest for which this class stores rotamer information
      FragmentComplete                                                  m_Fragment;

      //!< isomorphism object of molecule and fragment
      storage::Vector< storage::Vector< size_t> >                       m_FragmentCSI;

      //!< rotmaer information or dihedral angle information for the fragment of interest
      storage::Vector< linal::Vector< double> >                         m_Rotamers;

      //!< counts for each rotamer
      storage::Vector< double>                                          m_RotamerCounts;

      //!< rotamer bin signature for each rotamer of the fragment
      storage::Vector< linal::Vector< int> >                            m_RotamerBins;

      //! Index map for rotamers, to enable faster lookup
      storage::Map< linal::Vector< int>, size_t>                        m_RotamerIndexMap;

      //! Free energy for each rotamer
      storage::Vector< double>                                          m_RotamerFreeEnergy;

      //! Free energy for each rotamer, assuming isometry changes are forbidden
      storage::Vector< double>                                          m_RotamerFreeEnergyIsometryChangesForbidden;

      //!< total number of rotamers for the fragment of interest
      size_t                                                            m_TotalRotamers;

      //!< total number of rotatable bonds in the fragment of interest
      size_t                                                            m_RotatableBondCount;

      //!< total number of rotatable bonds in the fragment of interest
      size_t                                                            m_RingBondCount;

      //!< total number of rotatable bonds in the fragment of interest
      size_t                                                            m_ChainBondCount;

      //!< total number of rotatable bonds in the fragment of interest
      size_t                                                            m_DihedralChainBondCount;

      //!< Maximum counts of any particular rotamer
      double                                                            m_MaxRotamerCounts;

      //!< coordinates of the cluster center of each rotamer of fragment of interest
      storage::Vector< storage::Vector< linal::Vector3D> >              m_RotamerCoordinates;

      //!< dihedral angle edges
      storage::Vector< storage::VectorND< 4, size_t> >                  m_DihedralAngleIndices;

      //!< dihedral bonds
      storage::Vector< graph::UndirectedEdge< ConfigurationalBondType> > m_DihedralBonds;

      //!< number of times fragment is seen in structure database
      double                                                            m_TotalFragmentCounts;

      //!< bin size used for binning of dihedral angles for the fragment
      double                                                            m_BinSize;

      //!< contains rings
      bool                                                              m_ContainsRings;

      //!< contains different ring conformations
      bool                                                              m_ContainsRingConformations;

      //! Maximum possible different rotamer combinations
      double                                                            m_TheoreticalMaxRotamers;

      //! Maximum possible different rotamer combinations if isometry is ignored
      double                                                            m_TheoreticalMaxRotamersNoIsometry;

      //! Free energy of unseen rotamer
      double                                                            m_FreeEnergyUnseenRotamer;

      //! Free energy of unseen rotamer
      double                                                            m_FreeEnergyUnseenRotamerNoIsometry;

      //! Free energy of best rotamer
      double                                                            m_FreeEnergyBestRotamer;

      //! Free energy of best rotamer
      double                                                            m_FreeEnergyBestRotamerNoIsometry;

      //! Average free energy with isometry
      math::RunningAverageSD< double>                                   m_ExpectedFreeEnergy;

      //! Average free energy without isometry
      math::RunningAverageSD< double>                                   m_ExpectedFreeEnergyNoIsometry;

      //! whether the mapping contains incomplete ring systems
      linal::Vector< int>                                               m_ContainsIncompleteRings;

      //! whether the mapping contains incomplete hydrogenation for non-terminal, non-ring atoms
      linal::Vector< int>                                               m_ContainsIncompleteHydrogenation;

      //! whether the small molecule fragment isomorphism only consists of incomplete hydrogenation
      bool                                                              m_AllIncompleteHydrogenation;

      //! Number of ring systems
      size_t                                                            m_NumberRingSystems;

      //! Number of ring systems
      size_t                                                            m_NumberAromaticRingSystems;

      //! Number of chain systems
      size_t                                                            m_NumberChainSystems;

      //! number of simulated states that haven't been observed
      double                                                            m_FractionSimulatedNotObserved;

      //! simulated count
      double                                                            m_SimulatedCount;

      //! Number simulated for each state
      linal::Vector< float>                                             m_RotamerFractionSimulated;

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
      SmallMoleculeFragmentIsomorphism();

      //! @brief Constructor from data members
      //! @param FRAGMENT fragment whose rotamer information this class stores
      //! @param FRAGMENTCSI isomorphism object of fragment and molecule
      SmallMoleculeFragmentIsomorphism
      (
        const FragmentComplete &MOLECULE,
        const FragmentComplete &FRAGMENT,
        const graph::SubgraphIsomorphism< size_t, size_t> &FRAGMENT_CSI,
        bool REMOVE_RING_ROTAMER_INFO_UNLESS_FRAGMENT_IS_JUST_A_RING = true
      );

      //! @brief Clone function
      //! @return pointer to new SmallMoleculeFragmentIsomorphism
      SmallMoleculeFragmentIsomorphism *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the fragment
      //! @return fragment
      const FragmentComplete &GetMolecule() const
      {
        return *m_Molecule;
      }

      //! @brief returns the fragment
      //! @return fragment
      const FragmentComplete &GetFragment() const
      {
        return m_Fragment;
      }

      //! @brief returns the isomorphism object of molecule and scaffold
      //! @return isomorphism object of molecule and scaffold
      const storage::Vector< storage::Vector< size_t> > &GetFragmentCSI() const
      {
        return m_FragmentCSI;
      }

      //! @brief returns dihedral bond angles of rotamers of the fragment of interest
      //! @return dihedral bond angles of rotamers of the fragment of interest
      const storage::Vector< linal::Vector< double> > &GetRotamers() const
      {
        return m_Rotamers;
      }

      //! @brief returns the counts for each rotamer of fragment of interest
      //! @return the counts for each rotamer of fragment of interest
      const size_t &GetNumberRingSystems() const
      {
        return m_NumberRingSystems;
      }

      //! @brief returns the counts for each rotamer of fragment of interest
      //! @return the counts for each rotamer of fragment of interest
      const size_t &GetNumberAromaticRingSystems() const
      {
        return m_NumberAromaticRingSystems;
      }

      //! @brief returns the counts for each rotamer of fragment of interest
      //! @return the counts for each rotamer of fragment of interest
      const size_t &GetNumberChainSystems() const
      {
        return m_NumberChainSystems;
      }

      //! @brief returns the counts for each rotamer of fragment of interest
      //! @return the counts for each rotamer of fragment of interest
      const storage::Vector< double> &GetRotamerCounts() const
      {
        return m_RotamerCounts;
      }

      //! @brief returns the counts for each rotamer of fragment of interest
      //! @return the counts for each rotamer of fragment of interest
      const double &GetMaxRotamerCount() const
      {
        return m_MaxRotamerCounts;
      }

      //! @brief returns dihedral bins of rotamers of the fragment of interest
      //! @return dihedral bins of rotamers of the fragment of interest
      const storage::Vector< linal::Vector< int> > &GetRotamerBins() const
      {
        return m_RotamerBins;
      }

      //! @brief rerturn number of rotamers of the fragment of interest
      //! @return number of rotamers of the fragment of interest
      const size_t &GetRotamerNumbers() const
      {
        return m_TotalRotamers;
      }

      //! @brief return number of rotatable dihedral bonds of the fragment of interest
      //! @return number of rotatable dihedral bonds of the fragment of interest
      const size_t &GetRotatableBondNumber() const
      {
        return m_RotatableBondCount;
      }

      //! @brief return number of rotatable dihedral bonds of the fragment of interest
      //! @return number of rotatable dihedral bonds of the fragment of interest
      const size_t &GetNumberRingBonds() const
      {
        return m_RingBondCount;
      }

      //! @brief return number of rotatable dihedral bonds of the fragment of interest
      //! @return number of rotatable dihedral bonds of the fragment of interest
      const size_t &GetNumberChainBonds() const
      {
        return m_ChainBondCount;
      }

      //! @brief return number of rotatable dihedral bonds of the fragment of interest
      //! @return number of rotatable dihedral bonds of the fragment of interest
      size_t GetNumberDihedralChainBonds() const
      {
        return m_DihedralChainBondCount;
      }

      //! @brief return the number of atoms in the fragment
      //! @return the number of atoms in the fragment
      const size_t &GetFragmentSize() const
      {
        return m_Fragment.GetNumberAtoms();
      }

      //! @brief return the fragment count in structure database
      //! @return the fragment count in structure database
      const double &GetFragmentCounts() const
      {
        return m_TotalFragmentCounts;
      }

      //! @brief return the rotamer coordinates for fragment of interest
      //! @return the rotamer coordinates for fragment of interest
      const storage::Vector< storage::Vector< linal::Vector3D> > &GetRotamerCoordinates() const
      {
        return m_RotamerCoordinates;
      }

      //! @brief return the dihedral angle indices
      const storage::Vector< storage::VectorND< 4, size_t> > &GetDihedralAngleIndices() const
      {
        return m_DihedralAngleIndices;
      }

      //! @brief return the dihedral angle indices
      const storage::Vector< graph::UndirectedEdge< ConfigurationalBondType> > &GetDihedralEdges() const
      {
        return m_DihedralBonds;
      }

      //! @brief return bin size used for binning of dihedral angles
      //! @return the bin size used for binning of dihedral angles
      const double &GetBinSize() const
      {
        return m_BinSize;
      }

      //! @brief return bin size used for binning of dihedral angles
      //! @return the bin size used for binning of dihedral angles
      const bool IsDefined() const
      {
        return m_Fragment.GetNumberAtoms();
      }

      //! @brief return bin size used for binning of dihedral angles
      //! @return the bin size used for binning of dihedral angles
      const bool ContainsRings() const
      {
        return m_ContainsRings;
      }

      //! @brief return bin size used for binning of dihedral angles
      //! @return the bin size used for binning of dihedral angles
      const bool ContainsRingConformations() const
      {
        return m_ContainsRingConformations;
      }

      //! @brief return true if the isomorphism has incomplete rings
      //! @return true if there are ring systems in the molecule that are incompletely represented in the fragment
      const int &ContainsIncompleteRings( const size_t &ISOMORPHISM_NUMBER) const
      {
        return m_ContainsIncompleteRings( ISOMORPHISM_NUMBER);
      }

      //! @brief return true if the isomorphism has incomplete hydrogenation, except at rings and otherwise terminal heavy atoms
      //! @return true if the isomorphism has incomplete hydrogenation, except at rings and otherwise terminal heavy atoms
      const int &ContainsIncompleteHydrogenation( const size_t &ISOMORPHISM_NUMBER) const
      {
        return m_ContainsIncompleteHydrogenation( ISOMORPHISM_NUMBER);
      }

      //! @brief return true if all isomorphisms have incomplete hydrogenation, except at rings and otherwise terminal heavy atoms
      //! @return true if all isomorphism has incomplete hydrogenation, except at rings and otherwise terminal heavy atoms
      const bool &ContainsOnlyIncompleteHydrogenation() const
      {
        return m_AllIncompleteHydrogenation;
      }

      //! @brief get the theoretical maximum # of rotamers
      double GetBestFreeEnergy() const
      {
        return m_FreeEnergyBestRotamer;
      }

      //! @brief get the theoretical maximum # of rotamers assuming isometry changes are forbidden
      double GetBestFreeEnergyNoIsometry() const
      {
        return m_FreeEnergyBestRotamerNoIsometry;
      }

      //! @brief Get the expected free energy
      const math::RunningAverageSD< double> &GetExpectedFreeEnergy( const bool &CONSIDER_ISOMETRY) const
      {
        return CONSIDER_ISOMETRY ? m_ExpectedFreeEnergy : m_ExpectedFreeEnergyNoIsometry;
      }

      //! @brief get the free energy of a particular rotamer
      //! @param ROTAMER rotamer signature of interest
      //! @param CONSIDER_ISOMETRY_CHANGES true - consider changes in isometry
      double GetRotamerFreeEnergy( const linal::Vector< int> &ROTAMER, const bool &CONSIDER_ISOMETRY_CHANGES) const;

      //! @brief get the free energy of a particular rotamer relative to others with the same dihedral bin
      //! @param ROTAMER rotamer signature of interest
      double GetRotamerDihedralFineEnergy
      (
        const linal::Vector< double> &ROTAMER_ANGLES,
        const linal::Vector< int> &ROTAMER
      ) const;

      //! @brief Get estimated maximum number of rotamers
      double GetEstimatedMaxRotamers() const
      {
        return m_TheoreticalMaxRotamers;
      }

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

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

    private:

      //! @brief get all stored properties of the fragment of interest
      //! @return rotamer data of the fragment
      void RetrieveStoredProperties( const bool &REMOVE_RING_ROTAMER_INFO_UNLESS_FRAGMENT_IS_JUST_A_RING);

    }; // class SmallMoleculeFragmentIsomorphism

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_SMALL_MOLECULE_FRAGMENT_ISOMORPHISM_H_

