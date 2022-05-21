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

#ifndef BCL_CHEMISTRY_ROTAMER_DIHEDRAL_BOND_DATA_H_
#define BCL_CHEMISTRY_ROTAMER_DIHEDRAL_BOND_DATA_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_small_molecule_fragment_isomorphism.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RotamerDihedralBondData
    //! @brief Storage class for storing rotamer information for a single fragment that is contained in the molecule of
    //!        interest
    //! @details
    //!
    //! @see @link example_chemistry_rotamer_dihedral_bond_data.cpp @endlink
    //! @author kothiwsk
    //! @date Jul 10, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RotamerDihedralBondData :
      public util::ObjectInterface
    {

    private:

      //!< pointer to PriorityDihedralAngle object of the molecule that needs to be sampled
      storage::Vector< graph::UndirectedEdge< ConfigurationalBondType> >      m_DihedralEdges;

      //!< pointer to the fragment for which this class stores information
      util::SiPtr< const SmallMoleculeFragmentIsomorphism>                    m_Fragment;

      //!< all the isomorphisms of this fragment with the same set of atoms of the molecule of interest
      storage::Vector< storage::Vector< size_t> >                               m_Isomorphisms;

      //!< priority dihedral angle of fragment in terms of atom count of molecule of interest
      storage::Vector< storage::Vector< storage::VectorND< 4, size_t> > >       m_DihedralAtomIndices;

      //!< edges that represent central bonds of fragment of interest in terms of molecule edges
      storage::Vector< graph::UndirectedEdge< ConfigurationalBondType> >      m_CenterBonds;

      //!< fragment dihedral bonds that are not contained in a ring
      storage::Vector< size_t>                                                m_NonRingBonds;

      //!< mapping between molecule dihedral bond and fragment dihedral bond
      storage::Vector< storage::Vector< size_t> >                             m_CenterBondIsomorphisms;

      //! Weights associated with each isomorphism
      linal::Vector< double>                                                  m_IsomorphismWeights;

      //!< whether fragment contains rings
      bool                                                                    m_ContainsRings;

      //!< whether fragment contains different conformations of rings
      bool                                                                    m_HasDifferentRingConformations;

      //! Whether different ring rotamers exist (e.g. due to different isomorphisms)
      bool                                                                    m_HasDifferentRingRotamers;

      //! Whether fragment has incomplete rings
      bool                                                                    m_IncompleteRings;

      //! Minimum number of times any central bond of this fragment is represented in the fragment ensemble
      double                                                                  m_Uniqueness;

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
      RotamerDihedralBondData()
      {
      }

      //! @brief constructor given molecule and dihedral angle map
      //! @param MOLECULE_PRIORITY PriorityDihedralAngles object of the molecule whose conformations need to be sampled
      //! @param FRAGMENT the fragment for which this class stores the dihedral bond information
      //! @param ISOMORPHISMS all possible of isomorphisms for the same set of atoms in the parent molecule
      RotamerDihedralBondData
      (
        const ConformationInterface &MOLECULE,
        const SmallMoleculeFragmentIsomorphism &FRAGMENT,
        const storage::List< storage::Vector< size_t> > &ISOMORPHISMS
      );

      //! @brief Clone function
      //! @return pointer to new RotamerDihedralBondData
      RotamerDihedralBondData *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return center bonds of the molecule of interest
      //! @return center bonds of the molecule of interest
      const storage::Vector< graph::UndirectedEdge< ConfigurationalBondType> > &GetMoleculeCenterBonds() const;

      //! @brief return fragment for which this class stores dihedral bond information
      //! @return fragment for which this class stores dihedral bond information
      const SmallMoleculeFragmentIsomorphism &GetFragment() const;

      //! @brief returns the total number of times the rotatable element has been seen in the CSD
      //! @return the total number of times the rotatable element is seen in the CSD
      const double &GetFragmentCounts() const;

      //! @brief returns rotamers of the fragment whose dihedral information is stored by this class
      //! @return the rotamers of the fragment whose dihedral information is stored by this class
      const storage::Vector< linal::Vector< double> > &GetRotamers() const;

      //! @brief returns counts for each rotamer of the fragment
      //! @return counts for each rotamer
      const storage::Vector< double> &GetRotamerCounts() const;

      //! @brief returns all isomorphism of this fragment for the same set of atoms in the parent molecule
      //! @return isomorphism of this fragment for the same set of atoms in the parent molecule
      const storage::Vector< storage::Vector< size_t> > &GetIsomorphisms() const;

      //! @brief returns dihedral bonds of the molecule of interest that this fragment isomorphism represents
      //! @return dihedral bonds of the molecule of interest that this fragment represents
      const storage::Vector< storage::Vector< storage::VectorND< 4, size_t> > > &GetRotamerBonds() const;

      //! @brief returns center bonds that this fragment represents in the molecule of interest
      //! @return center bonds of that this fragment represents in the molecule of interest
      //! @note This specifically corresponds to the first isomorphism -- others may have a different central bond isomorphism
      //!       that covers the same set of bonds, but potentially in different orders
      const storage::Vector< graph::UndirectedEdge< ConfigurationalBondType> > &GetCenterBonds() const;

      //! @brief returns weights for each isomorphism
      const linal::Vector< double> &GetIsomorphismWeights() const;

      //! @brief returns the MoleculeCenterBonds that are not contained in rings
      //! @return the MoleculeCenterBonds that are not contained in rings
      //! @note This specifically corresponds to the first isomorphism -- others may have a different central bond isomorphism
      //!       that covers the same set of bonds, but potentially in different orders
      const storage::Vector< size_t> &GetNonRingBonds() const;

      //! @brief returns isomorphism between fragment center bond and molecule bond
      //! @return isomorphism between fragment center bond and molecule bond
      //! @note This specifically corresponds to the first isomorphism -- others may have a different central bond isomorphism
      //!       that covers the same set of bonds, but potentially in different orders
      const storage::Vector< size_t> &GetCenterBondIsomorphism() const;

      //! @brief returns isomorphism between fragment center bond and molecule bond
      //! @return isomorphism between fragment center bond and molecule bond
      const storage::Vector< storage::Vector< size_t> > &GetCenterBondIsomorphisms() const;

      //! @brief returns true if fragment contains rings
      //! @return true if fragment contains ring otherwise false
      bool ContainsRings() const;

      //! @brief returns true if fragment can be used to create rings with different conformations of the overall molecule
      //! @return true if fragment can be used to create rings with different conformations of the overall molecule
      bool ContainsRingRotamers() const;

      //! @brief returns true if fragment contains rings with different conformations
      //! @return true if fragment contains rings with different conformations else false
      bool HasIncompleteRings() const
      {
        return m_IncompleteRings;
      }

      //! @brief returns the bin size strategy used for representing conformation of fragment
      //! @return the bin size strategy used for representing conformation of fragment
      const double &GetBinSize() const;

      //! @brief determine the imputed counts of a given rotamer
      //! @return imputed counts of a given rotamer
      double GetImputedRotamerCounts( const linal::Vector< int> &BINS) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief check whether elements are fully contained in the container of interest
      //! @param CONTAINER container which needs to be checked to see if it contains elements
      //! @param ELEMENTS elements that need to be checked to see if they are fully contained in the CONTAINER
      //! @return true if elements completely contained in the container otherwise false.
      static bool IfFullIntersection
      (
        const storage::Vector< size_t> &CONTAINER,
        const linal::Vector< size_t> &ELEMENTS
      );

    ///////////////
    // operators //
    ///////////////

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

      //! @brief determine atoms involved in rotatable bonds and atoms involved in center bonds
      //! @return void
      void CalculateBondData();

      //! @brief determine rings and bonds contained in the fragment
      //! @return void
      void DetermineBondsAndRings();

    }; // class RotamerDihedralBondData

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_ROTAMER_DIHEDRAL_BOND_DATA_H_
