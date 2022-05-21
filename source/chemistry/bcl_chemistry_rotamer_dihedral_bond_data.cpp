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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "chemistry/bcl_chemistry_rotamer_dihedral_bond_data.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_bond_dihedral_steric_weight.h"
#include "chemistry/bcl_chemistry_rotamer_cluster_center.h"
#include "chemistry/bcl_chemistry_small_molecule_fragment_mapping.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> RotamerDihedralBondData::s_Instance
    (
      GetObjectInstances().AddInstance( new RotamerDihedralBondData())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor given molecule and dihedral angle map
    //! @param MOLECULE_PRIORITY PriorityDihedralAngles object of the molecule whose conformations need to be sampled
    //! @param FRAGMENT the fragment for which this class stores the dihedral bond information
    //! @param ISOMORPHISMS all possible of isomorphisms for the same set of atoms in the parent molecule
    RotamerDihedralBondData::RotamerDihedralBondData
    (
      const ConformationInterface &MOLECULE,
      const SmallMoleculeFragmentIsomorphism &FRAGMENT,
      const storage::List< storage::Vector< size_t> > &ISOMORPHISMS
    ) :
      m_DihedralEdges( PriorityDihedralAngles::GetDihedralEdges( MOLECULE)),
      m_Fragment( FRAGMENT),
      m_Isomorphisms( ISOMORPHISMS.Begin(), ISOMORPHISMS.End()),
      m_ContainsRings( FRAGMENT.ContainsRings()),
      m_HasDifferentRingRotamers( FRAGMENT.ContainsRingConformations()),
      m_Uniqueness( 1.0)
    {
      if( FRAGMENT.GetFragment().IsPropertyStored( "RotamerCounts"))
      {
//        if( FRAGMENT.GetFragmentSize() == size_t( 4))
//        {
//          FRAGMENT.GetFragment().WriteMDL( util::GetLogger());
//        }
//        BCL_Debug( ISOMORPHISMS);
        CalculateBondData();
        DetermineBondsAndRings();
      }
      bool any( false), all( true);
      for( auto itr_iso( ISOMORPHISMS.Begin()), itr_iso_end( ISOMORPHISMS.End()); itr_iso != itr_iso_end; ++itr_iso)
      {
        const size_t position( m_Fragment->GetFragmentCSI().Find( *itr_iso));
        bool result( FRAGMENT.ContainsIncompleteRings( position));
        any = any || result;
        all = all && result;
      }
      if( any != all && !m_Fragment->GetNumberChainBonds())
      {
        m_Fragment->GetFragment().WriteMDL( util::GetLogger());
        BCL_Exit( "For a given set of central bonds, either it should have incomplete rings or it should not...not both!", -1);
      }
      m_IncompleteRings = any;
      if
      (
        !m_Fragment->ContainsRingConformations() && !m_Fragment->GetNumberChainBonds()
        && m_Fragment->GetNumberRingBonds() && ISOMORPHISMS.GetSize() > size_t( 1)
        && !m_IncompleteRings
        && m_Fragment->GetFragment().GetSize() > size_t( 3)
        && FRAGMENT.GetRotamerBins().GetSize()
      )
      {
        m_HasDifferentRingRotamers = false;
        for( auto itr( FRAGMENT.GetRotamerBins()( 0).Begin()), itr_end( FRAGMENT.GetRotamerBins()( 0).End()); itr != itr_end; ++itr)
        {
          if( *itr != 12 && *itr != 6)
          {
            // check for non-planar ring region
            m_HasDifferentRingRotamers = true;
            break;
          }
        }
      }
      m_IsomorphismWeights
        = RotamerClusterCenter::GetIsomorphismWeights( m_Fragment->GetFragmentCSI(), MOLECULE.GetSize());
      m_IsomorphismWeights.SetToSum( 1.0);
    }

    //! @brief Clone function
    //! @return pointer to new RotamerDihedralBondData
    RotamerDihedralBondData *RotamerDihedralBondData::Clone() const
    {
      return new RotamerDihedralBondData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &RotamerDihedralBondData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return center bonds of the molecule of interest
    //! @return center bonds of the molecule of interest
    const storage::Vector< graph::UndirectedEdge< ConfigurationalBondType> > &RotamerDihedralBondData::GetMoleculeCenterBonds() const
    {
      return m_DihedralEdges;
    }

    //! @brief return fragment for which this class stores dihedral bond information
    //! @return fragment for which this class stores dihedral bond information
    const SmallMoleculeFragmentIsomorphism &RotamerDihedralBondData::GetFragment() const
    {
      return *m_Fragment;
    }

    //! @brief returns rotamers of the fragment whose dihedral information is stored by this class
    //! @return the rotamers of the fragment whose dihedral information is stored by this class
    const storage::Vector< linal::Vector< double> > &RotamerDihedralBondData::GetRotamers() const
    {
      return m_Fragment->GetRotamers();
    }

    //! @brief returns counts for each rotamer of the fragment
    //! @return counts for each rotamer of the fragment
    const storage::Vector< double> &RotamerDihedralBondData::GetRotamerCounts() const
    {
      return m_Fragment->GetRotamerCounts();
    }

    //! @brief returns all isomorphism of this fragment for the same set of atoms in the parent molecule
    //! @return isomorphism of this fragment for the same set of atoms in the parent molecule
    const storage::Vector< storage::Vector< size_t> > &RotamerDihedralBondData::GetIsomorphisms() const
    {
      return m_Isomorphisms;
    }

    //! @brief returns dihedral bonds of the molecule of interest that this fragment isomorphism represents
    //! @return dihedral bonds of the molecule of interest that this fragment represents
    const storage::Vector< storage::Vector< storage::VectorND< 4, size_t> > > &RotamerDihedralBondData::GetRotamerBonds() const
    {
      return m_DihedralAtomIndices;
    }

    //! @brief returns center bonds that this fragment represents in the molecule of interest
    //! @return center bonds of that this fragment represents in the molecule of interest
    const storage::Vector< graph::UndirectedEdge< ConfigurationalBondType> > &RotamerDihedralBondData::GetCenterBonds() const
    {
      return m_CenterBonds;
    }

    //! @brief returns the CenterBonds that are not contained in rings
    //! @return the CenterBonds that are not contained in rings
    const storage::Vector< size_t> &RotamerDihedralBondData::GetNonRingBonds() const
    {
      return m_NonRingBonds;
    }

    //! @brief returns isomorphism between fragment center bond and molecule bond
    //! @return isomorphism between fragment center bond and molecule bond
    const storage::Vector< size_t> &RotamerDihedralBondData::GetCenterBondIsomorphism() const
    {
      return m_CenterBondIsomorphisms( 0);
    }

    //! @brief returns isomorphism between fragment center bond and molecule bond
    //! @return isomorphism between fragment center bond and molecule bond
    const storage::Vector< storage::Vector< size_t> > &RotamerDihedralBondData::GetCenterBondIsomorphisms() const
    {
      return m_CenterBondIsomorphisms;
    }

    //! @brief returns true if fragment contains rings
    //! @return true if fragment contains ring otherwise false
    bool RotamerDihedralBondData::ContainsRings() const
    {
      return m_ContainsRings;
    }

    //! @brief returns true if fragment can be used to create rings with different conformations of the overall molecule
    //! @return true if fragment can be used to create rings with different conformations of the overall molecule
    bool RotamerDihedralBondData::ContainsRingRotamers() const
    {
      return m_HasDifferentRingRotamers;
    }

    //! @brief returns the bin size strategy used for representing conformation of fragment
    //! @return the bin size strategy used for representing conformation of fragment
    const double &RotamerDihedralBondData::GetBinSize() const
    {
      return m_Fragment->GetBinSize();
    }

    //! @brief returns the total number of times the rotatable element has been seen in the CSD
    //! @return the total number of times the rotatable element is seen in the CSD
    const double &RotamerDihedralBondData::GetFragmentCounts() const
    {
      return m_Fragment->GetFragmentCounts();
    }

    //! @brief returns weights for each isomorphism
    const linal::Vector< double> &RotamerDihedralBondData::GetIsomorphismWeights() const
    {
      return m_IsomorphismWeights;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief determine atoms involved in rotatable bonds and atoms involved in center bonds
    //! @return void
    void RotamerDihedralBondData::CalculateBondData()
    {
      // get the central bond of dihedral angles of the fragment of interest
      const storage::Vector< graph::UndirectedEdge< ConfigurationalBondType> > &fragment_dihedral_edge
      (
        m_Fragment->GetDihedralEdges()
      );

      // all isomorphisms will yield the same central bonds so take the first one
      storage::Vector< storage::Vector< size_t> >::const_iterator itr_isomorphisms( m_Isomorphisms.Begin());

      // get priority dihedral angles of fragment
      const storage::Vector< storage::VectorND< 4, size_t> > &dihedral_angle_info( m_Fragment->GetDihedralAngleIndices());

      bool first_iso( true);
      // store priority dihedral angle information of fragment in terms of fragment isomorphism with the molecule of interest
      for
      (
        storage::Vector< storage::Vector< size_t> >::const_iterator itr_isomorphisms_end( m_Isomorphisms.End());
        itr_isomorphisms != itr_isomorphisms_end;
        ++itr_isomorphisms
      )
      {
        // get isomorphism between fragment and the molecule of interest
        const storage::Vector< size_t> &isomorphism( *itr_isomorphisms);

        // central bonds of fragment in terms of molecule indices
        m_CenterBondIsomorphisms.PushBack();
        for
        (
          storage::Vector< graph::UndirectedEdge< ConfigurationalBondType> >::const_iterator
            itr( fragment_dihedral_edge.Begin()), itr_end( fragment_dihedral_edge.End());
          itr != itr_end;
          ++itr
        )
        {
          graph::UndirectedEdge< ConfigurationalBondType> bond_iso
          (
            isomorphism( itr->GetIndexLow()),
            isomorphism( itr->GetIndexHigh()),
            itr->GetEdgeData()
          );
          if( first_iso)
          {
            m_CenterBonds.PushBack( bond_iso);
          }

          // determine isomorphism between center bonds of fragment and those of molecule, Needed in calculating score in
          // class FragmentCorrelationProbability
          size_t bond_index( 0);
          for
          (
            storage::Vector< graph::UndirectedEdge< ConfigurationalBondType> >::const_iterator
              itr_mol_center( m_DihedralEdges.Begin()), itr_mol_center_end( m_DihedralEdges.End());
            itr_mol_center != itr_mol_center_end;
            ++itr_mol_center, ++bond_index
          )
          {
            if( itr_mol_center->GetIndexLow() == bond_iso.GetIndexLow() && itr_mol_center->GetIndexHigh() == bond_iso.GetIndexHigh())
            {
              m_CenterBondIsomorphisms.LastElement().PushBack( bond_index);
              break;
            }
          }
        }
        first_iso = false;

        // create a vector to store  dihedral bond atom indices of molecule ( this is isomorphic to fragment dihedral)
        storage::Vector< storage::VectorND< 4, size_t> > molecule_dihedral;

        // get atom indices from the frame of reference of the molecule of interest
        for
        (
          storage::Vector< storage::VectorND< 4, size_t> >::const_iterator
            itr( dihedral_angle_info.Begin()), itr_end( dihedral_angle_info.End());
          itr != itr_end;
          ++itr
        )
        {
          const storage::VectorND< 4, size_t> &atom_indices_fragment( *itr);
          storage::VectorND< 4, size_t> atom_indices_molecule
          (
            isomorphism( atom_indices_fragment.First()),
            isomorphism( atom_indices_fragment.Second()),
            isomorphism( atom_indices_fragment.Third()),
            isomorphism( atom_indices_fragment.Fourth())
          );
          molecule_dihedral.PushBack( atom_indices_molecule);
        }
        m_DihedralAtomIndices.PushBack( molecule_dihedral);
      }
      return;
    }

    //! @brief determine rings and bonds contained in the fragment
    //! @return void
    void RotamerDihedralBondData::DetermineBondsAndRings()
    {
      // get list of vertices of fragment that are in a ring and contained in the molecule of interest
      storage::List< storage::Vector< size_t> > ring_vertices
      (
        SmallMoleculeFragmentMapping::GetRingVertices( m_Fragment->GetFragment())
      );

      // get the isomorphism between molecule and fragment
      const storage::Vector< size_t> &isomorphism( *m_Isomorphisms.Begin());

      storage::List< storage::Vector< size_t> > corresponding_ring_mol;
      storage::Set< size_t> ring_bonds;

      // iterate through the rings contained in the fragment of interest
      for
      (
        storage::List< storage::Vector< size_t> >::const_iterator
          itr_ring( ring_vertices.Begin()), itr_ring_end( ring_vertices.End());
        itr_ring != itr_ring_end;
        ++itr_ring
      )
      {
        // get vertices contained in the current ring
        const storage::Vector< size_t> &cur_ring( *itr_ring);

        // corresponding ring in the molecule whose conformations needs to be sampled
        storage::Vector< size_t> cur_corresponding_ring;

        // a container to store bonds that belong to current ring
        storage::List< size_t> cur_ring_bonds;

        // get ring in molecule of interest that corresponds to current ring being iterated
        for
        (
          storage::Vector< size_t>::const_iterator itr_cur( cur_ring.Begin()), itr_cur_end( cur_ring.End());
          itr_cur != itr_cur_end;
          ++itr_cur
        )
        {
          cur_corresponding_ring.PushBack( isomorphism( *itr_cur));
        }

        // go through the center bonds of the molecule that the fragment contains
        size_t bond_index( 0);
        for
        (
          storage::Vector< graph::UndirectedEdge< ConfigurationalBondType> >::const_iterator
            itr_center_bonds( m_CenterBonds.Begin()), itr_center_bonds_end( m_CenterBonds.End());
          itr_center_bonds != itr_center_bonds_end;
          ++itr_center_bonds, ++bond_index
        )
        {
          // if ring vertices contains center bond atoms then the bond is contained in the current ring
          if
          (
            IfFullIntersection
            (
              cur_corresponding_ring,
              storage::Vector< size_t>::Create( itr_center_bonds->GetIndexLow(), itr_center_bonds->GetIndexHigh())
            )
          )
          {
            cur_ring_bonds.PushBack( bond_index);
            ring_bonds.Insert( bond_index);
          }
        }
      }

      size_t number_of_bonds( m_CenterBonds.GetSize());

      // now determine which bonds are not in rings
      if( !ring_bonds.IsEmpty())
      {
        for( size_t count( 0); count < number_of_bonds; ++count)
        {
          if( !ring_bonds.Contains( count))
          {
            m_NonRingBonds.PushBack( count);
          }
        }
      }
      else
      {
        for( size_t count( 0); count < number_of_bonds; ++count)
        {
          m_NonRingBonds.PushBack( count);
        }
      }
    }

    //! @brief check whether elements are fully contained in the container of interest
    //! @param CONTAINER container which needs to be checked to see if it contains elements
    //! @param ELEMENTS elements that need to be checked to see if they are fully contained in the CONTAINER
    //! @return true if elements completely contained in the container otherwise false.
    bool RotamerDihedralBondData::IfFullIntersection
    (
      const storage::Vector< size_t> &CONTAINER,
      const linal ::Vector< size_t> &ELEMENTS
    )
    {
      std::vector< size_t> intersec_vector( ELEMENTS.GetSize());
      std::vector< size_t>::iterator itr_instersec;
      storage::Set< size_t> container_set( CONTAINER.Begin(), CONTAINER.End());
      storage::Set< size_t> element_set( ELEMENTS.Begin(), ELEMENTS.End());

      itr_instersec = std::set_intersection
                      (
                        container_set.Begin(),
                        container_set.End(),
                        element_set.Begin(),
                        element_set.End(),
                        intersec_vector.begin()
                      );
      size_t intersection_size( itr_instersec - intersec_vector.begin());
      if( intersection_size == ELEMENTS.GetSize())
      {
        return true;
      }
      else
      {
        return false;
      }
    }

    namespace
    {
      float GetSimilarity( const ConfigurationalBondType &BOND_TYPE, const size_t &ROTAMER_BIN, const size_t MODEL_BIN)
      {
        if( BOND_TYPE->IsBondInRing() || BOND_TYPE->GetNumberOfElectrons() == size_t( 6))
        {
          return ROTAMER_BIN == MODEL_BIN ? 1.0 : 0.0;
        }
        else if( ROTAMER_BIN == MODEL_BIN)
        {
          return 1.0;
        }
        static const float s_NonConjSimilarity[ 12][ 12] =
        {
          { 1.00, 0.80, 0.40, 1.00, 0.81, 1.00, 0.12, 0.50, 0.73, 1.00, 1.00, 0.78},
          { 0.02, 1.00, 0.01, 0.03, 0.02, 0.34, 0.01, 0.02, 0.02, 0.59, 0.02, 0.02},
          { 0.95, 1.00, 1.00, 1.00, 1.00, 0.92, 0.40, 0.88, 1.00, 1.00, 1.00, 1.00},
          { 0.64, 0.65, 0.38, 1.00, 0.67, 0.57, 0.12, 0.47, 0.66, 0.64, 0.64, 0.62},
          { 0.84, 0.81, 0.46, 1.00, 1.00, 0.75, 0.18, 0.57, 0.84, 0.81, 0.87, 0.88},
          { 0.08, 0.76, 0.03, 0.06, 0.05, 1.00, 0.01, 0.03, 0.05, 0.29, 0.19, 0.04},
          { 0.51, 0.76, 0.75, 0.83, 0.76, 0.82, 1.00, 1.00, 0.95, 0.84, 0.78, 0.77},
          { 0.71, 0.84, 0.54, 1.00, 0.78, 0.72, 0.33, 1.00, 0.92, 0.73, 0.82, 0.89},
          { 0.76, 0.82, 0.54, 1.00, 0.84, 0.77, 0.22, 0.67, 1.00, 0.85, 0.84, 0.83},
          { 0.06, 0.92, 0.02, 0.05, 0.03, 0.20, 0.01, 0.02, 0.04, 1.00, 0.13, 0.03},
          { 0.46, 0.31, 0.17, 0.40, 0.32, 1.00, 0.07, 0.22, 0.31, 1.00, 1.00, 0.33},
          { 0.80, 0.79, 0.44, 1.00, 0.87, 0.70, 0.18, 0.64, 0.82, 0.76, 0.87, 1.00}
        };
        static const float s_ConjSingleSimilarity[ 12][ 12] =
        {
          { 1.00, 0.25, 0.11, 0.26, 0.82, 0.66, 0.58, 0.09, 0.13, 0.19, 0.89, 0.75},
          { 0.73, 1.00, 0.37, 0.79, 0.76, 0.74, 0.68, 0.33, 0.66, 0.67, 0.50, 0.64},
          { 0.89, 0.98, 1.00, 0.89, 0.99, 0.74, 0.89, 0.28, 0.95, 1.00, 0.77, 1.00},
          { 0.62, 0.64, 0.27, 1.00, 0.68, 0.85, 0.60, 0.29, 0.37, 0.45, 0.49, 0.49},
          { 0.35, 0.11, 0.05, 0.12, 1.00, 0.39, 0.36, 0.05, 0.06, 0.08, 0.92, 0.11},
          { 0.21, 0.08, 0.03, 0.11, 0.28, 1.00, 0.17, 0.04, 0.03, 0.05, 0.29, 0.12},
          { 1.00, 0.40, 0.20, 0.44, 1.00, 0.95, 1.00, 0.17, 0.23, 0.28, 1.00, 0.83},
          { 1.00, 1.00, 0.44, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.97, 0.98, 1.00},
          { 0.52, 0.90, 0.49, 0.63, 0.60, 0.42, 0.53, 0.35, 1.00, 0.76, 0.39, 0.76},
          { 0.72, 0.87, 0.50, 0.72, 0.72, 0.66, 0.61, 0.31, 0.71, 1.00, 0.55, 0.79},
          { 0.41, 0.08, 0.04, 0.09, 0.98, 0.42, 0.36, 0.04, 0.04, 0.07, 1.00, 0.17},
          { 0.11, 0.03, 0.02, 0.03, 0.04, 0.06, 0.07, 0.01, 0.03, 0.03, 0.05, 1.00}
        };
        static const float s_ConjDoubleSimilarity[ 12][ 12] =
        {
          { 1.00, 0.24, 0.19, 0.18, 0.13, 0.79, 0.14, 0.09, 0.20, 0.16, 0.68, 0.10},
          { 0.03, 1.00, 0.30, 0.03, 0.34, 0.38, 0.01, 0.01, 0.04, 0.55, 0.03, 0.00},
          { 0.07, 0.77, 1.00, 0.05, 0.03, 0.71, 0.02, 0.02, 0.20, 0.14, 0.08, 0.00},
          { 0.16, 0.16, 0.11, 1.00, 0.27, 0.16, 0.17, 0.31, 0.16, 0.16, 0.12, 0.40},
          { 0.04, 0.68, 0.02, 0.09, 1.00, 0.34, 0.10, 0.02, 0.03, 0.80, 0.06, 0.15},
          { 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.01},
          { 0.11, 0.06, 0.05, 0.15, 0.29, 0.50, 1.00, 0.06, 0.06, 0.55, 0.05, 0.14},
          { 0.24, 0.17, 0.12, 0.90, 0.17, 0.16, 0.18, 1.00, 0.17, 0.17, 0.13, 0.46},
          { 0.25, 0.35, 0.67, 0.22, 0.11, 0.19, 0.09, 0.08, 1.00, 0.35, 0.29, 0.02},
          { 0.03, 0.78, 0.08, 0.04, 0.56, 0.34, 0.14, 0.01, 0.06, 1.00, 0.03, 0.00},
          { 0.41, 0.14, 0.13, 0.08, 0.13, 1.00, 0.04, 0.03, 0.14, 0.08, 1.00, 0.16},
          { 0.00, 0.00, 0.00, 0.00, 0.00, 0.01, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00}
        };
        if( BOND_TYPE->GetNumberOfElectrons() == size_t( 4) || BOND_TYPE->GetConjugation() == ConstitutionalBondTypeData::e_Amide)
        {
          return s_ConjDoubleSimilarity[ ROTAMER_BIN - 1][ MODEL_BIN - 1];
        }
        else if( BOND_TYPE->GetConjugation() == ConstitutionalBondTypeData::e_Conjugated)
        {
          return s_ConjSingleSimilarity[ ROTAMER_BIN - 1][ MODEL_BIN - 1];
        }
        return s_NonConjSimilarity[ ROTAMER_BIN - 1][ MODEL_BIN - 1];
      }
    }

    //! @brief determine the imputed counts of a given rotamer
    //! @return imputed counts of a given rotamer
    double RotamerDihedralBondData::GetImputedRotamerCounts( const linal::Vector< int> &BINS) const
    {
      double max_count( 0.0);
      double total_count( 0.0);
      const size_t fnd( m_Fragment->GetRotamerBins().Find( BINS));
      if( fnd < m_Fragment->GetRotamerNumbers())
      {
        return m_Fragment->GetRotamerCounts()( fnd);
      }
      for( size_t rotamer_number( 0), n_rotamers( m_Fragment->GetRotamerNumbers()); rotamer_number < n_rotamers; ++rotamer_number)
      {
        auto itr_bonds( m_CenterBonds.Begin());
        double est_counts( m_Fragment->GetRotamerCounts()( rotamer_number));
        auto itr_bins( BINS.Begin());
        const linal::Vector< size_t> &rot_bins( m_Fragment->GetRotamerBins()( rotamer_number));
        for( auto itr_b( rot_bins.Begin()), itr_b_end( rot_bins.End()); itr_b != itr_b_end; ++itr_b, ++itr_bonds, ++itr_bins)
        {
          est_counts *= GetSimilarity( itr_bonds->GetEdgeData(), *itr_b, *itr_bins);
        }
        max_count = std::max( est_counts, max_count);
        total_count += est_counts;
      }
      return max_count;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RotamerDihedralBondData::Read( std::istream &ISTREAM)
    {

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &RotamerDihedralBondData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      if( m_Fragment.IsDefined())
      {
        io::Serialize::Write( m_DihedralEdges, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( *m_Fragment, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Isomorphisms, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_DihedralAtomIndices, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_CenterBonds, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_NonRingBonds, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_CenterBondIsomorphisms, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_ContainsRings, OSTREAM, INDENT) << '\t';
        io::Serialize::Write( m_HasDifferentRingRotamers, OSTREAM, INDENT) << '\t';
        io::Serialize::Write( m_IncompleteRings, OSTREAM, INDENT) << '\t';
        io::Serialize::Write( m_Uniqueness, OSTREAM, INDENT) << '\n';
      }

      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
