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
#include "chemistry/bcl_chemistry_priority_dihedral_angles.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_substituent_conformational.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_angle.h"
#include "storage/bcl_storage_triplet.h"

//#define BCL_PROFILE_PriorityDihedralAngles
#ifdef BCL_PROFILE_PriorityDihedralAngles
#include "util/bcl_util_stopwatch.h"
#endif

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    //! @brief get the wrap around angle; used for determining when to wrap angles in CalculateMinimumDihedralAngle
    //! @return the wrapping angle
    double &PriorityDihedralAngles::GetChangeableWrappingAngle()
    {
      static double s_wrapping_angle( 15.0);
      return s_wrapping_angle;
    }

    //! @brief Clone function
    //! @return pointer to new PriorityDihedralAngles
    PriorityDihedralAngles *PriorityDihedralAngles::Clone() const
    {
      return new PriorityDihedralAngles( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &PriorityDihedralAngles::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the wrap around angle; used for determining when to wrap angles in CalculateMinimumDihedralAngle
    //! @return the wrapping angle
    const double &PriorityDihedralAngles::GetWrappingAngle()
    {
      return GetChangeableWrappingAngle();
    }

    //! @brief set the wrap around angle; used for determining when to wrap angles in CalculateMinimumDihedralAngle
    //! @param WRAPPING_ANGLE the wrapping angle
    void PriorityDihedralAngles::SetWrappingAngle( const double &WRAPPING_ANGLE)
    {
      static sched::Mutex s_Mutex;
      s_Mutex.Lock();
      GetChangeableWrappingAngle() = WRAPPING_ANGLE;
      s_Mutex.Unlock();
    }

    //! @brief calculate prioritized dihedral angle and indices of atoms that make the angle
    //! @param MOLECULE the molecule of interest
    //! @return a pair of vector containing wrapped around dihedral angles in degrees, in ascending order of atom positions
    storage::Pair< storage::Vector< double>, storage::Vector< storage::VectorND< 4, size_t> > >
    PriorityDihedralAngles::operator()( const ConformationInterface &MOLECULE) const
    {
      util::SiPtr< const linal::Vector< float> > cached_priorities( MOLECULE.FindInCache( util::ObjectDataLabel( "PriorityAtoms")));
      if( !cached_priorities.IsDefined())
      {
        RecalculatePriority( MOLECULE);
        cached_priorities = MOLECULE.FindInCache( util::ObjectDataLabel( "PriorityAtoms"));
      }

      // begin a vector of dihedral angles
      storage::Pair< storage::Vector< double>, storage::Vector< storage::VectorND< 4, size_t> > > dihedral_angles;

      // compute the angle A-B-C-D, where
      // B is the position of the first atom in the bond
      // C is the position of the second atom in the bond
      // A is the positions of the set of atoms that the first atom in the bond is connected to, except C
      // D is the positions of the set of atoms that the second atom in the bond is connected to, except B

      // iterate through atoms to find the an atom involved in the middle bond of dihedral bond
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atom_b( MOLECULE.GetAtomsIterator());
        itr_atom_b.NotAtEnd();
        ++itr_atom_b
      )
      {
        const AtomConformationalInterface &atom_b( *itr_atom_b);
        const storage::Vector< BondConformational> &connected_atoms_b( atom_b.GetBonds());
        if( connected_atoms_b.GetSize() < 2)
        {
          continue;
        }
        // get atoms connected to atom b
        for
        (
          storage::Vector< BondConformational>::const_iterator
            itr_bond_b( connected_atoms_b.Begin()), itr_bond_b_end( connected_atoms_b.End());
          itr_bond_b != itr_bond_b_end;
          ++itr_bond_b
        )
        {
          const bool bond_in_ring( itr_bond_b->GetBondType()->IsBondInRing());
          const AtomConformationalInterface &atom_c( itr_bond_b->GetTargetAtom());
          // ensure that atom c is before atom b; this orders the dihedral angles so that we don't see each
          // dihedral angle twice

          if( &atom_c > &atom_b)
          {
            continue;
          }

          const storage::Vector< BondConformational> &connected_atoms_c( atom_c.GetBonds());
          if( connected_atoms_c.GetSize() == size_t( 1))
          {
            continue;
          }
          // get priority atoms connected to atom b
          storage::Vector< size_t> atom_priorities_a
          (
            DetermineSubstituentPriority( MOLECULE, *cached_priorities, connected_atoms_b, atom_c, bond_in_ring)
          );

          // get priority atoms connected to atom c
          storage::Vector< size_t> atom_priorities_d
          (
            DetermineSubstituentPriority( MOLECULE, *cached_priorities, connected_atoms_c, atom_b, bond_in_ring)
          );

          // calculate priority dihedral angle from the priority atom and lowest angle measure
          storage::Pair< double, storage::VectorND< 4, size_t> > dihedral_angle_info
          (
            CalculateMinimumDihedralAngle( MOLECULE, atom_priorities_a, atom_b, atom_c, atom_priorities_d, bond_in_ring)
          );
          dihedral_angles.First().PushBack( dihedral_angle_info.First());
          dihedral_angles.Second().PushBack( dihedral_angle_info.Second());
        }
      }
      return dihedral_angles;
    }

    //! @brief get the multiplicity of the priorities, to identify symmetric positions
    //! @param MOLECULE the molecule of interest
    storage::Vector< size_t> PriorityDihedralAngles::CalculateMultiplicity( const ConformationInterface &MOLECULE)
    {
      util::SiPtr< const linal::Vector< float> > cached_priorities; //( MOLECULE.FindInCache( util::ObjectDataLabel( "PriorityAtoms")));
      if( !cached_priorities.IsDefined())
      {
        RecalculatePriority( MOLECULE);
        cached_priorities = MOLECULE.FindInCache( util::ObjectDataLabel( "PriorityAtoms"));
      }
      storage::Vector< size_t> priority_counts( MOLECULE.GetSize(), size_t( 0));
      for( size_t priorities( 0), sz( MOLECULE.GetSize()); priorities < sz; ++priorities)
      {
        ++priority_counts( size_t( cached_priorities->operator ()( priorities)));
      }
      storage::Vector< size_t> priority_multiplicities( priority_counts.GetSize(), size_t( 0));
      for( size_t priorities( 0), sz( MOLECULE.GetSize()); priorities < sz; ++priorities)
      {
        priority_multiplicities( priorities) = priority_counts( size_t( cached_priorities->operator ()( priorities)));
      }
      return priority_multiplicities;
    }

    //! @brief get the priorities, to identify symmetric positions
    //! @param MOLECULE the molecule of interest
    storage::Vector< size_t> PriorityDihedralAngles::GetPriority( const ConformationInterface &MOLECULE)
    {
      util::SiPtr< const linal::Vector< float> > cached_priorities
      (
        MOLECULE.FindInCache( util::ObjectDataLabel( "PriorityAtoms"))
      );
      if( !cached_priorities.IsDefined())
      {
        RecalculatePriority( MOLECULE);
        cached_priorities = MOLECULE.FindInCache( util::ObjectDataLabel( "PriorityAtoms"));
      }
      return storage::Vector< size_t>( cached_priorities->Begin(), cached_priorities->End());
    }

    //! @brief calculate maximal dihedral angle differences
    //! @param MOLECULE_A, MOLECULE_B the molecules to compare
    //! @return maximal dihedral angle differences
    storage::Vector< double> PriorityDihedralAngles::operator()
    (
      const ConformationInterface &MOLECULE_A,
      const ConformationInterface &MOLECULE_B
    ) const
    {
      util::SiPtr< const linal::Vector< float> > cached_priorities_a( MOLECULE_A.FindInCache( util::ObjectDataLabel( "PriorityAtoms")));
      if( !cached_priorities_a.IsDefined())
      {
        RecalculatePriority( MOLECULE_A);
        cached_priorities_a = MOLECULE_A.FindInCache( util::ObjectDataLabel( "PriorityAtoms"));
      }

      util::SiPtr< const linal::Vector< float> > cached_priorities_b( MOLECULE_B.FindInCache( util::ObjectDataLabel( "PriorityAtoms")));
      if( !cached_priorities_b.IsDefined())
      {
        RecalculatePriority( MOLECULE_B);
        cached_priorities_b = MOLECULE_B.FindInCache( util::ObjectDataLabel( "PriorityAtoms"));
      }

      #ifdef BCL_PROFILE_PriorityDihedralAngles
      static util::Stopwatch s_angles( "PDA: Comparing angles", util::Time( 1, 0), util::Message::e_Standard, true, false);
      s_angles.Start();
      #endif

      // begin a vector of dihedral angles
      storage::Vector< double> dihedral_angles;

      // compute the angle A-B-C-D, where
      // B is the position of the first atom in the bond
      // C is the position of the second atom in the bond
      // A is the positions of the set of atoms that the first atom in the bond is connected to, except C
      // D is the positions of the set of atoms that the second atom in the bond is connected to, except B

      // iterate through atoms to find the an atom involved in the middle bond of dihedral bond
      for
      (
        iterate::Generic< const AtomConformationalInterface>
          itr_atom_b_a( MOLECULE_A.GetAtomsIterator()), itr_atom_b_b( MOLECULE_B.GetAtomsIterator());
        itr_atom_b_a.NotAtEnd();
        ++itr_atom_b_a, ++itr_atom_b_b
      )
      {
        const AtomConformationalInterface &atom_b_a( *itr_atom_b_a);
        const storage::Vector< BondConformational> &connected_atoms_b_a( atom_b_a.GetBonds());
        if( connected_atoms_b_a.GetSize() < 2)
        {
          continue;
        }
        const AtomConformationalInterface &atom_b_b( *itr_atom_b_b);
        const storage::Vector< BondConformational> &connected_atoms_b_b( atom_b_b.GetBonds());
        // get atoms connected to atom b
        for
        (
          storage::Vector< BondConformational>::const_iterator
            itr_bond_b_a( connected_atoms_b_a.Begin()), itr_bond_b_end( connected_atoms_b_a.End()),
            itr_bond_b_b( connected_atoms_b_b.Begin());
          itr_bond_b_a != itr_bond_b_end;
          ++itr_bond_b_a, ++itr_bond_b_b
        )
        {
          const bool bond_in_ring( itr_bond_b_a->GetBondType()->IsBondInRing());
          const AtomConformationalInterface &atom_c_a( itr_bond_b_a->GetTargetAtom());
          // ensure that atom c is before atom b; this orders the dihedral angles so that we don't see each
          // dihedral angle twice

          if( &atom_c_a > &atom_b_a)
          {
            continue;
          }

          const storage::Vector< BondConformational> &connected_atoms_c_a( atom_c_a.GetBonds());
          if( connected_atoms_c_a.GetSize() == size_t( 1))
          {
            continue;
          }
          const AtomConformationalInterface &atom_c_b( itr_bond_b_b->GetTargetAtom());
          const storage::Vector< BondConformational> &connected_atoms_c_b( atom_c_b.GetBonds());
          // get priority atoms connected to atom b
          storage::Vector< size_t> atom_priorities_a_a
          (
            DetermineSubstituentPriority( MOLECULE_A, *cached_priorities_a, connected_atoms_b_a, atom_c_a, bond_in_ring)
          );
          storage::Vector< size_t> atom_priorities_a_b
          (
            DetermineSubstituentPriority( MOLECULE_B, *cached_priorities_b, connected_atoms_b_b, atom_c_b, bond_in_ring)
          );

          // get priority atoms connected to atom c
          storage::Vector< size_t> atom_priorities_d_a
          (
            DetermineSubstituentPriority( MOLECULE_A, *cached_priorities_a, connected_atoms_c_a, atom_b_a, bond_in_ring)
          );
          storage::Vector< size_t> atom_priorities_d_b
          (
            DetermineSubstituentPriority( MOLECULE_B, *cached_priorities_b, connected_atoms_c_b, atom_b_b, bond_in_ring)
          );

          storage::Vector< double> dihedrals_a
          (
            CalculateOrderedDihedralAngles
            (
              MOLECULE_A,
              atom_priorities_a_a,
              atom_b_a,
              atom_c_a,
              atom_priorities_d_a,
              bond_in_ring
            )
          );
          storage::Vector< double> dihedrals_b
          (
            CalculateOrderedDihedralAngles
            (
              MOLECULE_B,
              atom_priorities_a_b,
              atom_b_b,
              atom_c_b,
              atom_priorities_d_b,
              bond_in_ring
            )
          );
          dihedral_angles.PushBack( CalculateMaximumDihedralAngleDifference( dihedrals_a, dihedrals_b));
        }
      }
      #ifdef BCL_PROFILE_PriorityDihedralAngles
      s_angles.Stop();
      #endif
      return dihedral_angles;
    }

    //! @brief calculate priority of all atoms in the molecule and populate the member variable that stores priority
    //! @param MOLECULE the molecule of interest
    void PriorityDihedralAngles::RecalculatePriority( const ConformationInterface &MOLECULE)
    {
      #ifdef BCL_PROFILE_PriorityDihedralAngles
      static util::Stopwatch s_initial( "Recalculating priority", util::Time( 1, 0), util::Message::e_Standard, true, false);
      s_initial.Start();
      #endif
      storage::Vector< AtomType> atom_types( MOLECULE.GetAtomTypesVector());
      storage::Vector< sdf::BondInfo> bond_info( MOLECULE.GetBondInfo());

      static sched::Mutex s_mutex;
      s_mutex.Lock();
      static size_t s_priority_queue_size( 0);
      typedef storage::Triplet
      <
        storage::Vector< AtomType>,
        storage::Vector< sdf::BondInfo>,
        descriptor::CacheMap::value_type
      > t_PriorityTriplet;
      static storage::List< t_PriorityTriplet> s_priority_queue;

      for
      (
        storage::List< t_PriorityTriplet>::iterator
          itr( s_priority_queue.Begin()), itr_end( s_priority_queue.End());
        itr != itr_end;
        ++itr
      )
      {
        if( itr->First() == atom_types && itr->Second() == bond_info)
        {
          MOLECULE.Cache( util::ObjectDataLabel( "PriorityAtoms"), itr->Third());
          if( itr != s_priority_queue.Begin())
          {
            s_priority_queue.InternalData().splice( s_priority_queue.Begin(), s_priority_queue.InternalData(), itr);
          }
          s_mutex.Unlock();
          #ifdef BCL_PROFILE_PriorityDihedralAngles
          s_initial.Stop();
          #endif
          return;
        }
      }
      s_mutex.Unlock();

      // create and allocate memory for vector to store  priority atoms
      storage::Vector< SubstituentConformational> holders;
      holders.AllocateMemory( MOLECULE.GetNumberAtoms());

      // push back atoms into the vector
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atoms( MOLECULE.GetAtomsIterator());
        itr_atoms.NotAtEnd();
        ++itr_atoms
      )
      {
        holders.PushBack( SubstituentConformational( *itr_atoms));
      }

      // sort the vector of substituents
      holders.Sort( std::less< SubstituentConformational>());

      linal::Vector< float> atom_priority( MOLECULE.GetNumberAtoms(), float( 0));

      MOLECULE.Uncache( util::ObjectDataLabel( "PriorityAtoms"));
      if( holders.IsEmpty())
      {
        #ifdef BCL_PROFILE_PriorityDihedralAngles
        s_initial.Stop();
        #endif
        MOLECULE.Cache( util::ObjectDataLabel( "PriorityAtoms"), atom_priority);
        return;
      }

      storage::Vector< SubstituentConformational>::const_iterator last_unique_substituent_itr( holders.Begin());

      size_t current_priority( 0);

      // make an iterator to walk ahead of last_unique_substituent_itr
      storage::Vector< SubstituentConformational>::const_iterator itr( last_unique_substituent_itr);
      ++itr;

      // store atom index and its priority in m_AtomPriority
      for
      (
        storage::Vector< SubstituentConformational>::const_iterator
          itr_end( holders.End());
        itr != itr_end;
        ++itr
      )
      {
        if( *last_unique_substituent_itr < *itr) // is *itr a unique substituent?
        {
          ++current_priority;
          last_unique_substituent_itr = itr;
        }
        atom_priority( MOLECULE.GetAtomIndex( *itr->GetRootAtom())) = current_priority;
      }

      MOLECULE.Cache( util::ObjectDataLabel( "PriorityAtoms"), atom_priority);
      s_mutex.Lock();
      if( s_priority_queue_size >= size_t( 32))
      {
        s_priority_queue.PopBack();
      }
      s_priority_queue.PushFront( t_PriorityTriplet ( atom_types, bond_info, atom_priority));
      s_mutex.Unlock();
      #ifdef BCL_PROFILE_PriorityDihedralAngles
      s_initial.Stop();
      #endif
    }

    //! @brief calculate the center atoms of dihedral edges
    //! @param MOLECULE the molecule of interest
    //! @return a vector of edges that represent center bonds of dihedral bonds
    storage::Vector< graph::UndirectedEdge< ConfigurationalBondType> > PriorityDihedralAngles::GetDihedralEdges
    (
      const ConformationInterface &MOLECULE
    )
    {
      storage::Vector< graph::UndirectedEdge< ConfigurationalBondType> > dihedral_edges;
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atom_b( MOLECULE.GetAtomsIterator());
        itr_atom_b.NotAtEnd();
        ++itr_atom_b
      )
      {
        const AtomConformationalInterface &atom_b( *itr_atom_b);
        const storage::Vector< BondConformational> &connected_atoms_b( atom_b.GetBonds());
        if( connected_atoms_b.GetSize() < 2)
        {
          continue;
        }
        // get atoms connected to atom b
        for
        (
          storage::Vector< BondConformational>::const_iterator
            itr_bond_b( connected_atoms_b.Begin()), itr_bond_b_end( connected_atoms_b.End());
          itr_bond_b != itr_bond_b_end;
          ++itr_bond_b
        )
        {
          const AtomConformationalInterface &atom_c( itr_bond_b->GetTargetAtom());
          // ensure that atom c is before atom b; this orders the dihedral angles so that we don't see each
          // dihedral angle twice

          if( &atom_c > &atom_b)
          {
            continue;
          }

          const storage::Vector< BondConformational> &connected_atoms_c( atom_c.GetBonds());
          if( connected_atoms_c.GetSize() == size_t( 1))
          {
            continue;
          }
          dihedral_edges.PushBack
          (
            graph::UndirectedEdge< ConfigurationalBondType>
            (
              MOLECULE.GetAtomIndex( atom_b),
              MOLECULE.GetAtomIndex( atom_c),
              itr_bond_b->GetBondType()
            )
          );
        }
      }
      return dihedral_edges;
    }
    //! @brief return priority atoms among all atoms connected to a single atom
    //! @param MOLECULE the molecule of interest
    //! @param BONDS all the bonds of the atom for which substituent priority has to be determined
    //! @param DIHEDRAL_BONDED_ATOM atom which needs to be excluded from list of priorities
    //! @return priority atoms among all atoms connected to a single atom
    const storage::Vector< size_t> PriorityDihedralAngles::DetermineSubstituentPriority
    (
      const ConformationInterface &MOLECULE,
      const linal::Vector< float> &ATOM_PRIORITY,
      const storage::Vector< BondConformational> &BONDS,
      const AtomConformationalInterface &DIHEDRAL_BONDED_ATOM,
      const bool &IS_RING_BOND
    )
    {
      // create a vector for storing atom index and its priority
      storage::Vector< storage::VectorND< 2, size_t> > atom_index_priorities;
      // go over all the bonds of the atom for whose substituent priority has to be determined
      for
      (
        storage::Vector< BondConformational>::const_iterator itr_atom( BONDS.Begin()), itr_atom_end( BONDS.End());
         itr_atom != itr_atom_end;
        ++itr_atom
      )
      {
        // get index of target atom of bond which are substituents of the atom of interest
        const size_t atom_index( MOLECULE.GetAtomIndex( itr_atom->GetTargetAtom()));
        // check whether the target atom needs to be skipped
        if( &itr_atom->GetTargetAtom() == &DIHEDRAL_BONDED_ATOM)
        {
          continue;
        }

        // if atom is in ring, then use the atom so that we dont choose substituents of the ring
        if( BONDS.GetSize() > 2 && IS_RING_BOND && ( !itr_atom->GetBondType()->IsBondInRing()))
        {
          continue;
        }

        // get global priority of target atom
        const size_t atom_priority( ATOM_PRIORITY( atom_index));
        // store the index and priority and insert it in appropriate place in the atom_index_priorities
        storage::VectorND< 2, size_t> new_element( atom_index, atom_priority);
        // if atom_index_priorities is empty then we can insert directly
        if( atom_index_priorities.IsEmpty())
        {
          atom_index_priorities.PushBack( new_element);
        }
        // if atom_index_priorities not empty then compare priority with existing elements and insert at the right place
        else if( atom_priority < atom_index_priorities.FirstElement().Second())
        {
          atom_index_priorities.InsertElement( atom_index_priorities.Begin(), new_element);
        }
        else if( atom_priority >= atom_index_priorities.LastElement().Second())
        {
          atom_index_priorities.PushBack( new_element);
        }
        else
        {
          atom_index_priorities.InsertElement( atom_index_priorities.Begin() + 1, new_element);
        }
      }

      // return the indices of highest priority atoms, multiple if there are multiple atoms with high prioirity
      if( atom_index_priorities.GetSize() == size_t( 1))
      {
        return storage::Vector< size_t>( 1, atom_index_priorities.FirstElement().First());
      }
      else if( atom_index_priorities.GetSize() == size_t( 2))
      {
        if( atom_index_priorities.FirstElement().Second() != atom_index_priorities.LastElement().Second())
        {
          return storage::Vector< size_t>( 1, atom_index_priorities.FirstElement().First());
        }
        else
        {
          return storage::Vector< size_t>::Create
                (
                  atom_index_priorities.FirstElement().First(),
                  atom_index_priorities.LastElement().First()
                );
        }
      }
      else
      {
        if( atom_index_priorities.FirstElement().Second() != atom_index_priorities.LastElement().Second())
        {
          if( atom_index_priorities( 0).Second() == atom_index_priorities( 1).Second())
          {
            return storage::Vector< size_t>( 1, atom_index_priorities.LastElement().First());
          }
          else
          {
            return storage::Vector< size_t>( 1, atom_index_priorities.FirstElement().First());
          }
        }
        else
        {
          return storage::Vector< size_t>::Create
                (
                  atom_index_priorities( 0).First(),
                  atom_index_priorities( 1).First(),
                  atom_index_priorities( 2).First()
                );
        }
      }
      return storage::Vector< size_t>();
    }

    //! @brief calculate prioritized dihedral angle
    //! @param MOLECULE the molecule of interest
    //! @param ATOM_A_PRIORITIES vector of indices of priority atoms connected to atom B
    //! @param ATOM_B atom involved in the central bond of the dihedral bond
    //! @param ATOM_C atom involved in the central bond of the dihedral bond
    //! @param ATOM_D_PRIORITIES vector of indices of priority atoms connected to atom C
    //! @param IS_RING_BOND true if B-C is in a ring
    storage::Pair< double, storage::VectorND< 4, size_t> > PriorityDihedralAngles::CalculateMinimumDihedralAngle
    (
      const ConformationInterface &MOLECULE,
      const storage::Vector< size_t> &ATOM_A_PRIORITIES,
      const AtomConformationalInterface &ATOM_B,
      const AtomConformationalInterface &ATOM_C,
      const storage::Vector< size_t> &ATOM_D_PRIORITIES,
      const bool &IS_RING_BOND
    ) const
    {
      // get iterators to the two central atoms involved in the dihedral bond
      iterate::Generic< const AtomConformationalInterface> itr_atoms_a( MOLECULE.GetAtomsIterator());
      iterate::Generic< const AtomConformationalInterface> itr_atoms_d( MOLECULE.GetAtomsIterator());

      // create a variable to store the lowest angle for all possible dihedral angle values
      double lowest_angle_measure( std::numeric_limits< double>::infinity());
      size_t atom_index_a( 0);
      size_t atom_index_d( 0);

      // get measurement value for all possible dihedral angles about a central bond
      for
      (
        storage::Vector< size_t>::const_iterator itr_a( ATOM_A_PRIORITIES.Begin()), itr_a_end( ATOM_A_PRIORITIES.End());
        itr_a != itr_a_end;
        ++itr_a
      )
      {
        itr_atoms_a.GotoPosition( *itr_a);
        for
        (
          storage::Vector< size_t>::const_iterator itr_d( ATOM_D_PRIORITIES.Begin()), itr_d_end( ATOM_D_PRIORITIES.End());
          itr_d != itr_d_end;
          ++itr_d
        )
        {
          itr_atoms_d.GotoPosition( *itr_d);

          double angle_measure
          (
            math::Angle::Degree
            (
              linal::Dihedral
              (
                itr_atoms_a->GetPosition(),
                ATOM_B.GetPosition()      ,
                ATOM_C.GetPosition()      ,
                itr_atoms_d->GetPosition()
              )
            )
          );

          if( angle_measure < GetWrappingAngle())
          {
            angle_measure += 360;
          }
          // if the new angle is less than previously measured lowest angle, then store it and the atom indices of
          // atom_a and atom_d which give rise to the lowest value
          if( angle_measure < lowest_angle_measure)
          {
            lowest_angle_measure = angle_measure;
            atom_index_a = MOLECULE.GetAtomIndex( *itr_atoms_a);
            atom_index_d = MOLECULE.GetAtomIndex( *itr_atoms_d);
          }
        }
      }
      if( lowest_angle_measure > 180)
      {
        lowest_angle_measure = lowest_angle_measure - 360;
      }

      // if either central atom's bond angles have to be essentially linear and are not in a ring,
      // the dihedral angle will be ill-defined, so return 180
      const bool is_linear
      (
        ATOM_B.GetAtomType()->GetFormsOnlyLinearBonds()
        || ATOM_C.GetAtomType()->GetFormsOnlyLinearBonds()
      );
      if( is_linear && !IS_RING_BOND)
      {
        lowest_angle_measure = 180.0;
      }
      return storage::Pair< double, storage::VectorND< 4, size_t> >
             (
               lowest_angle_measure,
               storage::VectorND< 4, size_t>
               (
                 atom_index_a,
                 MOLECULE.GetAtomIndex( ATOM_B),
                 MOLECULE.GetAtomIndex( ATOM_C),
                 atom_index_d
               )
             );
    }

    //! @brief calculate prioritized dihedral angle
    //! @param MOLECULE the molecule of interest
    //! @param ATOM_A_PRIORITIES vector of indices of priority atoms connected to atom B
    //! @param ATOM_B atom involved in the central bond of the dihedral bond
    //! @param ATOM_C atom involved in the central bond of the dihedral bond
    //! @param ATOM_D_PRIORITIES vector of indices of priority atoms connected to atom C
    //! @param IS_RING_BOND true if B-C is in a ring
    storage::Vector< double> PriorityDihedralAngles::CalculateOrderedDihedralAngles
    (
      const ConformationInterface &MOLECULE,
      const storage::Vector< size_t> &ATOM_A_PRIORITIES,
      const AtomConformationalInterface &ATOM_B,
      const AtomConformationalInterface &ATOM_C,
      const storage::Vector< size_t> &ATOM_D_PRIORITIES,
      const bool &IS_RING_BOND
    ) const
    {
      // if either central atom's bond angles have to be essentially linear and are not in a ring,
      // the dihedral angle will be ill-defined, so return 180
      const bool is_linear
      (
        ATOM_B.GetAtomType()->GetFormsOnlyLinearBonds()
        || ATOM_C.GetAtomType()->GetFormsOnlyLinearBonds()
      );
      storage::Vector< double> ordered_dihedral_angles;
      // get iterators to the two central atoms involved in the dihedral bond
      iterate::Generic< const AtomConformationalInterface> itr_atoms_a( MOLECULE.GetAtomsIterator());
      iterate::Generic< const AtomConformationalInterface> itr_atoms_d( MOLECULE.GetAtomsIterator());

      // get measurement value for all possible dihedral angles about a central bond
      for
      (
        storage::Vector< size_t>::const_iterator itr_a( ATOM_A_PRIORITIES.Begin()), itr_a_end( ATOM_A_PRIORITIES.End());
        itr_a != itr_a_end;
        ++itr_a
      )
      {
        itr_atoms_a.GotoPosition( *itr_a);
        for
        (
          storage::Vector< size_t>::const_iterator itr_d( ATOM_D_PRIORITIES.Begin()), itr_d_end( ATOM_D_PRIORITIES.End());
          itr_d != itr_d_end;
          ++itr_d
        )
        {
          itr_atoms_d.GotoPosition( *itr_d);

          double angle_measure
          (
            math::Angle::Degree
            (
              linal::Dihedral
              (
                itr_atoms_a->GetPosition(),
                ATOM_B.GetPosition()      ,
                ATOM_C.GetPosition()      ,
                itr_atoms_d->GetPosition()
              )
            )
          );
          if( is_linear && !IS_RING_BOND)
          {
            angle_measure = 180.0;
          }
          ordered_dihedral_angles.PushBack( angle_measure);
        }
      }

      return ordered_dihedral_angles;
    }

    //! @brief calculate maximum dihedral angle difference, given ordered equivalent sets of dihedral angles
    //! @param DIHEDRALS_A dihedrals about a single bond, ordered by priority, for the first molecule
    //! @param DIHEDRALS_B dihedrals about a single bond, ordered by priority, for the second molecule
    double PriorityDihedralAngles::CalculateMaximumDihedralAngleDifference
    (
      const storage::Vector< double> &DIHEDRALS_A,
      const storage::Vector< double> &DIHEDRALS_B
    ) const
    {
      double maximal_difference( 0.0);
      for
      (
        storage::Vector< double>::const_iterator itr_a( DIHEDRALS_A.Begin()), itr_a_end( DIHEDRALS_A.End());
        itr_a != itr_a_end;
        ++itr_a
      )
      {
        double minimal_difference( 360.0);
        const double dihedral_a( *itr_a);
        for
        (
          storage::Vector< double>::const_iterator itr_b( DIHEDRALS_B.Begin()), itr_b_end( DIHEDRALS_B.End());
          itr_b != itr_b_end;
          ++itr_b
        )
        {
          const double difference( math::Absolute( dihedral_a - *itr_b));
          minimal_difference = std::min( minimal_difference, std::min( difference, 360.0 - difference));
        }
        maximal_difference = std::max( minimal_difference, maximal_difference);
      }
      return maximal_difference;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PriorityDihedralAngles::Read( std::istream &ISTREAM)
    {
      BCL_Exit( "Cannot read " + GetClassIdentifier(), -1);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &PriorityDihedralAngles::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      BCL_Exit( "Cannot write " + GetClassIdentifier(), -1);
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
