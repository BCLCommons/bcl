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
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_hydrogens_handler.h"
#include "chemistry/bcl_chemistry_substituent_conformational.h"
#include "graph/bcl_graph_edge_cover_ring_perception.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a small molecule to be standardized
    //! @param ATOMS the atoms to standardize
    //! @param ID identifier for the molecule
    //! @param FORCE_RECALCULATION true if all atom types should be recalculated
    AtomsCompleteStandardizer::AtomsCompleteStandardizer
    (
      AtomVector< AtomComplete> &ATOMS,
      const std::string &ID,
      const bool &FORCE_RECALCULATION
    ) :
      m_Atoms( RemoveObviousIonicBonds( ATOMS)),
      m_ID( ID),
      m_Graph // convert the molecule into a graph
      (
        storage::Vector< size_t>( ATOMS.GetSize(), 1),
        ATOMS.GetAdjacencyList( ConfigurationalBondTypeData::e_BondOrderOrAromatic),
        GetConfigurationalBondTypes().e_Undefined
      ),
      m_Rings( graph::EdgeCoverRingPerception( m_Graph).GetRings()),
      m_SmallestRingSize( ATOMS.GetSize()),
      m_NeedToRecomputeOtherAtoms( false)
    {
      AddRingInformationToBondTypes();

      // initialize atom type recalculation
      if( Initialize( FORCE_RECALCULATION))
      {
        // always detect aromatic ring systems
        SplitRingsByAromaticity();

        // set the atom types
        // only assign single-double bond orders if all atom types are valid
        SetAtomTypes();

        //BCL_MessageStd("enters SetConjugationOfBondTypes");
        // add information about whether a bond was conjugated to the bond types
        SetConjugationOfBondTypes();

        // assign single-double bond orders to aromatic systems
        DetermineUnknownBondOrders();
      }
    }

    //! @brief Clone function
    //! @return pointer to new AtomsCompleteStandardizer
    AtomsCompleteStandardizer *AtomsCompleteStandardizer::Clone() const
    {
      return new AtomsCompleteStandardizer( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AtomsCompleteStandardizer::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AtomsCompleteStandardizer::Read( std::istream &ISTREAM)
    {
      BCL_Exit( "Small molecule standardizer cannot be read because it contains references!", -1);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &AtomsCompleteStandardizer::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
       io::Serialize::Write( m_Graph, OSTREAM, INDENT);
       io::Serialize::Write( m_Rings, OSTREAM, INDENT);
       io::Serialize::Write( m_AromaticRings, OSTREAM, INDENT);
       io::Serialize::Write( m_ConjugatedRings, OSTREAM, INDENT);
       io::Serialize::Write( m_PossibleAtomTypes, OSTREAM, INDENT);

      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief add ring information to the bond types
    void AtomsCompleteStandardizer::AddRingInformationToBondTypes()
    {
      m_NumberRings = storage::Vector< size_t>( m_Atoms.GetSize(), size_t( 0));
      m_NumberDeclaredAromaticBonds = storage::Vector< size_t>( m_Atoms.GetSize(), size_t( 0));
      m_SmallestRingSize.SetAllElements( util::GetUndefined< size_t>());

      {
        // count aromatic bonds for each atom
        size_t index( 0);
        for
        (
          AtomVector< AtomComplete>::const_iterator itr( m_Atoms.Begin()), itr_end( m_Atoms.End());
          itr != itr_end;
          ++itr, ++index
        )
        {
          m_NumberDeclaredAromaticBonds( index) =
            itr->CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsAromatic, size_t( 1));
        }
      }

      for
      (
        storage::List< graph::Ring>::const_iterator itr_ring( m_Rings.Begin()), itr_ring_end( m_Rings.End());
        itr_ring != itr_ring_end;
        ++itr_ring
      )
      {
        // record the index of the previous/last atom
        size_t prev_index( itr_ring->LastElement());

        if( RingHasAllAromaticBondTypes( *itr_ring))
        {
          for( graph::Ring::const_iterator itr( itr_ring->Begin()), itr_end( itr_ring->End()); itr != itr_end; ++itr)
          {
            ++m_NumberRings( *itr);
          }
        }
        else
        {
          // iterate over all pairs of consecutive atoms in the ring
          for( graph::Ring::const_iterator itr( itr_ring->Begin()), itr_end( itr_ring->End()); itr != itr_end; ++itr)
          {
            ++m_NumberRings( *itr);

            // get the old bond type between these atoms
            const ConfigurationalBondType &old_bond_type( m_Atoms( prev_index).GetBondTypeTo( m_Atoms( *itr)));

            // get the bond type between these atoms and find the corresponding bond type in a ring
            const ConfigurationalBondType &new_bond_type( old_bond_type->WithInRing());
            if( old_bond_type != new_bond_type)
            {
              // set the bond type up bidirectionally
              m_Atoms( prev_index).SetBondTypeTo( m_Atoms( *itr), new_bond_type);
            }

            // update smallest ring size
            m_SmallestRingSize( *itr) = std::min( m_SmallestRingSize( *itr), itr_ring->GetSize());

            // update previous index
            prev_index = *itr;
          }
        }
        // update smallest ring size
        for( graph::Ring::const_iterator itr( itr_ring->Begin()), itr_end( itr_ring->End()); itr != itr_end; ++itr)
        {
          m_SmallestRingSize( *itr) = std::min( m_SmallestRingSize( *itr), itr_ring->GetSize());
        }
      }

      {
        // Find "Aromatic" bonds outside rings, change them to Conjugated instead (sdf type 5)
        size_t index( 0);
        for
        (
          AtomVector< AtomComplete>::iterator itr( m_Atoms.Begin()), itr_end( m_Atoms.End());
          itr != itr_end;
          ++itr, ++index
        )
        {
          if( m_NumberDeclaredAromaticBonds( index) && !util::IsDefined( m_SmallestRingSize( index)))
          {
            BCL_MessageCrt( "Warning: Converting Aromatic Bonds that were not in rings to conjugated types!");
            // translate any existing constitutional single bonds to the non-conjugated state
            for
            (
              storage::Vector< BondConformational>::const_iterator
                itr_bond( itr->GetBonds().Begin()), itr_bond_end( itr->GetBonds().End());
              itr_bond != itr_bond_end;
              ++itr_bond
            )
            {
              if( itr_bond->GetBondType() == GetConfigurationalBondTypes().e_AromaticBond)
              {
                itr->SetBondTypeMonoDirectional
                (
                  itr_bond->GetTargetAtom(),
                  GetConfigurationalBondTypes().e_ConjugatedBond
                );
              }
            }
            m_NumberDeclaredAromaticBonds( index) = 0;
          }
        }
      }
    }

    //! @brief initialize the standardization
    //! @param FORCE whether the force recalculation of atom typesConstitutionalBondTypeData::e_Conjugated
    //! @return true if all types could be assigned an atom type
    bool AtomsCompleteStandardizer::Initialize( const bool &FORCE)
    {
      // clean the vector of possible atom types
      m_PossibleAtomTypes.Reset();
      m_PossibleAtomTypes.AllocateMemory( m_Atoms.GetSize());

      // track whether all atoms have at least one possible type
      bool atom_types_are_defined( true);
      for
      (
        AtomVector< AtomComplete>::iterator itr_atom( m_Atoms.Begin()), itr_atom_end( m_Atoms.End());
        itr_atom != itr_atom_end;
        ++itr_atom
      )
      {
        // add a new possible type object for this atom
        m_PossibleAtomTypes.PushBack();

        // get a reference to the added possible atom type
        PossibleAtomTypesForAtom &possible_types( m_PossibleAtomTypes.LastElement());

        // if the type was known, add the known value to possible_types and continue
        if( !FORCE && itr_atom->GetAtomType()->IsGasteigerAtomType())
        {
          possible_types.SetToType( itr_atom->GetAtomType());
        }
        else
        {
          // remove any existing gasteiger atom type
          itr_atom->SetAtomType( itr_atom->GetAtomType()->GetBaseAtomType());
          // translate any existing constitutional single bonds to the non-conjugated state
          for
          (
            storage::Vector< BondConformational>::const_iterator
              itr_bond( itr_atom->GetBonds().Begin()), itr_bond_end( itr_atom->GetBonds().End());
            itr_bond != itr_bond_end;
            ++itr_bond
          )
          {
            if
            (
              itr_bond->GetBondType()->GetBondData( ConfigurationalBondTypeData::e_BondOrder) == size_t( 1)
              &&
              (
                itr_bond->GetBondType()->GetConjugation() == ConstitutionalBondTypeData::e_Conjugated
                || itr_bond->GetBondType()->GetConjugation() == ConstitutionalBondTypeData::e_Amide
              )
            )
            {
              itr_atom->SetBondTypeMonoDirectional
              (
                itr_bond->GetTargetAtom(),
                itr_bond->GetBondType()->WithConjugation( ConstitutionalBondTypeData::e_Nonconjugated)
              );
            }
            else if
            (
              itr_bond->GetBondType()->IsBondOrderKnown()
              && itr_bond->GetBondType()->GetConjugation() == ConstitutionalBondTypeData::e_Aromatic
            )
            {
              itr_atom->SetBondTypeMonoDirectional
              (
                itr_bond->GetTargetAtom(),
                itr_bond->GetBondType()->GetBondData( ConfigurationalBondTypeData::e_BondOrder) == size_t( 1)
                ? itr_bond->GetBondType()->WithConjugation( ConstitutionalBondTypeData::e_Nonconjugated)
                : itr_bond->GetBondType()->WithConjugation( ConstitutionalBondTypeData::e_Conjugated)
              );
            }
          }
          // determine initially the possible types for the atom and its bonds
          possible_types = GetPossibleTypesForAtom( *itr_atom);
        }

        // track whether any atom cannot be assigned a possible atom type
        atom_types_are_defined = atom_types_are_defined && possible_types.GetMostStableType().IsDefined();
      }
      if( m_NeedToRecomputeOtherAtoms)
      {
        m_NeedToRecomputeOtherAtoms = false;
        return Initialize( FORCE);
      }

      return atom_types_are_defined;
    }

    //! @brief split rings by conjugation (non-conjugated, conjugated, aromatic)
    void AtomsCompleteStandardizer::SplitRingsByAromaticity()
    {
      storage::Vector< size_t> is_conjugatable_ring_atom( m_Atoms.GetSize(), size_t( 0));

      // whether the atom is an oxygen or nitrogen with a double bond to an atom in a ring
      m_IsONDoubleBondedToRing = storage::Vector< size_t>( m_Atoms.GetSize(), size_t( 0));

      // vector of all 1s for each atom
      storage::Vector< size_t> ones( m_Atoms.GetSize(), size_t( 1));

      // make a vector of indices of conjugated atoms
      storage::Vector< size_t> conjugated_atom_vector;
      conjugated_atom_vector.AllocateMemory( m_Atoms.GetSize());

      // track the # of non-conjugated atoms in rings
      size_t non_conjugated_atoms_in_rings_count( 0);
      for( size_t atom_index( 0), n_atoms( m_Atoms.GetSize()); atom_index < n_atoms; ++atom_index)
      {
        if( m_NumberDeclaredAromaticBonds( atom_index))
        {
          is_conjugatable_ring_atom( atom_index) = true;
          conjugated_atom_vector.PushBack( atom_index);
        }
        else if( m_NumberRings( atom_index))
        {
          is_conjugatable_ring_atom( atom_index) = m_PossibleAtomTypes( atom_index).CouldBeConjugated();
          if( !is_conjugatable_ring_atom( atom_index))
          {
            ++non_conjugated_atoms_in_rings_count;
          }
          else
          {
            conjugated_atom_vector.PushBack( atom_index);
          }
        }
        else if
        (
          m_Atoms( atom_index).GetElementType() == GetElementTypes().e_Nitrogen
          || m_Atoms( atom_index).GetElementType() == GetElementTypes().e_Oxygen
        )
        {
          // check for double bond to an atom in a ring that is less electronegative (e.g. N=N doesn't apply).
          // The pi electrons in the double bond will localize near these highly electronegative species, and therefore
          // not break a pi-system, though they also do not donate electrons to the pi-system
          const storage::Vector< BondConformational> &bonds( m_Atoms( atom_index).GetBonds());
          for
          (
            storage::Vector< BondConformational>::const_iterator itr( bonds.Begin()), itr_end( bonds.End());
            itr != itr_end;
            ++itr
          )
          {
            if
            (
              itr->GetBondType()->GetNumberOfElectrons() == size_t( 4)
              && m_NumberRings( m_Atoms.GetAtomIndex( itr->GetTargetAtom()))
              // do not count tetravalent atoms with double bonds outside the ring, even if they are to N or O
              && itr->GetTargetAtom().GetBonds().GetSize() != size_t( 4)
            )
            {
              m_IsONDoubleBondedToRing( atom_index) =
                itr->GetTargetAtom().GetElementType()->GetProperty( ElementTypeData::e_ElectroNegativity)
                < m_Atoms( atom_index).GetElementType()->GetProperty( ElementTypeData::e_ElectroNegativity);
              is_conjugatable_ring_atom( atom_index) = !m_IsONDoubleBondedToRing( atom_index);
              break;
            }
          }
        }
      }

      if( non_conjugated_atoms_in_rings_count)
      {
        // detect additional conjugated rings where a conjugated ring is hidden by the presence of a non-conjugated bridge
        // get possible rings from list of conjugated atoms outside aromatic rings
        m_ConjugatedRings =
          graph::EdgeCoverRingPerception( m_Graph.GetSubgraph( conjugated_atom_vector)).GetRings();

        // translate the atom indices back into the whole-molecule frame-of-reference
        for
        (
          storage::List< graph::Ring>::iterator itr_ring( m_ConjugatedRings.Begin()),
                                                itr_ring_end( m_ConjugatedRings.End());
          itr_ring != itr_ring_end;
          ++itr_ring
        )
        {
          storage::Vector< size_t> translated_ring_indices;
          translated_ring_indices.AllocateMemory( itr_ring->GetSize());
          for
          (
            graph::Ring::const_iterator itr( itr_ring->Begin()), itr_end( itr_ring->End());
            itr != itr_end;
            ++itr
          )
          {
            translated_ring_indices.PushBack( conjugated_atom_vector( *itr));
          }
          *itr_ring = graph::Ring( m_Atoms.GetSize(), translated_ring_indices);
        }
      }
      else
      {
        // all atoms in rings are conjugated (a common case)
        m_ConjugatedRings = m_Rings;
      }

      for
      (
        storage::List< graph::Ring>::iterator itr_ring( m_ConjugatedRings.Begin()),
                                              itr_ring_end( m_ConjugatedRings.End());
        itr_ring != itr_ring_end;
        // iteration in loop
      )
      {
        storage::List< graph::Ring>::iterator old_ring_itr( itr_ring);
        ++itr_ring;
        if( IsAromatic( *old_ring_itr, is_conjugatable_ring_atom))
        {
          SetAromaticBondTypes( *old_ring_itr);
          m_AromaticRings.InternalData().splice( m_AromaticRings.End(), m_ConjugatedRings.InternalData(), old_ring_itr);
        }
      }

      const size_t number_basic_aromatic_rings( m_AromaticRings.GetSize());

      storage::Set< util::SiPtr< const graph::Ring> > problematic_rings;

      // fuse conjugated rings to check for aromaticity e.g. azulene
      size_t itr_a_id( 0);
      for
      (
        storage::List< graph::Ring>::iterator itr_a( m_ConjugatedRings.Begin());
        itr_a != m_ConjugatedRings.End();
        ++itr_a, ++itr_a_id
      )
      {
        // try combining with all other conjugated rings
        bool fusing_produced_aromatic_ring( false);

        if
        (
          !RingHasAllAromaticBondTypes( *itr_a)
          && TestBondTypesAlreadySetToAromatic( *itr_a)
        )
        {
          storage::List< graph::Ring>::iterator old_ring_itr( itr_a);
          --itr_a;
          m_AromaticRings.InternalData().splice( m_AromaticRings.End(), m_ConjugatedRings.InternalData(), old_ring_itr);
          continue;
        }

        size_t itr_b_id( 0);
        for
        (
          storage::List< graph::Ring>::iterator itr_b( m_ConjugatedRings.Begin());
          itr_b != itr_a;
          ++itr_b, ++itr_b_id
        )
        {
          size_t nfs( graph::Ring::GetNominalFusionSize( *itr_a, *itr_b));
          if( nfs <= std::max( itr_a->GetSize(), itr_b->GetSize()))
          {
            continue;
          }
          graph::Ring fused_path( graph::Ring::FuseRings( *itr_a, *itr_b));
          if( !fused_path.GetSize())
          {
            continue;
          }
          if( TestBondTypesAlreadySetToAromatic( fused_path) && !RingHasAllAromaticBondTypes( fused_path))
          {
            fusing_produced_aromatic_ring = true;
            continue;
          }
          // if the rings could not be fused, the resulting ring is empty, empty rings do not satisfy hueckel's rule (4n+2)
          if( IsAromatic( fused_path, is_conjugatable_ring_atom))
          {
            m_AromaticRings.PushBack( fused_path);
            SetAromaticBondTypes( fused_path);
            fusing_produced_aromatic_ring = true;
            storage::Set< util::SiPtr< const graph::Ring> >::iterator
              itr_problematic_ring( problematic_rings.Find( *itr_b));
            if( itr_problematic_ring != problematic_rings.End())
            {
              problematic_rings.RemoveElement( itr_problematic_ring);
            }
            if
            (
              !RingHasAllAromaticBondTypes( *itr_a)
              && TestBondTypesAlreadySetToAromatic( *itr_a)
            )
            {
              storage::List< graph::Ring>::iterator old_ring_itr( itr_a);
              --itr_a;
              m_ConjugatedRings.InternalData().erase( old_ring_itr);
              break;
            }
          }
          else
          {
            m_ConjugatedRings.PushBack( fused_path);
          }
        }

        // try combining with all the original aromatic rings, if combining with other conjugated rings failed
        if( fusing_produced_aromatic_ring)
        {
          continue;
        }
        size_t ring_count_b( 0);
        for
        (
          storage::List< graph::Ring>::const_iterator
            itr_b( m_AromaticRings.Begin()), itr_b_end( m_AromaticRings.End());
          ring_count_b < number_basic_aromatic_rings;
          ++itr_b, ++ring_count_b
        )
        {
          size_t nfs( graph::Ring::GetNominalFusionSize( *itr_a, *itr_b));
          if( nfs <= std::max( itr_a->GetSize(), itr_b->GetSize()))
          {
            continue;
          }
          graph::Ring fused_path( graph::Ring::FuseRings( *itr_a, *itr_b));
          // if the rings could not be fused, the resulting ring is empty, empty rings do not satisfy hueckel's rule (4n+2)
          if( IsAromatic( fused_path, is_conjugatable_ring_atom))
          {
            if( !TestBondTypesAlreadySetToAromatic( fused_path) || RingHasAllAromaticBondTypes( *itr_a))
            {
              SetAromaticBondTypes( fused_path);
              m_AromaticRings.PushBack( fused_path);
            }
            fusing_produced_aromatic_ring = true;
          }
        }

        if
        (
          !fusing_produced_aromatic_ring
          && RingHasAllAromaticBondTypes( *itr_a)
          && !IsAromatic( *itr_a, ones)
        )
        {
          problematic_rings.Insert( *itr_a);
        }
      }
      for
      (
        storage::Set< util::SiPtr< const graph::Ring> >::const_iterator
          itr_ring( problematic_rings.Begin()), itr_ring_end( problematic_rings.End());
        itr_ring != itr_ring_end;
        ++itr_ring
      )
      {
        RefineAtomTypesForAromaticRing( **itr_ring);
        if( IsAromatic( **itr_ring, is_conjugatable_ring_atom))
        {
          m_AromaticRings.PushBack( **itr_ring);
        }
        else if
        (
          !IsAromatic( **itr_ring, ones)
          && ( CountNominalEInPiSystem( **itr_ring) % 2) == size_t( 1)
        )
        {
          // run over the atoms; find any atoms next to a charged atom
          storage::Map< std::string, size_t> atom_type_count;
          for
          (
            graph::Ring::const_iterator itr( ( *itr_ring)->Begin()), itr_end( --( *itr_ring)->End());
            itr != itr_end;
            ++itr
          )
          {
            atom_type_count[ m_PossibleAtomTypes( *itr).GetMostStableType()]++;
          }
          // if so, then it could be that the initial choices for possible atom types are incorrect
          std::ostringstream msg;
          msg << "Molecule " << m_ID << " aromaticity detection disagreed with input file! # pi e- == "
              << CountNominalEInPiSystem( **itr_ring) << " but ring was declared aromatic nominal types: ";
          for
          (
            storage::Map< std::string, size_t>::const_iterator
              itr( atom_type_count.Begin()), itr_end( atom_type_count.End());
            itr != itr_end;
            ++itr
          )
          {
            msg << itr->first << "(" << itr->second << ") ";
          }
          BCL_MessageStd( msg.str());
        }
      }

      // set the bond types in all aromatic rings
      for
      (
        storage::List< graph::Ring>::iterator itr_ring( m_AromaticRings.Begin()), itr_ring_end( m_AromaticRings.End());
        itr_ring != itr_ring_end;
        ++itr_ring
      )
      {
        SetAromaticBondTypes( *itr_ring);
      }
    }

    //! @brief Check a ring for aromaticity
    //! @param RING ring which has to be detected for aromaticity
    //! @param IS_CONJUGATED_ATOM_IN_RING  0 for atoms that are not conjugated atoms in rings, 1 otherwise
    //! @return true if ring is aromatic; false otherwise
    bool AtomsCompleteStandardizer::IsAromatic
    (
      const graph::Ring &RING,
      const storage::Vector< size_t> &IS_CONJUGATED_ATOM_IN_RING
    )
    {
      if( RING.GetSize() < size_t( 3))
      {
        return false;
      }

      // get the max # of pi electrons in the pi system
      const size_t max_pi_electrons( CountNominalEInPiSystem( RING));

      // take #pi electrons mod 4, to see if pi-electrons == 4 * n + 2
      const size_t max_pi_electrons_mod_4( ( max_pi_electrons) % size_t( 4));

      if( max_pi_electrons_mod_4 == size_t( 2)) // does the ring satisfy hueckel's rule?
      {
        // check that the number of pi electrons is not excessive. Rings with pi-electron excess do not appear to exhibit
        // aromatic stabilization, e.g. compare 1-2 dioxin to the chain structure. Ditto for tetrazine. See discussion on
        // https://www.chemaxon.com/forum/ftopic319.html. The note about ring systems does not apply for the bcl, since
        // this function only considers individual rings, possibly those fused from a larger ring system.
        if( max_pi_electrons >= RING.GetSize() + 2)
        {
          return false;
        }
        // Yes, so check that the pi-ring system is uninterrupted by exocyclic double bonds to atoms that
        // are not part of any conjugated ring system and are not to highly-electronegative species such as O or N
        // Such bonds cause, e.g. 2-pyridone, the oxygen-containing ring of coumarin, and 2,4-cyclopentadiene-1-one
        // to be non-aromatic
        if( DoubleBondsAreToAtomsInSet( RING, IS_CONJUGATED_ATOM_IN_RING))
        {
          return true;
        }
      }
      return false;
    }

    //! @brief test whether a given ring contains an atom with a lone pair
    //! @param RING the path to examine
    //! @return the number of atoms that can donate a lone pair to the pi system
    size_t AtomsCompleteStandardizer::CountAtomsWithLonePairForPiSystem( const graph::Ring &RING)
    {
      size_t lone_pair_count( 0);
      for
      (
        graph::Ring::const_iterator itr( RING.Begin()), itr_end( RING.End());
        itr != itr_end;
        ++itr
      )
      {
        if
        (
          m_PossibleAtomTypes( *itr).GetMostStableType()->GetPiElectronContributionType() == AtomTypeData::e_ZeroOrTwo
        )
        {
          ++lone_pair_count;
        }
      }
      return lone_pair_count;
    }

    //! Set bond types of all bonds in a ring to include aromatic character
    void AtomsCompleteStandardizer::SetAromaticBondTypes( const graph::Ring &RING)
    {
      // record the index of the previous/last atom
      size_t prev_index( RING.LastElement());

      // iterate over all pairs of consecutive atoms in the ring
      for( graph::Ring::const_iterator itr( RING.Begin()), itr_end( RING.End()); itr != itr_end; ++itr)
      {
        // get the current bond type
        const ConfigurationalBondType &bond_type( m_Atoms( *itr).GetBondTypeTo( m_Atoms( prev_index)));

        // set the bond type bidirectionally to be aromatic
        m_Atoms( *itr).SetBondTypeTo
        (
          m_Atoms( prev_index),
          bond_type->WithConjugation( ConstitutionalBondTypeData::e_Aromatic)
        );

        prev_index = *itr;
      }
    }

    //! Test whether all the bond types in a ring were already set to being aromatic
    bool AtomsCompleteStandardizer::TestBondTypesAlreadySetToAromatic( const graph::Ring &RING) const
    {
      // record the index of the previous/last atom
      size_t prev_index( RING.LastElement());

      // iterate over all pairs of consecutive atoms in the ring
      for( graph::Ring::const_iterator itr( RING.Begin()), itr_end( RING.End()); itr != itr_end; ++itr)
      {
        if( m_Atoms( *itr).GetBondTypeTo( m_Atoms( prev_index))->GetConjugation() != ConstitutionalBondTypeData::e_Aromatic)
        {
          return false;
        }
        prev_index = *itr;
      }
      return true;
    }

    //! @brief test whether a particular path contains any non-conjugatable atoms
    //! @param RING the path to check
    //! @return true if the path contains any non-conjugatable atoms
    bool AtomsCompleteStandardizer::RingContainsNonConjugableAtoms( const graph::Ring &RING)
    {
      for
      (
        graph::Ring::const_iterator itr( RING.Begin()), itr_end( RING.End());
        itr != itr_end;
        ++itr
      )
      {
        if( !m_PossibleAtomTypes( *itr).CouldBeConjugated())
        {
          return true;
        }
      }
      return false;
    }

    //! @brief test whether any atoms in a path contain double bonds to atoms outside the conjugated atoms list
    //! @param RING the path whose atoms to check
    //! @param IS_CONJUGATED_ATOM_IN_RING 0 for atoms that are not conjugated atoms in rings
    bool AtomsCompleteStandardizer::DoubleBondsAreToAtomsInSet
    (
      const graph::Ring &RING,
      const storage::Vector< size_t> &IS_CONJUGATED_ATOM_IN_RING
    ) const
    {
      static const size_t double_bond_type( 2);
      for
      (
        graph::Ring::const_iterator itr( RING.Begin()), itr_end( RING.End());
        itr != itr_end;
        ++itr
      )
      {
        // check for conjugated valences
        if( m_Atoms( *itr).GetNumberElectronsInValenceBonds() > m_Atoms( *itr).GetNumberValenceBonds())
        {
          return false;
        }

        for
        (
          graph::ConstGraph< size_t, size_t>::t_EdgeTargetsOfVertex::const_iterator
            itr_neighbor( m_Graph.GetNeighborIndices( *itr).Begin()),
            itr_neighbor_end( m_Graph.GetNeighborIndices( *itr).End());
          itr_neighbor != itr_neighbor_end;
          ++itr_neighbor
        )
        {
          if
          (
            m_Graph.GetEdgeData( *itr, *itr_neighbor) == double_bond_type
            && !IS_CONJUGATED_ATOM_IN_RING( *itr_neighbor)
          )
          {
            return false;
          }
        }
      }
      return true;
    }

    //! @brief test whether a ring appears to have purely aromatic bond types
    //! @param RING the path whose atoms to check
    bool AtomsCompleteStandardizer::RingHasAllAromaticBondTypes( const graph::Ring &RING) const
    {
      static const size_t aromatic_bond_type( 4);
      size_t prev_index( RING.LastElement());
      for
      (
        graph::Ring::const_iterator itr( RING.Begin()), itr_end( RING.End());
        itr != itr_end;
        ++itr
      )
      {
        if( m_Graph.GetEdgeData( *itr, prev_index) != aromatic_bond_type)
        {
          return false;
        }
        prev_index = *itr;
      }
      return true;
    }

    //! @brief refine atom types for a ring with all aromatic declared bonds
    //! @param RING the ring, for which RingHasAllAromaticBondTypes must have returned true
    void AtomsCompleteStandardizer::RefineAtomTypesForAromaticRing( const graph::Ring &RING)
    {
      // take #pi electrons mod 4, to see if pi-electrons == 4 * n + 2
      const size_t max_pi_electrons_mod_4( CountNominalEInPiSystem( RING) % size_t( 4));

      // if pi electrons modulo 4 is already 0 or 2, then the IsAromatic function will sort it out
      if( max_pi_electrons_mod_4 == size_t( 2) || max_pi_electrons_mod_4 == size_t( 0))
      {
        return;
      }

      // Consider all ways of getting to 2 pi electrons while only creating one charge
      //  From 3 pi electrons:
      //    Remove a pi-electron from bonds
      //    Lone pair demotion
      //    If lone pair present:
      //      Add pi electron (to form a 2x bond)
      //  From 1 pi electron:
      //    Lone pair creation
      //    Add pi electron (to form a 2x bond)
      //    If lone pair present:
      //      Lone pair demotion (if at least one other lone pair is present)
      //      Remove a pi-electron from bonds

      // iterate over ring,
      // find any atom with an alternative atom type that is connected to an atom type with a charge

      util::SiPtr< AtomComplete> chosen_atom;
      AtomType chosen_atom_type;
      util::SiPtr< PossibleAtomTypesForAtom> chosen_possible_atom_type;

      for
      (
        graph::Ring::const_iterator itr( RING.Begin()), itr_end( RING.End());
        itr != itr_end;
        ++itr
      )
      {
        AtomComplete &atom( m_Atoms( *itr));
        PossibleAtomTypesForAtom &possible_atom_types( m_PossibleAtomTypes( *itr));
        const AtomType &current_pick( possible_atom_types.GetMostStableType());
        if
        (
          possible_atom_types.GetMostStableType()->GetFormalCharge() != short( 0)
          || possible_atom_types.GetNumberPossibleTypes() <= size_t( 1)
          || m_NumberDeclaredAromaticBonds( *itr) >= size_t( 3)
        )
        {
          continue;
        }

        // get the bonds of this atom
        const storage::Vector< BondConformational> &bonds( atom.GetBonds());

        for
        (
          storage::Vector< BondConformational>::const_iterator
            itr_bonds( bonds.Begin()), itr_bonds_end( bonds.End());
          itr_bonds != itr_bonds_end;
          ++itr_bonds
        )
        {
          AtomComplete &target_atom( m_Atoms( itr_bonds->GetTargetAtom()));
          const short target_atom_charge
          (
            m_PossibleAtomTypes( m_Atoms.GetAtomIndex( target_atom)).GetMostStableType()->GetFormalCharge()
          );
          if( target_atom_charge == short( 0))
          {
            continue;
          }

          AtomType alternate_type( possible_atom_types.GetAlternateTypeWithCharge( -target_atom_charge));
          if( !alternate_type.IsDefined())
          {
            continue;
          }

          if( current_pick->GetMaxEContributionToPiSystem() > alternate_type->GetMaxEContributionToPiSystem())
          {
            if( current_pick->GetMaxEContributionToPiSystem() == alternate_type->GetMaxEContributionToPiSystem() + 1)
            {
              if( max_pi_electrons_mod_4 == size_t( 1))
              {
                continue;
              }
            }
            else
            {
              continue;
            }
          }
          else if( current_pick->GetMaxEContributionToPiSystem() + 1 == alternate_type->GetMaxEContributionToPiSystem())
          {
            if( max_pi_electrons_mod_4 == size_t( 3))
            {
              continue;
            }
          }
          else
          {
            continue;
          }

          if
          (
            !chosen_atom.IsDefined()
            ||
            (
              &*chosen_atom == &atom
              && alternate_type->GetStabilityMetric() > chosen_atom_type->GetStabilityMetric()
            )
            ||
            (
              &*chosen_atom != &atom
              && SubstituentConformational( atom) < SubstituentConformational( *chosen_atom)
            )
          )
          {
            chosen_atom = atom;
            chosen_possible_atom_type = possible_atom_types;
            chosen_atom_type = alternate_type;
          }
        }
      }

      if( chosen_atom.IsDefined())
      {
        // if pi electrons modulo 4 is already 0 or 2, then the IsAromatic function will sort it out
        BCL_MessageStd
        (
          "Changed nominal type: " + util::Format()( chosen_possible_atom_type->GetMostStableType().GetName())
          + " into " + util::Format()( chosen_atom_type.GetName()) + " to fix aromaticity"
        );
        chosen_possible_atom_type->SetToType( chosen_atom_type);
        return;
      }

      for
      (
        graph::Ring::const_iterator itr( RING.Begin()), itr_end( RING.End());
        itr != itr_end;
        ++itr
      )
      {
        AtomComplete &atom( m_Atoms( *itr));
        PossibleAtomTypesForAtom &possible_atom_types( m_PossibleAtomTypes( *itr));
        const AtomType &current_pick( possible_atom_types.GetMostStableType());
        if
        (
          possible_atom_types.GetMostStableType()->GetFormalCharge() != short( 0)
          || possible_atom_types.GetNumberPossibleTypes() <= size_t( 1)
          || m_NumberDeclaredAromaticBonds( *itr) >= size_t( 3)
        )
        {
          continue;
        }

        storage::Vector< AtomType> alternate_types( possible_atom_types.GetAlternateTypes());

        for
        (
          storage::Vector< AtomType>::const_iterator
            itr_alternate( alternate_types.Begin()), itr_alternate_end( alternate_types.End());
          itr_alternate != itr_alternate_end;
          ++itr_alternate
        )
        {
          AtomType alternate_type( *itr_alternate);
          if( current_pick->GetMaxEContributionToPiSystem() > alternate_type->GetMaxEContributionToPiSystem())
          {
            if( current_pick->GetMaxEContributionToPiSystem() == alternate_type->GetMaxEContributionToPiSystem() + 1)
            {
              if( max_pi_electrons_mod_4 == size_t( 1))
              {
                continue;
              }
            }
            else
            {
              continue;
            }
          }
          else if( current_pick->GetMaxEContributionToPiSystem() + 1 == alternate_type->GetMaxEContributionToPiSystem())
          {
            if( max_pi_electrons_mod_4 == size_t( 3))
            {
              continue;
            }
          }
          else
          {
            continue;
          }

          if
          (
            !chosen_atom.IsDefined()
            ||
            (
              &*chosen_atom != &atom
              && SubstituentConformational( atom) < SubstituentConformational( *chosen_atom)
            )
          )
          {
            chosen_atom = atom;
            chosen_atom_type = alternate_type;
            chosen_possible_atom_type = possible_atom_types;
          }
        }
      }

      if( chosen_atom.IsDefined())
      {
        // if pi electrons modulo 4 is already 0 or 2, then the IsAromatic function will sort it out
        BCL_MessageStd
        (
          "Changed nominal type: " + util::Format()( chosen_possible_atom_type->GetMostStableType().GetName())
          + " into " + util::Format()( chosen_atom_type.GetName()) + " to fix aromaticity"
        );
        chosen_possible_atom_type->SetToType( chosen_atom_type);
        return;
      }
    }

    //! @brief test whether a ring appears to have purely aromatic bond types
    //! @param RING the path whose atoms to check
    bool AtomsCompleteStandardizer::CheckAllAtomTypeValid( const graph::Ring &RING)
    {
      for
      (
        graph::Ring::const_iterator itr( RING.Begin()), itr_end( RING.End());
        itr != itr_end;
        ++itr
      )
      {
        if( m_PossibleAtomTypes( *itr).GetNumberPossibleTypes() == size_t( 0))
        {
          return false;
        }
      }
      return true;
    }

    //! @brief count the electrons in a pi-system along atoms in RING
    //! @param RING the path to examine
    //! @return the number of electrons in a pi-system along atoms in RING
    size_t AtomsCompleteStandardizer::CountNominalEInPiSystem
    (
      const graph::Ring &RING
    )
    {
      size_t count( 0);
      for( graph::Ring::const_iterator itr( RING.Begin()), itr_end( RING.End()); itr != itr_end; ++itr)
      {
        count += m_PossibleAtomTypes( *itr).GetMaxElectronsParticipatingInPiSystem();
        size_t neigh( 0);
        for
        (
          storage::Vector< size_t>::const_iterator
            itr_neighbor( m_Graph.GetNeighborIndices( *itr).Begin()),
            itr_neighbor_end( m_Graph.GetNeighborIndices( *itr).End());
          itr_neighbor != itr_neighbor_end;
          ++itr_neighbor, ++neigh
        )
        {
          if( m_IsONDoubleBondedToRing( *itr_neighbor) && m_Graph.GetNeighborData( *itr)( neigh) == size_t( 2))
          {
            --count;
          }
        }
      }
      return count;
    }

    namespace
    {
      //! @brief test whether a particular bond of this atoms is an amide bond (which are generally non-rotatable)
      //! @param BOND bond to test
      //! @return true if the bond is an amide bond
      //! We base our definition of amide bonds on an internal investigation where we identified the dihedral bond
      //! tuples that are generally guaranteed to be planar in crystallographic structures, primarily from the
      //! crystallographic open database, with the only real constraint being that one of the terminii of the dihedral bond
      //! had to be a double bond to an oxygen or sulfur
      bool IsNonRotatableAmideLikeBond( const AtomConformationalInterface &ATOM, const BondConformational &BOND)
      {
        if( BOND.GetBondType()->GetNumberOfElectrons() != size_t( 2))
        {
          return false;
        }

        const AtomConformationalInterface &target_atom( BOND.GetTargetAtom());
        if( ATOM.GetAtomType()->GetNumberBonds() != size_t( 3) || target_atom.GetAtomType()->GetNumberBonds() != size_t( 3))
        {
          return false;
        }

        if
        (
          ATOM.GetAtomType()->GetHybridOrbitalType() != GetHybridOrbitalTypes().e_SP3
          && target_atom.GetElementType() == GetElementTypes().e_Nitrogen
          && target_atom.GetAtomType()->GetHybridOrbitalType() == GetHybridOrbitalTypes().e_SP2
          // && target_atom.GetAtomType()->GetNumberBonds() == 3 (already checked)
        )
        {
          for
          (
            storage::Vector< BondConformational>::const_iterator
              itr_bonds_c( ATOM.GetBonds().Begin()), itr_bonds_c_end( ATOM.GetBonds().End());
            itr_bonds_c != itr_bonds_c_end;
            ++itr_bonds_c
          )
          {
            // check for a double bond
            if
            (
              itr_bonds_c->GetTargetAtom().GetAtomType() == GetAtomTypes().O_Tr2Tr2TrPi
              || itr_bonds_c->GetTargetAtom().GetAtomType() == GetAtomTypes().S_Tr2Tr2TrPi
            )
            {
              return true;
            }
          }
        }

        if
        (
          target_atom.GetAtomType()->GetHybridOrbitalType() != GetHybridOrbitalTypes().e_SP3
          && ATOM.GetElementType() == GetElementTypes().e_Nitrogen
          && ATOM.GetAtomType()->GetHybridOrbitalType() == GetHybridOrbitalTypes().e_SP2
          // && ATOM.GetAtomType()->GetNumbe
        )
        {
          for
          (
            storage::Vector< BondConformational>::const_iterator
              itr_bonds_c( target_atom.GetBonds().Begin()), itr_bonds_c_end( target_atom.GetBonds().End());
            itr_bonds_c != itr_bonds_c_end;
            ++itr_bonds_c
          )
          {
            // check for a double bond
            if
            (
              itr_bonds_c->GetTargetAtom().GetAtomType() == GetAtomTypes().O_Tr2Tr2TrPi
              || itr_bonds_c->GetTargetAtom().GetAtomType() == GetAtomTypes().S_Tr2Tr2TrPi
            )
            {
              return true;
            }
          }
        }
        return false;
      }
    }

    //! @brief set bond type conjugations to conjugated/aromatic/or non-conjugated
    void AtomsCompleteStandardizer::SetConjugationOfBondTypes( AtomVector< AtomComplete> &ATOMS)
    {
      // set the bond types to be conjugated between pairs of conjugated atoms
      for
      (
        AtomVector< AtomComplete>::iterator itr_atom( ATOMS.Begin()), itr_atom_end( ATOMS.End());
        itr_atom != itr_atom_end;
        ++itr_atom
      )
      {
        // make a reference to the current atom
        AtomComplete &atom_a( *itr_atom);

        // skip non-conjugated atoms
        if( !atom_a.GetAtomType()->IsConjugated())
        {
          continue;
        }

        // make the bond type between all conjugated atoms conjugated properly
        const storage::Vector< BondConformational> &bonds_a( atom_a.GetBonds());
        for
        (
          storage::Vector< BondConformational>::const_iterator
            itr_bonds( bonds_a.Begin()), itr_bonds_end( bonds_a.End());
          itr_bonds != itr_bonds_end;
          ++itr_bonds
        )
        {
          // skip this bond if the target atom of the bond is non-conjugated
          if( !itr_bonds->GetTargetAtom().GetAtomType()->IsConjugated())
          {
            continue;
          }

          // get the current bond type
          const ConfigurationalBondType &bond_type( itr_bonds->GetBondType());

          // skip bonds that have bond order > 1; they are already conjugated or aromatic
          if( bond_type->GetNumberOfElectrons() > size_t( 2))
          {
            continue;
          }

          // skip bonds that are already labeled aromatic or amide
          if
          (
            bond_type->GetConjugation() == ConstitutionalBondTypeData::e_Aromatic
            || bond_type->GetConjugation() == ConstitutionalBondTypeData::e_Amide
          )
          {
            continue;
          }

          // set the bond type (iterating over all bonds, so a bidirectional set is unnecessary)
          atom_a.SetBondTypeMonoDirectional
          (
            itr_bonds->GetTargetAtom(),
            bond_type->WithConjugation
            (
              itr_bonds->GetBondType()->IsBondInRing() || !IsNonRotatableAmideLikeBond( atom_a, *itr_bonds)
              ? ConstitutionalBondTypeData::e_Conjugated
              : ConstitutionalBondTypeData::e_Amide
            )
          );
        }
      }
    }

    //! @brief set bond type conjugations to conjugated/aromatic/or non-conjugated
    void AtomsCompleteStandardizer::SetConjugationOfBondTypes()
    {
      AtomsCompleteStandardizer::SetConjugationOfBondTypes( m_Atoms);
    }

    //! @brief set the selected atom types for all atoms in the molecule
    //! @return true if all types were defined
    bool AtomsCompleteStandardizer::SetAtomTypes()
    {
      bool all_were_defined( true);
      storage::Vector< PossibleAtomTypesForAtom>::iterator itr_atom_types( m_PossibleAtomTypes.Begin());
      storage::Vector< size_t>::const_iterator itr_smallest_rings( m_SmallestRingSize.Begin());
      for
      (
        AtomVector< AtomComplete>::iterator itr_atom( m_Atoms.Begin()), itr_atom_end( m_Atoms.End());
        itr_atom != itr_atom_end;
        ++itr_atom, ++itr_atom_types, ++itr_smallest_rings
      )
      {
        if( itr_atom_types->GetMostStableType().IsDefined())
        {
          // determine final atom type, taking into account neighboring atom types
          itr_atom_types->Finalize( *itr_atom, *itr_smallest_rings);
          itr_atom->SetAtomType( itr_atom_types->GetMostStableType());
        }
        else
        {
          // atom type not found; remove any existing gasteiger atom type
          all_were_defined = false;
        }
      }
      return all_were_defined;
    }

    //! @brief get all atom types matching a given atom considering its bonds
    //! @param ATOM atom of interest
    //! @return PossibleAtomTypesForAtom
    PossibleAtomTypesForAtom AtomsCompleteStandardizer::GetPossibleTypesForAtom( AtomComplete &ATOM)
    {
      // get a reference to the bonds
      const storage::Vector< BondConformational> &bonds( ATOM.GetBonds());

      // determine number of sigma orbitals
      int nr_bonds( bonds.GetSize());

      // count conjugated bonds, that is, that had a partial bond order
      int number_aromatic_bonds( 0);
      int number_bonds_unknown_order( 0);

      // count electrons in bonds
      int nr_e_in_bonds( 0);
      for
      (
        storage::Vector< BondConformational>::const_iterator itr_bonds( bonds.Begin()), itr_bonds_end( bonds.End());
        itr_bonds != itr_bonds_end;
        ++itr_bonds
      )
      {
        nr_e_in_bonds += itr_bonds->GetBondType()->GetNumberOfElectrons();
        number_bonds_unknown_order += !itr_bonds->GetBondType()->IsBondOrderKnown();
        number_aromatic_bonds += itr_bonds->GetBondType()->GetBondData( ConfigurationalBondTypeData::e_IsAromatic);
      }

      if( number_bonds_unknown_order == size_t( 1))
      {
        // discontinuous aromatic ring, nr_e_in_bonds by 1 to compensate
        ++nr_e_in_bonds;
      }

      // divide by 2 to count electrons in orbitals
      nr_e_in_bonds /= 2;

      // handle implicit hydrogens and charge
      const ElectronConfiguration &e_config( ATOM.GetElementType()->GetElectronConfiguration());
      const int valence_e( e_config.ValenceElectronsSP() - ATOM.GetCharge());

      // Find the total number of unpaired valence electrons of *itr_atom
      const int unpaired_valence_e
      (
        nr_e_in_bonds > 4
        ? std::max( valence_e, int( e_config.MaxValenceElectronsSP()) - valence_e)
        : std::min( valence_e, int( e_config.MaxValenceElectronsSP()) - valence_e)
      );
      const int implicit_hydrogens( std::min( int( 4 - nr_bonds), unpaired_valence_e - nr_e_in_bonds));

      // add implicit hydrogens, unless the number of hydrogens is unclear because
      // hydrogens are absent and there were bonds with unknown order
      if( implicit_hydrogens > 0)
      {
        nr_e_in_bonds += implicit_hydrogens;
        nr_bonds += implicit_hydrogens;
      }

      PossibleAtomTypesForAtom possible_types
      (
        ATOM.GetElementType(),
        nr_e_in_bonds,
        nr_bonds,
        ATOM.GetCharge(),
        number_aromatic_bonds
      );

      // check that the possible type was defined
      if( !possible_types.GetMostStableType().IsDefined())
      {
        // handle special cases of common sdf canonicalization errors induced by some software *cough* openbabel *cough*
        if( ATOM.GetElementType() == GetElementTypes().e_Nitrogen && nr_bonds == size_t( 2) && nr_e_in_bonds == size_t( 5))
        {
          // demote the triple bond to double bond
          const size_t trip_bond( bonds( 1).GetBondType()->GetBondData( ConfigurationalBondTypeData::e_BondOrder) == size_t( 3) ? 1 : 0);
          const size_t atom_id( m_Atoms.GetAtomIndex( bonds( trip_bond).GetTargetAtom()));
          ATOM.SetBondTypeTo( m_Atoms( atom_id), GetConfigurationalBondTypes().e_ConjugatedDoubleBond_X);
          m_NeedToRecomputeOtherAtoms = true;
          return GetPossibleTypesForAtom( ATOM);
        }

        if( ATOM.GetElementType() == GetElementTypes().e_Sulfur && nr_bonds == size_t( 4) && nr_e_in_bonds == nr_bonds && implicit_hydrogens <= 0)
        {
          size_t oxygen_pos( util::GetUndefinedSize_t());
          // find an oxygen that this sulfur is bonded to that is itself terminal or has an implicit hydrogen
          for( size_t bond_nr( 0), n_bonds( bonds.GetSize()); bond_nr < n_bonds; ++bond_nr)
          {
            if( bonds( bond_nr).GetTargetAtom().GetElementType() == GetElementTypes().e_Oxygen)
            {
              if( bonds( bond_nr).GetTargetAtom().GetBonds().GetSize() == size_t( 1))
              {
                oxygen_pos = bond_nr;
                break;
              }
            }
          }
          if( util::IsDefined( oxygen_pos))
          {
            const size_t atom_id( m_Atoms.GetAtomIndex( bonds( oxygen_pos).GetTargetAtom()));
            ATOM.SetBondTypeTo( m_Atoms( atom_id), GetConfigurationalBondTypes().e_ConjugatedDoubleBond);
            m_NeedToRecomputeOtherAtoms = true;
            return GetPossibleTypesForAtom( ATOM);
          }
        }
        if( ATOM.GetElementType() == GetElementTypes().e_Carbon && nr_bonds == size_t( 4) && nr_e_in_bonds == size_t( 5) && implicit_hydrogens <= 0)
        {
          size_t double_pos( util::GetUndefinedSize_t());
          // find an oxygen that this sulfur is bonded to that is itself terminal or has an implicit hydrogen
          for( size_t bond_nr( 0), n_bonds( bonds.GetSize()); bond_nr < n_bonds; ++bond_nr)
          {
            if( bonds( bond_nr).GetBondType()->GetNumberOfElectrons() == size_t( 4))
            {
              double_pos = bond_nr;
            }
          }
          if( util::IsDefined( double_pos))
          {
            const size_t atom_id( m_Atoms.GetAtomIndex( bonds( double_pos).GetTargetAtom()));
            ATOM.SetBondTypeTo( m_Atoms( atom_id), bonds( double_pos).GetBondType()->WithOrder( 1));
            m_NeedToRecomputeOtherAtoms = true;
            return GetPossibleTypesForAtom( ATOM);
          }
        }

        storage::Map< std::string, size_t> bond_type_counts;

        for
        (
            storage::Vector< BondConformational>::const_iterator itr( bonds.Begin()), itr_end( bonds.End());
            itr != itr_end;
            ++itr
        )
        {
          ++bond_type_counts[ itr->GetBondType().GetName()];
        }

        std::string connections;
        for
        (
            storage::Map< std::string, size_t>::const_iterator
            itr( bond_type_counts.Begin()),
            itr_end( bond_type_counts.End());
            itr != itr_end;
            ++itr
        )
        {
          connections += util::Format()( itr->second) + " " + itr->first + "s,";
        }

        BCL_MessageVrb
        (
          " Found no types for element "
          + ATOM.GetElementType()->GetChemicalSymbol() + " with "
          + std::string( " # electrons in bonds: ") + util::Format()( nr_e_in_bonds) + ","
          + std::string( " # bonds: ") + util::Format()( nr_bonds) + ","
          + std::string( " Charge: ") + util::Format()( ATOM.GetCharge()) + ","
          + std::string( " # implicit H: ") + util::Format()( std::max( implicit_hydrogens, 0)) + ","
          + std::string( " Bonds: ") + connections
        );
      }
      else
      {
        if( possible_types.GetMostStableType()->GetFormalCharge() != ATOM.GetCharge())
        {
          // update the charge
          ATOM.SetCharge( possible_types.GetMostStableType()->GetFormalCharge());
        }
      }

      return possible_types;
    } // GetPossibleTypesForAtom

    //! @brief using known atom types, determine a possible bond type for bonds that have unknown order
    void AtomsCompleteStandardizer::DetermineUnknownBondOrders()
    {
      // make a set of atom ids with any bond orders that are unknown
      // we use ids rather than ShPtr's to atoms with configuration, id's are ordered by position in the sdf file, which
      // allows the algorithm below to make a unique choice of bond order for a given small molecule

      // track the # of bonds with unknown order for each atom
      storage::Vector< size_t> bonds_with_unknown_order( m_Atoms.GetSize(), size_t( 0));

      // track the # of electrons in bonds with unknown order for each atom
      storage::Vector< size_t> e_with_unknown_order( m_Atoms.GetSize(), size_t( 0));

      // make an adjacency list containing indices of atoms connected to atom i by a bond with unknown order
      storage::Vector< storage::List< size_t> > adjacency_list( m_Atoms.GetSize());

      // make a list of all indices of atoms that have unknown bond orders
      storage::List< size_t> unknown_bond_order_atom_ids;

      size_t atom_index( 0);

      for
      (
        AtomVector< AtomComplete>::iterator itr_atom( m_Atoms.Begin()), itr_atom_end( m_Atoms.End());
        itr_atom != itr_atom_end;
        ++itr_atom, ++atom_index
      )
      {
        const AtomType atom_type( itr_atom->GetAtomType());

        size_t number_electrons_in_bonds_with_known_order( 0);
        size_t number_bonds_with_known_order( 0);

        // get the atoms this atom is bonded to
        const storage::Vector< BondConformational> &connectivity( itr_atom->GetBonds());

        // walk over them, adding up electrons in bonds with known order
        for
        (
          storage::Vector< BondConformational>::const_iterator
            itr_bond( connectivity.Begin()),
            itr_bond_end( connectivity.End());
          itr_bond != itr_bond_end;
          ++itr_bond
        )
        {
          // add up electrons in bonds with known order
          if( itr_bond->GetBondType()->IsBondOrderKnown())
          {
            number_electrons_in_bonds_with_known_order += itr_bond->GetBondType()->GetNumberOfElectrons() / 2;
            ++number_bonds_with_known_order;
          }
          else
          {
            adjacency_list( atom_index).PushBack( m_Atoms.GetAtomIndex( itr_bond->GetTargetAtom()));
          }
        }

        // skip atoms which have all bonds with a defined order
        if( number_bonds_with_known_order == connectivity.GetSize())
        {
          continue;
        }

        // update list of atoms with bonds of unknown order
        unknown_bond_order_atom_ids.PushBack( atom_index);

        // compute # bonds with unknown order
        const size_t nr_bonds_with_unknown_order( connectivity.GetSize() - number_bonds_with_known_order);
        bonds_with_unknown_order( atom_index) = nr_bonds_with_unknown_order;

        // compute number of implicit H bonds
        const size_t number_implicit_h_bonds( atom_type->GetNumberBonds() - connectivity.GetSize());

        // update number of electrons in bonds with known order; since implicit h bonds have known order
        number_electrons_in_bonds_with_known_order += number_implicit_h_bonds;

        e_with_unknown_order( atom_index) =
          atom_type->GetNumberElectronsInBonds() - number_electrons_in_bonds_with_known_order;

        BCL_MessageVrb
        (
          util::Format()( atom_index) + " atom type " + atom_type.GetName()
          + " bonds unknown: " + util::Format()( bonds_with_unknown_order( atom_index))
          + " e unknown: " + util::Format()( e_with_unknown_order( atom_index))
        );
      }

      // keep track of the last arbitrary choice to enable backtracking
      // backtracking is used successfully in only one case in the whole CSD (55k molecules); in that case the
      // ring system is identified incorrectly as aromatic
      storage::List< sdf::BondInfo> last_bonds_changed;
      size_t last_arbitrary_bond_order_choice( 0);

      while( !unknown_bond_order_atom_ids.IsEmpty())
      {
        // handle all cases where the bond order is determined locally (no more than two bonds away from this bond)
        bool change_was_made( true);
        bool change_success( true);
        while( change_was_made)
        {
          change_was_made = false;
          for
          (
            storage::List< size_t>::iterator
              itr_atom_id( unknown_bond_order_atom_ids.Begin()), itr_atom_id_end( unknown_bond_order_atom_ids.End());
            itr_atom_id != itr_atom_id_end;
            // iteration in loop
          )
          {
            const size_t atom_index( *itr_atom_id);
            if( !adjacency_list( atom_index).GetSize()) // skip atoms that have all bonds with well-defined types
            {
              storage::List< size_t>::iterator old_itr( itr_atom_id++);
              unknown_bond_order_atom_ids.RemoveElement( old_itr);
              continue;
            }

            // get the atoms this atom is bonded to
            storage::List< size_t> &connectivity( adjacency_list( atom_index));

            size_t number_electrons_in_unknown_order_bonds_neighbors( 0);
            size_t number_unknown_order_bonds_neighbors( 0);

            // walk over them, adding up electrons in bonds with known order
            for
            (
              storage::List< size_t>::const_iterator itr_bond( connectivity.Begin()), itr_bond_end( connectivity.End());
              itr_bond != itr_bond_end;
              ++itr_bond
            )
            {
              // add up electrons in bonds with known order
              number_electrons_in_unknown_order_bonds_neighbors += e_with_unknown_order( *itr_bond);
              number_unknown_order_bonds_neighbors += bonds_with_unknown_order( *itr_bond);
            }

            // check whether the unknown bonds have to be single bonds or whether there is only
            // one unknown bond, in which case the bond order is just the difference between the known orders and the
            // atom types order
            if
            (
              number_electrons_in_unknown_order_bonds_neighbors ==
                e_with_unknown_order( atom_index) + number_unknown_order_bonds_neighbors - bonds_with_unknown_order( atom_index)
            )
            {
              // walk over them, adding up electrons in bonds with known order
              for
              (
                storage::List< size_t>::const_iterator itr_bond( connectivity.Begin()), itr_bond_end( connectivity.End());
                itr_bond != itr_bond_end;
                ++itr_bond
              )
              {
                const ConfigurationalBondType old_bond_type( m_Atoms( atom_index).GetBondTypeTo( m_Atoms( *itr_bond)));
                last_bonds_changed.PushBack( sdf::BondInfo( atom_index, *itr_bond, old_bond_type));

                const ConfigurationalBondType new_bond_type
                (
                  old_bond_type->WithOrder( e_with_unknown_order( *itr_bond) + 1 - bonds_with_unknown_order( *itr_bond))
                );
                m_Atoms( atom_index).SetBondTypeTo( m_Atoms( *itr_bond), new_bond_type);
                e_with_unknown_order( atom_index) -= new_bond_type->GetBondData( ConfigurationalBondTypeData::e_BondOrder);
                --bonds_with_unknown_order( atom_index);
                e_with_unknown_order( *itr_bond) -= new_bond_type->GetBondData( ConfigurationalBondTypeData::e_BondOrder);
                --bonds_with_unknown_order( *itr_bond);
                storage::List< size_t>::iterator itr_bond_bonds( adjacency_list( *itr_bond).Begin());
                while( *itr_bond_bonds != atom_index)
                {
                  ++itr_bond_bonds;
                }
                adjacency_list( *itr_bond).RemoveElement( itr_bond_bonds);
              }
              connectivity.Reset();
              change_was_made = true;
              storage::List< size_t>::iterator old_itr( itr_atom_id++);
              unknown_bond_order_atom_ids.RemoveElement( old_itr);
              continue;
            }
            else if
            (
              number_electrons_in_unknown_order_bonds_neighbors <
              e_with_unknown_order( atom_index) + number_unknown_order_bonds_neighbors - bonds_with_unknown_order( atom_index)
            )
            {
              change_success = false;
              change_was_made = false;
              last_arbitrary_bond_order_choice = 3;
              BCL_MessageVrb
              (
                "Insufficient electrons around atom: " + util::Format()( atom_index) + " to serve kekulization"
              );
              break;
            }

            ++itr_atom_id;
          }
        }

        if( change_success)
        {
          if( last_arbitrary_bond_order_choice > size_t( 1))
          {
            BCL_MessageVrb
            (
              "Backtracking used to solve kekule structure on molecule " + m_ID
            );
          }

          last_bonds_changed.Reset();
          while( !unknown_bond_order_atom_ids.IsEmpty())
          {
            const size_t id( unknown_bond_order_atom_ids.FirstElement());
            if( bonds_with_unknown_order( id) == size_t( 0))
            {
              unknown_bond_order_atom_ids.PopFront();
              continue;
            }

            const size_t target_id( adjacency_list( id).FirstElement());
            BCL_MessageVrb
            (
              "Arbitrarily assigning the first bond with undefined type between atom:" + util::Format()( id)
              + " e-: " + util::Format()( e_with_unknown_order( id))
              + " #bonds: " + util::Format()( bonds_with_unknown_order( id))
              + " and atom: " + util::Format()( adjacency_list( id).FirstElement())
              + " bond order 1"
            );

            last_arbitrary_bond_order_choice = 1;
            const ConfigurationalBondType old_bond_type( m_Atoms( id).GetBondTypeTo( m_Atoms( target_id)));
            last_bonds_changed.PushBack( sdf::BondInfo( id, target_id, old_bond_type));
            const ConfigurationalBondType new_bond_type( old_bond_type->WithOrder( size_t( 1)));
            m_Atoms( id).SetBondTypeTo( m_Atoms( target_id), new_bond_type);
            e_with_unknown_order( id) -= new_bond_type->GetBondData( ConfigurationalBondTypeData::e_BondOrder);
            --bonds_with_unknown_order( id);
            e_with_unknown_order( target_id) -= new_bond_type->GetBondData( ConfigurationalBondTypeData::e_BondOrder);
            --bonds_with_unknown_order( target_id);
            storage::List< size_t>::iterator itr_bond_bonds( adjacency_list( target_id).Begin());
            while( *itr_bond_bonds != id)
            {
              ++itr_bond_bonds;
            }
            adjacency_list( target_id).RemoveElement( itr_bond_bonds);
            adjacency_list( id).PopFront();
            break;
          }
          continue;
        }
        if( last_arbitrary_bond_order_choice != 0)
        {
          ++last_arbitrary_bond_order_choice;
          // undo all the changes that were done
          for
          (
            storage::List< sdf::BondInfo>::const_reverse_iterator
              itr( last_bonds_changed.ReverseBegin()), itr_end( last_bonds_changed.ReverseEnd());
            itr != itr_end;
            ++itr
          )
          {
            const size_t &atom_a_index( itr->GetAtomIndexLow());
            const size_t &atom_b_index( itr->GetAtomIndexHigh());
            const ConfigurationalBondType &old_bond_type( itr->GetConfigurationalBondType());
            const ConfigurationalBondType new_bond_type( m_Atoms( atom_a_index).GetBondTypeTo( m_Atoms( atom_b_index)));

            m_Atoms( atom_a_index).SetBondTypeTo( m_Atoms( atom_b_index), old_bond_type);
            e_with_unknown_order( atom_a_index) += new_bond_type->GetBondData( ConfigurationalBondTypeData::e_BondOrder);
            if( ++bonds_with_unknown_order( atom_a_index) == 1)
            {
              unknown_bond_order_atom_ids.PushFront( atom_a_index);
            }
            e_with_unknown_order( atom_b_index) += new_bond_type->GetBondData( ConfigurationalBondTypeData::e_BondOrder);
            if( ++bonds_with_unknown_order( atom_b_index) == 1)
            {
              unknown_bond_order_atom_ids.PushFront( atom_b_index);
            }
            adjacency_list( atom_a_index).PushFront( atom_b_index);
            adjacency_list( atom_b_index).PushFront( atom_a_index);
          }
          unknown_bond_order_atom_ids.Sort( std::less< size_t>());

          if( last_arbitrary_bond_order_choice < size_t( 3))
          {
            ++last_arbitrary_bond_order_choice;
            const size_t atom_a_index( last_bonds_changed.FirstElement().GetAtomIndexLow());
            const size_t atom_b_index( last_bonds_changed.FirstElement().GetAtomIndexHigh());
            const ConfigurationalBondType &old_bond_type( last_bonds_changed.FirstElement().GetConfigurationalBondType());
            const size_t max_bond_order
            (
              std::min
              (
                e_with_unknown_order( atom_a_index) + size_t( 1) - bonds_with_unknown_order( atom_a_index),
                e_with_unknown_order( atom_b_index) + size_t( 1) - bonds_with_unknown_order( atom_b_index)
              )
            );
            last_bonds_changed.Reset();
            if( last_arbitrary_bond_order_choice <= max_bond_order)
            {
              m_Atoms( atom_a_index).SetBondTypeTo
              (
                m_Atoms( atom_b_index),
                old_bond_type->WithOrder( last_arbitrary_bond_order_choice)
              );
              last_bonds_changed.PushBack( sdf::BondInfo( atom_a_index, atom_b_index, old_bond_type));
              BCL_MessageVrb
              (
                "Arbitrarily assigning the first bond with undefined type between atom:" + util::Format()( atom_a_index)
                + " e-: " + util::Format()( e_with_unknown_order( atom_a_index))
                + " #bonds: " + util::Format()( bonds_with_unknown_order( atom_a_index))
                + " and atom: " + util::Format()( atom_b_index)
                + " bond order " + util::Format()( last_arbitrary_bond_order_choice)
              );
              e_with_unknown_order( atom_a_index) -= last_arbitrary_bond_order_choice;
              --bonds_with_unknown_order( atom_a_index);
              e_with_unknown_order( atom_b_index) -= last_arbitrary_bond_order_choice;
              --bonds_with_unknown_order( atom_b_index);
              adjacency_list( atom_a_index).PopFront();
              adjacency_list( atom_b_index).PopFront();
              continue;
            }
          }
        }
        break;
      }
      if( unknown_bond_order_atom_ids.GetSize())
      {
        BCL_MessageCrt
        (
          "Could not kekulize molecule " + m_ID
          + " setting atom types of un-kekulizable atoms to simple (non-gasteiger) types.\n"
          + "If this is an error, please kekulize (e.g. replace aromatic bonds with single/double bonds)"
        );
        BCL_MessageVrb( "Index\tAtomType\t# Bonds Unknown\t# E Unknown");
        for
        (
          storage::List< size_t>::const_iterator
            i( unknown_bond_order_atom_ids.Begin()), i_end( unknown_bond_order_atom_ids.End());
          i != i_end;
          ++i
        )
        {
          BCL_MessageVrb
          (
            util::Format()( *i) + '\t'
            + m_Atoms( *i).GetAtomType().GetName() + '\t'
            + util::Format()( bonds_with_unknown_order( *i)) + '\t'
            + util::Format()( e_with_unknown_order( *i))
          );
          m_Atoms( *i).SetAtomType( AtomTypes::GetAtomType( m_Atoms( *i).GetAtomType()->GetBaseAtomType()));
        }
      }
    }

    //! @brief remove all bonds to / from group 1 elements if the # of bonds from this atom exceeds 1
    //! In configuration compounds, sometimes group 1 elements form a configuration center; the group 1
    //! atom becomes charged and centered around other atoms with a partial charge
    //! The bcl only supports normal covalent bonds at this point, so these bonds are thrown out
    AtomVector< AtomComplete> &AtomsCompleteStandardizer::RemoveObviousIonicBonds( AtomVector< AtomComplete> &ATOMS)
    {
      // set the bond types to be conjugated between pairs of conjugated atoms
      for
      (
        AtomVector< AtomComplete>::iterator itr_atom( ATOMS.Begin()), itr_atom_end( ATOMS.End());
        itr_atom != itr_atom_end;
        ++itr_atom
      )
      {
        // get the number of valence sp electrons for this atom ( == 1 for group 1 )
        const size_t group( itr_atom->GetElementType()->GetElectronConfiguration().ValenceElectronsSP());
        if( group != size_t( 1))
        {
          continue;
        }

        // only consider atoms with more than group bonds
        if( itr_atom->GetBonds().GetSize() <= size_t( 1))
        {
          continue;
        }

        // remove all the bonds bidirectionally
        const storage::Vector< BondConformational> ionic_bonds( itr_atom->GetBonds());

        for
        (
          storage::Vector< BondConformational>::const_iterator
            itr_bond( ionic_bonds.Begin()), itr_bond_end( ionic_bonds.End());
          itr_bond != itr_bond_end;
          ++itr_bond
        )
        {
          ATOMS( itr_bond->GetTargetAtom()).AddValenceBondByRemovingBondTo( *itr_atom, false);
        }

        itr_atom->SetBonds( storage::Vector< BondConformational>());
      }

      // drop conjugation of bonds that already have a defined order
      for
      (
        AtomVector< AtomComplete>::iterator itr_atom( ATOMS.Begin()), itr_atom_end( ATOMS.End());
        itr_atom != itr_atom_end;
        ++itr_atom
      )
      {
        for
        (
          storage::Vector< BondConformational>::const_iterator
            itr_bond( itr_atom->GetBonds().Begin()), itr_bond_end( itr_atom->GetBonds().End());
          itr_bond != itr_bond_end;
          ++itr_bond
        )
        {
          if
          (
            itr_bond->GetBondType()->GetBondData( ConfigurationalBondTypeData::e_BondOrder) == size_t( 1)
            &&
            (
              itr_bond->GetBondType()->GetConjugation() == ConstitutionalBondTypeData::e_Conjugated
              || itr_bond->GetBondType()->GetConjugation() == ConstitutionalBondTypeData::e_Amide
            )
          )
          {
            itr_atom->SetBondTypeMonoDirectional
            (
              itr_bond->GetTargetAtom(),
              itr_bond->GetBondType()->WithConjugation( ConstitutionalBondTypeData::e_Nonconjugated)
            );
          }
          else if
          (
            itr_bond->GetBondType()->IsBondOrderKnown()
            && itr_bond->GetBondType()->GetConjugation() == ConstitutionalBondTypeData::e_Aromatic
          )
          {
            itr_atom->SetBondTypeMonoDirectional
            (
              itr_bond->GetTargetAtom(),
              itr_bond->GetBondType()->GetBondData( ConfigurationalBondTypeData::e_BondOrder) == size_t( 1)
              ? itr_bond->GetBondType()->WithConjugation( ConstitutionalBondTypeData::e_Nonconjugated)
              : itr_bond->GetBondType()->WithConjugation( ConstitutionalBondTypeData::e_Conjugated)
            );
          }
        }
      }
      return ATOMS;
    }

    //! @brief try to neutralize a given molecule
    void AtomsCompleteStandardizer::TryNeutralize( AtomVector< AtomComplete> &ATOMS, const sdf::NeutralizationPref &NEUTRALIZATION_PREF)
    {
      if( NEUTRALIZATION_PREF == sdf::e_None)
      {
        return;
      }

      // try removing small counter-ions
      storage::Vector< size_t> is_not_counterion;
      is_not_counterion.AllocateMemory( ATOMS.GetSize());
      for( size_t i( 0), na( ATOMS.GetSize()); i < na; ++i)
      {
        if( !ATOMS( i).GetAtomType().IsDefined() || !ATOMS( i).GetAtomType()->IsGasteigerAtomType())
        {
          return;
        }
        if
        (
          ( ATOMS( i).GetElementType()->GetMainGroup() != size_t( 1) || ATOMS( i).GetElementType() == GetElementTypes().e_Hydrogen)
          &&
          ( !ATOMS( i).GetAtomType().IsDefined() || ATOMS( i).GetAtomType()->GetNumberBonds() != size_t( 0))
        )
        {
          is_not_counterion.PushBack( i);
        }
      }
      if( is_not_counterion.GetSize() != ATOMS.GetSize())
      {
        ATOMS.Reorder( is_not_counterion);
      }

      bool allow_arom_chng( sdf::GetNeutralizationPrefAllowsAromaticityChange( NEUTRALIZATION_PREF));
      bool allow_bnd_order_chng( sdf::GetNeutralizationPrefAllowsBondOrderChange( NEUTRALIZATION_PREF));

      // make vectors for atoms on which to add electrons or remove them
      size_t n_atoms( ATOMS.GetSize()), number_could_neutralize( 0);
      storage::Vector< size_t> atoms_could_change( n_atoms, size_t( 0));
      bool had_h( false);
      bool had_valences( false);
      for( size_t i( 0); i < n_atoms; ++i)
      {
        const AtomType atom_type( ATOMS( i).GetAtomType());
        if( atom_type == GetAtomTypes().H_S)
        {
          had_h = true;
        }
        if( ATOMS( i).GetNumberValenceBonds())
        {
          had_valences = true;
        }
        // skip uncharged atoms
        if( !atom_type->GetFormalCharge())
        {
          continue;
        }
        // skip aromatic atoms unless the neutralization preference allows for changing aromaticity
        if
        (
          !allow_arom_chng
          && ATOMS( i).CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsAromatic, 1)
        )
        {
          continue;
        }
        atoms_could_change( i) = 1;
        ++number_could_neutralize;
      }

      // escape if no formal charges to neutralize
      if( !number_could_neutralize)
      {
        return;
      }

      const bool no_add_h( had_valences && !had_h);
      bool did_something( false);

      if( allow_bnd_order_chng)
      {
        // try to add electrons to bonds first (e.g. increase bond order)
        for( size_t i( 0); i < n_atoms; ++i)
        {
          if( !atoms_could_change( i))
          {
            continue;
          }

          PossibleAtomTypesForAtom new_type_i
          (
            ATOMS( i).GetElementType(),
            ATOMS( i).GetAtomType()->GetNumberElectronsInBonds() + 1,
            ATOMS( i).GetAtomType()->GetNumberBonds(),
            0,
            false
          );

          // test whether bond order could be increased for this type
          if( !new_type_i.GetNumberPossibleTypes() || new_type_i.GetMostStableType()->GetFormalCharge())
          {
            continue;
          }

          // check whether any pairs of atoms in add or change any bonded e are connected by anything other than a triple bond
          // if so, we can neutralize charge by increasing the bond order of the bond between them
          const storage::Vector< BondConformational> &bonds( ATOMS( i).GetBonds());
          size_t most_electronegative_candidate( util::GetUndefinedSize_t());

          for( size_t bonded_atom_id( 0), n_bonds( bonds.GetSize()); bonded_atom_id < n_bonds; ++bonded_atom_id)
          {
            const size_t bonded_atom_index( ATOMS.GetAtomIndex( bonds( bonded_atom_id).GetTargetAtom()));

            // skip triple bonds
            const ConfigurationalBondType bond_type( bonds( bonded_atom_id).GetBondType());
            if( bond_type->GetBondData( ConfigurationalBondTypeData::e_BondOrder) == size_t( 3))
            {
              continue;
            }

            PossibleAtomTypesForAtom new_type_bonded
            (
              ATOMS( bonded_atom_index).GetElementType(),
              ATOMS( bonded_atom_index).GetAtomType()->GetNumberElectronsInBonds() + 1,
              ATOMS( bonded_atom_index).GetAtomType()->GetNumberBonds(),
              0,
              false
            );
            if( new_type_bonded.GetNumberPossibleTypes() && !new_type_bonded.GetMostStableType()->GetFormalCharge())
            {
              if( !util::IsDefined( most_electronegative_candidate) || ATOMS( bonded_atom_index).GetElementType()->GetProperty( ElementTypeData::e_ElectroNegativity) > ATOMS( most_electronegative_candidate).GetElementType()->GetProperty( ElementTypeData::e_ElectroNegativity))
              {
                most_electronegative_candidate = bonded_atom_index;
              }
            }
          }
          if( util::IsDefined( most_electronegative_candidate))
          {
            PossibleAtomTypesForAtom new_type_bonded
            (
              ATOMS( most_electronegative_candidate).GetElementType(),
              ATOMS( most_electronegative_candidate).GetAtomType()->GetNumberElectronsInBonds() + 1,
              ATOMS( most_electronegative_candidate).GetAtomType()->GetNumberBonds(),
              0,
              false
            );
            // only perform the neutralization if the user requested it. If the user requested only pH-dependent
            // neutralization, we want to ignore atoms for which a bond order change would neutralize the species and
            // decrease the net charge of the molecule

            atoms_could_change( i) = 0;
            atoms_could_change( most_electronegative_candidate) = atoms_could_change( most_electronegative_candidate) = 0;
            BCL_MessageVrb
            (
              "NEUTRAL: upped bond order ... new types are: " + new_type_i.GetMostStableType().GetName()
              + " " + new_type_bonded.GetMostStableType().GetName()
            );
            const ConfigurationalBondType bond_type( ATOMS( i).GetBondTypeTo( ATOMS( most_electronegative_candidate)));
            ATOMS( i).SetBondTypeTo
                (
                  ATOMS( most_electronegative_candidate),
                  bond_type->WithOrder( bond_type->GetBondData( ConfigurationalBondTypeData::e_BondOrder) + 1)
                );
            did_something = true;
            ATOMS( i).SetAtomType( new_type_i.GetMostStableType());
            ATOMS( most_electronegative_candidate).SetAtomType( new_type_bonded.GetMostStableType());
          }
        }

        // check whether any pairs of atoms in remove or change any bonded e are connected by a non-single bond, if so, we
        // can neutralize charge by decreasing bond order of the bond between them.
        for( size_t i( 0); i < n_atoms; ++i)
        {
          if( !atoms_could_change( i))
          {
            continue;
          }

          const AtomType atom_type_i( ATOMS( i).GetAtomType());

          // skip atoms with only single-bonds
          if( atom_type_i->GetNumberElectronsInBonds() == atom_type_i->GetNumberBonds())
          {
            continue;
          }
          PossibleAtomTypesForAtom new_type_i
          (
            ATOMS( i).GetElementType(),
            atom_type_i->GetNumberElectronsInBonds() - 1,
            atom_type_i->GetNumberBonds(),
            0,
            false
          );

          // test whether bond order could be decreased
          if( !new_type_i.GetNumberPossibleTypes() || new_type_i.GetMostStableType()->GetFormalCharge())
          {
            continue;
          }

          // check whether any pairs of atoms in remove any bonded e are connected by anything other than a single bond
          // if so, we can neutralize charge by decreasing the bond order of the bond between them
          const storage::Vector< BondConformational> &bonds( ATOMS( i).GetBonds());
          for( size_t bonded_atom_id( 0), n_bonds( bonds.GetSize()); bonded_atom_id < n_bonds; ++bonded_atom_id)
          {
            const size_t bonded_atom_index( ATOMS.GetAtomIndex( bonds( bonded_atom_id).GetTargetAtom()));
            if( !atoms_could_change( bonded_atom_index))
            {
              continue; // skip atoms that don't need fewer e- for neutralization
            }

            // skip single bonds
            const ConfigurationalBondType bond_type( bonds( bonded_atom_id).GetBondType());
            if( bond_type->GetBondData( ConfigurationalBondTypeData::e_BondOrder) == size_t( 1))
            {
              continue;
            }

            PossibleAtomTypesForAtom new_type_bonded
            (
              ATOMS( bonded_atom_index).GetElementType(),
              ATOMS( bonded_atom_index).GetAtomType()->GetNumberElectronsInBonds() - 1,
              ATOMS( bonded_atom_index).GetAtomType()->GetNumberBonds(),
              0,
              false
            );
            if( new_type_bonded.GetNumberPossibleTypes() && !new_type_bonded.GetMostStableType()->GetFormalCharge())
            {
              // only perform the neutralization if the user requested it. If the user requested only pH-dependent
              // neutralization, we want to ignore atoms for which a bond order change would neutralize the species and
              // decrease the net charge of the molecule
              atoms_could_change( i) = 0;
              atoms_could_change( bonded_atom_index) = 0;

              BCL_MessageVrb
              (
                "Neutralization lowered bond order ... new types are: " + new_type_i.GetMostStableType().GetName() + " "
                + new_type_bonded.GetMostStableType().GetName()
              );
              did_something = true;
              ATOMS( i).SetBondTypeTo
              (
                ATOMS( bonded_atom_index),
                bond_type->WithOrder( bond_type->GetBondData( ConfigurationalBondTypeData::e_BondOrder) - 1)
              );

              ATOMS( i).SetAtomType( new_type_i.GetMostStableType());
              ATOMS( bonded_atom_index).SetAtomType( new_type_bonded.GetMostStableType());
            }
            break;
          }
        }
      }

      if( NEUTRALIZATION_PREF == sdf::e_BondOrder)
      {
        if( did_something)
        {
          AtomsCompleteStandardizer( ATOMS, "", true);
          if( !no_add_h)
          {
            HydrogensHandler::Saturate( ATOMS);
          }
        }
        return;
      }

      // pH-dependent removal/addition

      // try to remove hydrogens. Keep track of removed H
      storage::Vector< size_t> hs_to_remove( ATOMS.GetSize(), size_t( 0));

      bool need_to_remove_h_atoms( false);
      for( size_t i( 0); i < n_atoms; ++i)
      {
        if
        (
          !atoms_could_change( i)
          ||
          (
            ATOMS( i).GetNumberCovalentlyBoundHydrogens() == size_t( 0)
            && ATOMS( i).GetNumberofValenceBondsWithOrder( 1) == size_t( 0)
          )
        )
        {
          continue;
        }

        PossibleAtomTypesForAtom new_type_i
        (
          ATOMS( i).GetElementType(),
          ATOMS( i).GetAtomType()->GetNumberElectronsInBonds() - 1,
          ATOMS( i).GetAtomType()->GetNumberBonds() - 1,
          0,
          false
        );
        if( new_type_i.GetNumberPossibleTypes() && !new_type_i.GetMostStableType()->GetFormalCharge())
        {
          BCL_MessageVrb( "NEUTRAL: removed an H to form: " + new_type_i.GetMostStableType().GetName());
          ATOMS( i).SetAtomType( new_type_i.GetMostStableType());
          atoms_could_change( i) = 0;
          did_something = true;

          // find an H on this atom
          const storage::Vector< BondConformational> &bonds( ATOMS( i).GetBonds());
          for( size_t bonded_atom_id( 0), n_bonds( bonds.GetSize()); bonded_atom_id < n_bonds; ++bonded_atom_id)
          {
            if( bonds( bonded_atom_id).GetTargetAtom().GetElementType() == GetElementTypes().e_Hydrogen)
            {
              const size_t h_index( ATOMS.GetAtomIndex( bonds( bonded_atom_id).GetTargetAtom()));
              hs_to_remove( h_index) = 1; // this hydrogen will be removed
              need_to_remove_h_atoms = true;
              break;
            }
          }
        }
      }

      storage::Vector< size_t> atoms_add_h_to( ATOMS.GetSize(), size_t( 0));

      // check whether any atom in add bonded e could have another hydrogen, if so, do it
      for( size_t i( 0); i < n_atoms; ++i)
      {
        if( !atoms_could_change( i))
        {
          continue;
        }

        PossibleAtomTypesForAtom new_type_i
        (
          ATOMS( i).GetElementType(),
          ATOMS( i).GetAtomType()->GetNumberElectronsInBonds() + 1,
          ATOMS( i).GetAtomType()->GetNumberBonds() + 1,
          0,
          false
        );
        // Skip formation of N_TeTeTeTePi. This atom type is not legitimate by standard rules, and is only present to
        // allow reading a subset of CSD molecules where it appears.
        if( new_type_i.GetNumberPossibleTypes() && !new_type_i.GetMostStableType()->GetFormalCharge())
        {
          if
          (
            new_type_i.GetMostStableType() != GetAtomTypes().N_TeTeTeTePi
            && new_type_i.GetMostStableType() != GetAtomTypes().P_TeTeTeTePi
            && new_type_i.GetMostStableType() != GetAtomTypes().N_TrTrTrPiPi
          )
          {
            BCL_MessageVrb( "NEUTRAL: added H to form: " + new_type_i.GetMostStableType().GetName());
            ATOMS( i).SetAtomType( new_type_i.GetMostStableType());
            atoms_could_change( i) = 0;
            did_something = true;
            atoms_add_h_to( i) = 1;
            atoms_could_change.PushBack( 0);
            atoms_add_h_to.PushBack( 0);
            hs_to_remove.PushBack( 0);
          }
        }
      }

      // saturate partially if there were any atoms to add H's to
      if( !no_add_h && atoms_could_change.GetSize() > n_atoms)
      {
        HydrogensHandler::SaturatePartial( ATOMS, atoms_add_h_to);
        n_atoms = ATOMS.GetSize();
        atoms_could_change.Resize( n_atoms, 0);
        atoms_add_h_to.Resize( n_atoms, 0);
        hs_to_remove.Resize( n_atoms, 0);
      }

      // check whether it is necessary to remove explicit hydrogens
      if( need_to_remove_h_atoms)
      {
        storage::Vector< size_t> new_order;
        new_order.AllocateMemory( n_atoms);
        for( size_t i( 0); i < n_atoms; ++i)
        {
          if( !hs_to_remove( i))
          {
            new_order.PushBack( i);
          }
        }
        ATOMS.Reorder( new_order);
        atoms_could_change.Reorder( new_order);
        n_atoms = ATOMS.GetSize();
        need_to_remove_h_atoms = false;
      }

      // return if we don't need to do the inter-dependent pH/BondOrder-dependent stuff
      if( NEUTRALIZATION_PREF == sdf::e_pH)
      {
        if( did_something)
        {
          AtomsCompleteStandardizer( ATOMS, "", true);
          if( !no_add_h)
          {
            HydrogensHandler::Saturate( ATOMS);
          }
        }
        return;
      }
      hs_to_remove.SetAllElements( 0);
      hs_to_remove.Resize( n_atoms);

      // next, try to neutralize charges by finding isolated charged species that need more e-, removing H from an adjacent
      // atom (the least electronegative is best), then increasing bond order
      // try to add electrons to bonds first (e.g. increase bond order)
      for( size_t i( 0); i < n_atoms; ++i)
      {
        if( !atoms_could_change( i))
        {
          continue;
        }
        AtomType atom_type_i( ATOMS( i).GetAtomType());
        PossibleAtomTypesForAtom new_type_i
        (
          ATOMS( i).GetElementType(),
          atom_type_i->GetNumberElectronsInBonds() + 1,
          atom_type_i->GetNumberBonds(),
          0,
          false
        );
        // check whether bond order could be increased
        if( !new_type_i.GetNumberPossibleTypes() || new_type_i.GetMostStableType()->GetFormalCharge())
        {
          continue;
        }
        if
        (
          new_type_i.GetMostStableType() == GetAtomTypes().N_TeTeTeTePi
          || new_type_i.GetMostStableType() == GetAtomTypes().P_TeTeTeTePi
          || new_type_i.GetMostStableType() == GetAtomTypes().S_Te2TeTeTePi
          || new_type_i.GetMostStableType() == GetAtomTypes().N_TrTrTrPiPi
          || new_type_i.GetMostStableType() == GetAtomTypes().P_TrTrTrPiPi
        )
        {
          continue;
        }

        // check whether any pairs of atoms in add or change any bonded e are connected by anything other than a triple bond
        // if so, we can neutralize charge by increasing the bond order of the bond between them
        const storage::Vector< BondConformational> &bonds( ATOMS( i).GetBonds());

        size_t best_match( util::GetUndefined< size_t>());
        double best_electronegativity( 0.0);
        ConfigurationalBondType best_bond_type;
        AtomType best_atom_type;
        for( size_t bonded_atom_id( 0), n_bonds( bonds.GetSize()); bonded_atom_id < n_bonds; ++bonded_atom_id)
        {
          const size_t bonded_atom_index( ATOMS.GetAtomIndex( bonds( bonded_atom_id).GetTargetAtom()));

          if( bonded_atom_id >= n_atoms)
          {
            continue;
          }
          // skip other charged atoms
          if( atoms_could_change( bonded_atom_index))
          {
            continue;
          }

          // skip aromatic bonded atoms unless the preference allows changing it
          if
          (
            !allow_arom_chng
            && ATOMS( bonded_atom_index).CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsAromatic, 1)
          )
          {
            continue;
          }

          // skip triple bonds
          const ConfigurationalBondType bond_type( bonds( bonded_atom_id).GetBondType());
          if( bond_type->GetBondData( ConfigurationalBondTypeData::e_BondOrder) == size_t( 3))
          {
            continue;
          }

          // test whether H-removal followed by increase of bond order would yield a valid type
          PossibleAtomTypesForAtom new_type_bonded
          (
            ATOMS( bonded_atom_index).GetElementType(),
            ATOMS( bonded_atom_index).GetAtomType()->GetNumberElectronsInBonds(),
            ATOMS( bonded_atom_index).GetAtomType()->GetNumberBonds() - 1,
            0,
            false
          );

          // no valid type if H is removed
          if( !new_type_bonded.GetNumberPossibleTypes() || new_type_bonded.GetMostStableType()->GetFormalCharge())
          {
            continue;
          }

          if
          (
            new_type_bonded.GetMostStableType() == GetAtomTypes().N_TeTeTeTePi
            || new_type_bonded.GetMostStableType() == GetAtomTypes().P_TeTeTeTePi
            || new_type_bonded.GetMostStableType() == GetAtomTypes().S_Te2TeTeTePi
            || new_type_bonded.GetMostStableType() == GetAtomTypes().N_TrTrTrPiPi
            || new_type_bonded.GetMostStableType() == GetAtomTypes().P_TrTrTrPiPi
            || new_type_bonded.GetMostStableType() == GetAtomTypes().O_Di2DiPi2Pi
          )
          {
            continue;
          }

          // test whether the atom has any convalently-bound hydrogens
          if
          (
            !ATOMS( bonded_atom_index).GetNumberCovalentlyBoundHydrogens()
            && !ATOMS( bonded_atom_index).GetNumberofValenceBondsWithOrder( 1)
          )
          {
            continue;
          }

          const double eneg
          (
            bonds( bonded_atom_id).GetTargetAtom().GetElementType()->GetProperty( ElementTypeData::e_ElectroNegativity)
          );
          if( eneg > best_electronegativity)
          {
            best_electronegativity = eneg;
            best_match = bonded_atom_index;
            best_bond_type = bonds( bonded_atom_id).GetBondType();
            best_atom_type = new_type_bonded.GetMostStableType();
          }
        }
        if( util::IsDefined( best_match))
        {
          did_something = true;
          atoms_could_change( i) = 0;
          atoms_could_change( best_match) = 0;
          BCL_MessageVrb
          (
            "NEUTRAL: increased bond order and removed an H ... new types are: " + new_type_i.GetMostStableType().GetName()
            + " " + best_atom_type.GetName() + " old types were: " + atom_type_i.GetName() + " " + ATOMS( best_match).GetAtomType().GetName()
          );
          ATOMS( i).SetBondTypeTo
          (
            ATOMS( best_match),
            best_bond_type->WithOrder( best_bond_type->GetBondData( ConfigurationalBondTypeData::e_BondOrder) + 1)
          );
          ATOMS( i).SetAtomType( new_type_i.GetMostStableType());
          ATOMS( best_match).SetAtomType( best_atom_type);
          const storage::Vector< BondConformational> &bonds_b( ATOMS( best_match).GetBonds());
          for
          (
            size_t bonded_atom_b_id( 0), n_b_bonds( bonds_b.GetSize());
            bonded_atom_b_id < n_b_bonds;
            ++bonded_atom_b_id
          )
          {
            if( bonds_b( bonded_atom_b_id).GetTargetAtom().GetAtomType() == GetAtomTypes().H_S)
            {
              const size_t bonded_atom_h_index( ATOMS.GetAtomIndex( bonds_b( bonded_atom_b_id).GetTargetAtom()));
              hs_to_remove( bonded_atom_h_index) = 1; // this hydrogen will be removed
              need_to_remove_h_atoms = true;
              break;
            }
          }
        }
      }

      if( need_to_remove_h_atoms)
      {
        storage::Vector< size_t> new_order;
        new_order.AllocateMemory( n_atoms);
        for( size_t i( 0); i < n_atoms; ++i)
        {
          if( !hs_to_remove( i))
          {
            new_order.PushBack( i);
          }
        }
        ATOMS.Reorder( new_order);
        atoms_could_change.Reorder( new_order);
        n_atoms = ATOMS.GetSize();
        need_to_remove_h_atoms = false;
      }

      // finally, try to neutralize charges by finding isolated charged species and adding H to an adjacent atom
      // the more electronegative, the better, then decreasing bond order, if that would result in valid atom types.
      atoms_add_h_to.Resize( n_atoms);
      atoms_add_h_to.SetAllElements( 0);
      // next, try to neutralize charges by finding isolated charged species that need more e-, removing H from an adjacent
      // atom (the least electronegative is best), then increasing bond order
      // try to add electrons to bonds first (e.g. increase bond order)
      for( size_t i( 0); i < n_atoms; ++i)
      {
        if( !atoms_could_change( i))
        {
          continue;
        }
        PossibleAtomTypesForAtom new_type_i
        (
          ATOMS( i).GetElementType(),
          ATOMS( i).GetAtomType()->GetNumberElectronsInBonds() - 1,
          ATOMS( i).GetAtomType()->GetNumberBonds(),
          0,
          false
        );
        // check whether bond order could be increased
        if( !new_type_i.GetNumberPossibleTypes() || new_type_i.GetMostStableType()->GetFormalCharge())
        {
          continue;
        }

        // check whether any pairs of atoms in add or change any bonded e are connected by anything other than a triple bond
        // if so, we can neutralize charge by increasing the bond order of the bond between them
        const storage::Vector< BondConformational> &bonds( ATOMS( i).GetBonds());

        size_t best_match( util::GetUndefined< size_t>());
        double best_electronegativity( -100.0);
        ConfigurationalBondType best_bond_type;
        AtomType best_atom_type;
        for( size_t bonded_atom_id( 0), n_bonds( bonds.GetSize()); bonded_atom_id < n_bonds; ++bonded_atom_id)
        {
          const size_t bonded_atom_index( ATOMS.GetAtomIndex( bonds( bonded_atom_id).GetTargetAtom()));

          if( atoms_could_change( bonded_atom_index) + atoms_could_change( bonded_atom_index))
          {
            continue;
          }

          // skip aromatic bonded atoms unless the preference allows changing it
          if
          (
            !allow_arom_chng
            && ATOMS( bonded_atom_index).CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsAromatic, 1)
          )
          {
            continue;
          }

          // skip single bonds
          const ConfigurationalBondType bond_type( bonds( bonded_atom_id).GetBondType());
          if( bond_type->GetBondData( ConfigurationalBondTypeData::e_BondOrder) == size_t( 1))
          {
            continue;
          }

          // test whether H-removal followed by increase of bond order would yield a valid type
          PossibleAtomTypesForAtom new_type_bonded
          (
            ATOMS( bonded_atom_index).GetElementType(),
            ATOMS( bonded_atom_index).GetAtomType()->GetNumberElectronsInBonds(),
            ATOMS( bonded_atom_index).GetAtomType()->GetNumberBonds() + 1,
            0,
            false
          );

          // no valid type if H is removed
          if( !new_type_bonded.GetNumberPossibleTypes() || new_type_bonded.GetMostStableType()->GetFormalCharge())
          {
            continue;
          }

          const double eneg
          (
            bonds( bonded_atom_id).GetTargetAtom().GetElementType()->GetProperty( ElementTypeData::e_ElectroNegativity)
          );
          if( eneg > best_electronegativity)
          {
            best_electronegativity = eneg;
            best_match = bonded_atom_index;
            best_bond_type = bonds( bonded_atom_id).GetBondType();
            best_atom_type = new_type_bonded.GetMostStableType();
          }
        }
        if( util::IsDefined( best_match))
        {
          did_something = true;
          atoms_could_change( i) = 0;
          atoms_could_change( best_match) = 0;
          BCL_MessageVrb
          (
            "NEUTRAL: decreased bond order ... new types are: " + new_type_i.GetMostStableType().GetName()
            + " " + best_atom_type.GetName() + " old types were: " + ATOMS( i).GetAtomType().GetName() + " " +
            ATOMS( best_match).GetAtomType().GetName()
          );
          ATOMS( i).SetBondTypeTo
          (
            ATOMS( best_match),
            best_bond_type->WithOrder( best_bond_type->GetBondData( ConfigurationalBondTypeData::e_BondOrder) - 1)
          );
          ATOMS( i).SetAtomType( new_type_i.GetMostStableType());
          ATOMS( best_match).SetAtomType( best_atom_type);
          atoms_add_h_to( best_match) = 1;
          atoms_add_h_to.PushBack( 0);
        }
      }

      // saturate partially if there were any atoms to add H's to
      if( !no_add_h && atoms_add_h_to.GetSize() > n_atoms)
      {
        HydrogensHandler::SaturatePartial( ATOMS, atoms_add_h_to);
        n_atoms = ATOMS.GetSize();
      }

      if( did_something)
      {
        AtomsCompleteStandardizer( ATOMS, "", true);
        if( !no_add_h)
        {
          HydrogensHandler::Saturate( ATOMS);
        }
      }
    }

  } // namespace chemistry
} // namespace bcl
