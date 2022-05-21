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
#include "chemistry/bcl_chemistry_fragment_complete.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"
#include "chemistry/bcl_chemistry_bond_lengths.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_hydrogens_handler.h"
#include "chemistry/bcl_chemistry_stereocenters_handler.h"
#include "graph/bcl_graph_connectivity.h"
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_vector_3d_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> FragmentComplete::s_Instance
    (
      GetObjectInstances().AddInstance( new FragmentComplete())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    FragmentComplete::FragmentComplete()
    {
    }

    //! @brief constructor given atoms with conformation and molecule configuration
    //! @param ATOMS atom vector of atom completes
    //! @param NAME name of the fragment
    FragmentComplete::FragmentComplete
    (
      const AtomVector< AtomComplete> &ATOMS,
      const std::string &NAME
    ) :
      HasProperties< ConformationInterface>(),
      m_Atoms( ATOMS),
      m_Name( NAME)
    {
    }

    //! @brief constructor given atoms with conformation and molecule configuration
    //! @param ATOMS atom vector of atom completes
    //! @param NAME name of the fragment
    //! @param STORED_PROPERTIES stored properties of the fragment
    FragmentComplete::FragmentComplete
    (
      const AtomVector< AtomComplete> &ATOMS,
      const std::string &NAME,
      const storage::Map< std::string, std::string> &STORED_PROPERTIES
    ) :
      HasProperties< ConformationInterface>( STORED_PROPERTIES),
      m_Atoms( ATOMS),
      m_Name( NAME)
    {
      CacheNumeric( STORED_PROPERTIES);
    }

    //! @brief constructor with a constitution interface
    //! @params CONSTITUTION constitution interface
    FragmentComplete::FragmentComplete( const ConstitutionInterface &CONSTITUTION) :
      HasProperties< ConformationInterface>( CONSTITUTION),
      m_Name( CONSTITUTION.GetName())
    {
      m_Atoms = AtomVector< AtomComplete>( CONSTITUTION.GetAtomInfo(), CONSTITUTION.GetBondInfo());
      CacheNumeric( GetStoredProperties().GetMDLProperties());
    }

    //! @brief constructor with a configuration interface
    //! @params CONFIGURATION configuration interface
    FragmentComplete::FragmentComplete( const ConfigurationInterface &CONFIGURATION) :
      HasProperties< ConformationInterface>( CONFIGURATION),
      m_Name( CONFIGURATION.GetName())
    {
      m_Atoms = AtomVector< AtomComplete>( CONFIGURATION.GetAtomInfo(), CONFIGURATION.GetBondInfo());
      CacheNumeric( GetStoredProperties().GetMDLProperties());
    }

    //! @brief constructor with a conformation interface
    //! @params CONFORMATION conformation interface
    FragmentComplete::FragmentComplete( const ConformationInterface &CONFORMATION) :
      HasProperties< ConformationInterface>( CONFORMATION),
      m_Name( CONFORMATION.GetName())
    {
      m_Atoms = AtomVector< AtomComplete>( CONFORMATION.GetAtomInfo(), CONFORMATION.GetBondInfo());
      // no need to copy stored properties since the conformation interface's cache was already copied
    }

    //! @brief Clone function
    //! @return pointer to new FragmentComplete
    FragmentComplete *FragmentComplete::Clone() const
    {
      return new FragmentComplete( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &FragmentComplete::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns an iterator of atom conformations
    //! @return a generic iterator for atom conformations
    iterate::Generic< const AtomConformationalInterface> FragmentComplete::GetAtomsIterator() const
    {
      return iterate::Generic< const AtomConformationalInterface>( m_Atoms.Begin(), m_Atoms.End());
    }

    //! @brief return the number of atoms
    //! @return the number of atoms
    const size_t &FragmentComplete::GetNumberAtoms() const
    {
      return m_Atoms.GetSize();
    }

    //! @brief return the number of bonds
    //! @return the number of bonds
    const size_t &FragmentComplete::GetNumberBonds() const
    {
      return m_Atoms.GetNumberBonds();
    }

    //! @brief return the adjacency list
    //! @param BOND_SCHEME how to represent bond data as a size_t
    //! @return the adjacency list
    storage::Vector< graph::UndirectedEdge< size_t> >
      FragmentComplete::GetAdjacencyList( const ConfigurationalBondTypeData::Data &BOND_SCHEME) const
    {
      return m_Atoms.GetAdjacencyList( BOND_SCHEME);
    }

    //! @brief get bonds of molecule
    //! @return all the connectivities of a molecule
    storage::Vector< sdf::AtomInfo> FragmentComplete::GetAtomInfo() const
    {
      return m_Atoms.GetAtomInfo();
    }

    //! @brief get bonds of molecule
    //! @return all the connectivities of a molecule
    storage::Vector< sdf::BondInfo> FragmentComplete::GetBondInfo() const
    {
      return m_Atoms.GetBondInfo();
    }

    //! @brief get the index of a particular atom conformational interface
    //! @param ATOM the atom of interest
    //! @return the index of a particular atom conformational interface (undefined if atom is not in this molecule)
    size_t FragmentComplete::GetAtomIndex( const AtomConformationalInterface &ATOM) const
    {
      return m_Atoms.GetAtomIndex( ATOM);
    }

    //! @brief get the number of hydrogen atoms that preceed (i.e. occur sequentially prior to) the given atom
    //! @param ATOM the atom of interest
    //! @return the number of preceeding hydrogen atoms
    size_t FragmentComplete::GetNumberPreceedingHydrogenAtoms( const AtomConformationalInterface &ATOM) const
    {
      return GetNumberPreceedingHydrogenAtoms( GetAtomIndex( ATOM));
    }

    //! @brief get the number of hydrogen atoms that preceed (i.e. occur sequentially prior to) the given atom
    //! @param ATOM the atom index of interest
    //! @return the number of preceeding hydrogen atoms
    size_t FragmentComplete::GetNumberPreceedingHydrogenAtoms( const size_t ATOM) const
    {
      size_t h_count( 0);
      for( size_t a( 0); a < ATOM; ++a)
      {
        if( m_Atoms( a).GetElementType() == GetElementTypes().e_Hydrogen)
        {
          ++h_count;
        }
      }
      return h_count;
    }

    //! @brief get name of small molecule
    const std::string &FragmentComplete::GetName() const
    {
      return m_Name;
    }

    //! @brief set name of small molecule
    void FragmentComplete::SetName( const std::string &NAME)
    {
      m_Name = NAME;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FragmentComplete::Read( std::istream &ISTREAM)
    {

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &FragmentComplete::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {

      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief set coordinates of the molecule
    //! @param COORDINATES vector of coordinates of atoms( indices) of molecule
    //! @param INDICES vector of indices of molecule
    void FragmentComplete::SetAtomCoordinates
    (
      const storage::Vector< linal::Vector3D> &COORDINATES,
      const storage::Vector< size_t> &INDICES
    )
    {
      // iterate through the list of AtomsWithPosition and change their positions
      storage::Vector< linal::Vector3D>::const_iterator itr_coord( COORDINATES.Begin());
      // if the default empty vector of indices was given, check to see if there are enough coordinates
      // to set the positions of all atoms
      if( INDICES.IsEmpty() && COORDINATES.GetSize() == GetNumberAtoms())
      {
        for
        (
          AtomVector< AtomComplete>::iterator itr( m_Atoms.Begin()), itr_end( m_Atoms.End());
          itr != itr_end;
          ++itr, ++itr_coord
        )
        {
          // set the new position
          itr->SetPosition( *itr_coord);
        }
        return;
      }
      BCL_Assert( INDICES.GetSize() == COORDINATES.GetSize(), "Different number of sizes and indices given");
      for
      (
        storage::Vector< size_t>::const_iterator itr( INDICES.Begin()), itr_end( INDICES.End());
        itr != itr_end;
        ++itr, ++itr_coord
      )
      {
        // set the new position
        m_Atoms( *itr).SetPosition( *itr_coord);
      }
      ConformationInterface::GetChangeSignal().Emit( *this);
    }

    //! @brief try to generate an initial idealized geometry. This will fail miserably for rings but generate something
    //!        reasonable (but potentially clashing) otherwise. Assumes only a single molecule (all atoms reachable from all
    //!        other atoms via covalent bonds, not a complex)
    void FragmentComplete::IdealizeGeometry()
    {
      ResetCache();
      if( m_Atoms.IsEmpty())
      {
        return;
      }
      // set the first atom at the origin
      m_Atoms( 0).SetPosition( linal::Vector3D( 0.0));

      // set all other atoms up undefined positions
      static const linal::Vector3D s_undefined( util::GetUndefinedDouble());
      for
      (
        auto itr( m_Atoms.Begin() + 1), itr_end( m_Atoms.End());
        itr != itr_end;
        ++itr
      )
      {
        itr->SetPosition( s_undefined);
      }
      storage::Vector< size_t> have_seen( m_Atoms.GetSize(), size_t( 0));
      have_seen( 0) = 1;
      storage::List< size_t> queue;
      queue.PushBack( 0);

      // iterate over queue
      while( !queue.IsEmpty())
      {
        size_t el( queue.FirstElement());
        HydrogensHandler::SetUndefinedHCoordinates( m_Atoms( el), m_Atoms);
        for
        (
          auto itr_bnd( m_Atoms( el).GetBonds().Begin()), itr_bnd_end( m_Atoms( el).GetBonds().End());
          itr_bnd != itr_bnd_end;
          ++itr_bnd
        )
        {
          size_t ind( m_Atoms.GetAtomIndex( itr_bnd->GetTargetAtom()));
          if( !have_seen( ind))
          {
            have_seen( ind) = 1;
            queue.PushBack( ind);
          }
        }
        queue.PopFront();
      }
      // add isometry info
      BondIsometryHandler::AddIsometryInformation( m_Atoms, true);

      // add stereocenter information
      StereocentersHandler::AddChiralityFromConformation( m_Atoms);
      ConformationInterface::GetChangeSignal().Emit( *this);
    }

    //! @brief Remove all hydrogens from the molecule
    void FragmentComplete::RemoveH()
    {
      if( GetNumberHydrogens())
      {
        ResetCache();
        HydrogensHandler::Remove( m_Atoms);
        ConformationInterface::GetChangeSignal().Emit( *this);
      }
    }

    //! @brief Saturate all single bond valences with H
    void FragmentComplete::SaturateWithH()
    {
      if( GetNumberValences())
      {
        ResetCache();
        HydrogensHandler::Saturate( m_Atoms);
        ConformationInterface::GetChangeSignal().Emit( *this);
      }
    }

    //! @brief Encode resonance
    void FragmentComplete::EncodeTerminalResonance()
    {
      bool removed_h( false);
      if( GetNumberHydrogens())
      {
        BCL_MessageCrt( "Removing all hydrogens for encoding resonance. Hopefully that is fine for your use case");
        RemoveH();
        removed_h = true;
      }
      for( auto itr( m_Atoms.Begin()), itr_end( m_Atoms.End()); itr != itr_end; ++itr)
      {
        // consider only atoms with double or triple bond to neighboring atom
        if
        (
          itr->GetAtomType()->GetNumberBonds() != size_t( 1)
          || itr->GetAtomType()->GetNumberElectronsInBonds() != size_t( 2)
          || itr->GetBonds().IsEmpty()
          || itr->GetCharge()
        )
        {
          continue;
        }

        AtomComplete &atom_a( *itr);
        const AtomConformationalInterface &atom_b_const( itr->GetBonds().Begin()->GetTargetAtom());
        if( atom_b_const.GetBonds().GetSize() == size_t( 1))
        {
          continue;
        }
        auto itr_c( atom_b_const.GetBonds().Begin()), itr_c_end( atom_b_const.GetBonds().End());
        for( ; itr_c != itr_c_end; ++itr_c)
        {
          // ignore bonds of the same order. This ignores atom_a as a side-effect
          if( itr_c->GetBondType()->GetNumberOfElectrons() / 2 >= atom_a.GetAtomType()->GetNumberElectronsInBonds())
          {
            continue;
          }
          if( itr_c->GetTargetAtom().GetElementType() != atom_a.GetElementType())
          {
            continue;
          }
          if( itr_c->GetTargetAtom().GetBonds().GetSize() != size_t( 1))
          {
            continue;
          }
          break;
        }
        if( itr_c == itr_c_end)
        {
          continue;
        }
        AtomComplete &atom_b( m_Atoms( m_Atoms.GetAtomIndex( atom_b_const)));
        AtomComplete &atom_c( m_Atoms( m_Atoms.GetAtomIndex( itr_c->GetTargetAtom())));

        PossibleAtomTypesForAtom new_type_b_bond_promotion
        (
          atom_b.GetElementType(),
          atom_b.GetAtomType()->GetNumberElectronsInBonds() + 1,
          atom_b.GetAtomType()->GetNumberBonds(),
          0,
          false
        );
        if( new_type_b_bond_promotion.GetNumberPossibleTypes() == size_t( 0))
        {
          PossibleAtomTypesForAtom new_type_b_bond_demotion
          (
            atom_b.GetElementType(),
            atom_b.GetAtomType()->GetNumberElectronsInBonds() - 1,
            atom_b.GetAtomType()->GetNumberBonds(),
            0,
            false
          );
          if( new_type_b_bond_demotion.GetNumberPossibleTypes() == size_t( 0))
          {
            PossibleAtomTypesForAtom new_type_b_bond_demo_w_h
            (
              atom_b.GetElementType(),
              atom_b.GetAtomType()->GetNumberElectronsInBonds(),
              atom_b.GetAtomType()->GetNumberBonds() + 1,
              0,
              false
            );
            if( new_type_b_bond_demo_w_h.GetNumberPossibleTypes() == size_t( 0))
            {
              BCL_MessageCrt( "Could not encode terminal resonance!");
              continue;
            }
            else
            {
              // demote the double bond between atom_b to atom_a to a single bond
              BCL_MessageVrb
              (
                "RESONANCE: lowered bond order and added H... new types are: " + new_type_b_bond_demo_w_h.GetMostStableType().GetName()
                + " " + atom_c.GetAtomType().GetName()
              );
              atom_b.SetBondTypeTo
              (
                atom_a,
                itr_c->GetBondType()
              );
              atom_b.SetAtomType( new_type_b_bond_demo_w_h.GetMostStableType());
              atom_a.SetAtomType( atom_c.GetAtomType());
            }
          }
          else
          {
            // demote the double bond between atom_b to atom_a to a single bond
            BCL_MessageVrb
            (
              "RESONANCE: lowered bond order ... new types are: " + new_type_b_bond_demotion.GetMostStableType().GetName()
              + " " + atom_c.GetAtomType().GetName()
            );
            atom_b.SetBondTypeTo
            (
              atom_a,
              itr_c->GetBondType()
            );
            atom_b.SetAtomType( new_type_b_bond_demotion.GetMostStableType());
            atom_a.SetAtomType( atom_c.GetAtomType());
          }
        }
        else
        {
          // promote the single bond between atom_b to atom_c to a double bond
          BCL_MessageVrb
          (
            "RESONANCE: upped bond order ... new types are: " + new_type_b_bond_promotion.GetMostStableType().GetName()
            + " " + atom_a.GetAtomType().GetName()
          );
          atom_b.SetBondTypeTo
          (
            atom_c,
            itr_c->GetBondType()->WithOrder( itr_c->GetBondType()->GetBondData( ConfigurationalBondTypeData::e_BondOrder) + 1)
          );
          atom_b.SetAtomType( new_type_b_bond_promotion.GetMostStableType());
          atom_c.SetAtomType( atom_a.GetAtomType());
        }
      }
      if( removed_h)
      {
        SaturateWithH();
      }
    }

    //! @brief Canonicalize the order of atoms in the molecule, such that the atoms are arranged in descending CIP priority
    void FragmentComplete::Canonicalize()
    {
      storage::Vector< SubstituentConformational> all_substituents;
      all_substituents.AllocateMemory( m_Atoms.GetSize());
      for
      (
        AtomVector< AtomComplete>::const_iterator itr( m_Atoms.Begin()), itr_end( m_Atoms.End());
        itr != itr_end;
        ++itr
      )
      {
        all_substituents.PushBack( SubstituentConformational( *itr));
      }

      //Sort SubstituentConformational vector and check that they are properly sorted in descending priority
      std::stable_sort( all_substituents.Begin(), all_substituents.End());

      // get the new order
      storage::Vector< size_t> new_order;
      new_order.AllocateMemory( m_Atoms.GetSize());
      for
      (
        storage::Vector< SubstituentConformational>::const_iterator
          itr( all_substituents.Begin()), itr_end( all_substituents.End());
        itr != itr_end;
        ++itr
      )
      {
        new_order.PushBack( m_Atoms.GetAtomIndex( *itr->GetRootAtom()));
      }

      m_Atoms.Reorder( new_order);
    }

    //! @brief Standardize non-ring bond lengths according to covalent radii
    //! @param MOBILE_ATOMS atoms that can be perturbed during standardization
    void FragmentComplete::StandardizeBondLengths( const storage::Vector< size_t> &MOBILE_ATOMS)
    {
      // so we can reference atoms in mask
      const AtomVector< AtomComplete> &atom_v( GetAtomVector());

      // assign all atoms as mobile if nothing was passed
      storage::Vector< size_t> mobile_atoms;
      if( MOBILE_ATOMS.IsEmpty())
      {
        for( size_t i( 0); i < atom_v.GetSize(); ++i)
        {
          mobile_atoms.PushBack( i);
        }
      }
      else
      {
        mobile_atoms = MOBILE_ATOMS;
      }

      // bond graph
      auto boarom( ConfigurationalBondTypeData::e_BondOrderOrAromatic);
      graph::ConstGraph< size_t, size_t> molgraph
      (
        ConformationGraphConverter( ConformationGraphConverter::e_Identity, boarom)( *this)
      );

      // iterate over all atoms
      size_t i( 0);
      for( iterate::Generic< const AtomConformationalInterface> itr( GetAtomsIterator()); itr.NotAtEnd(); ++itr, ++i)
      {
        for
        (
            auto itr_bonds( itr->GetBonds().Begin()), itr_bonds_end( itr->GetBonds().End());
            itr_bonds != itr_bonds_end;
            ++itr_bonds
        )
        {
          // set condition on checking order
          // TODO: make more performant by replacing index vector with binary vector for masking
          bool check( false);
          if
          (
              mobile_atoms.Find( atom_v.GetAtomIndex( *itr)) < mobile_atoms.GetSize() && // is mobile
              mobile_atoms.Find( atom_v.GetAtomIndex( itr_bonds->GetTargetAtom())) < mobile_atoms.GetSize() // is mobile
          )
          {
            check = true;
          }

          if( itr_bonds->GetBondType()->IsBondInRing() || ( &*itr > &itr_bonds->GetTargetAtom() && check))
          {
            continue;
          }
          double desired_length
          (
            BondLengths::GetBondLength
            (
              itr->GetAtomType(),
              itr_bonds->GetBondType()->GetBondData( boarom),
              itr_bonds->GetTargetAtom().GetAtomType()
            )
          );
          if( !util::IsDefined( desired_length))
          {
            BCL_MessageStd( "Undefined bond length during standardization: " + itr->GetAtomType().GetName() + " " +
              itr_bonds->GetBondType().GetName() + " " + itr_bonds->GetTargetAtom().GetAtomType().GetName());
            continue;
          }
          linal::Vector3D existing_bond_direction( itr_bonds->GetTargetAtom().GetPosition() - itr->GetPosition());
          const double norm( existing_bond_direction.Norm());
          if( norm)
          {
            existing_bond_direction *= desired_length / norm - 1.0;
          }
          else
          {
            existing_bond_direction = linal::Vector3D( desired_length, 0.0, 0.0);
          }
          // existing bond direction is now the translation to add to everything on the other side of the bond
          auto shptr_vectices_reachable
          (
            graph::Connectivity::GetVerticesReachableFromDirectedEdge
            (
              molgraph,
              m_Atoms.GetAtomIndex( itr_bonds->GetTargetAtom()),
              i
            )
          );
          for
          (
              auto itr_downstream( shptr_vectices_reachable->Begin()), itr_downstream_end( shptr_vectices_reachable->End());
              itr_downstream != itr_downstream_end;
              ++itr_downstream
          )
          {
            if
            (
                *itr_downstream != i && // not self atom
                mobile_atoms.Find( *itr_downstream) < mobile_atoms.GetSize() // is mobile
            )
            {
              m_Atoms( *itr_downstream).SetPosition( m_Atoms( *itr_downstream).GetPosition() + existing_bond_direction);
            }
          }
        }
      }
      ConformationInterface::GetChangeSignal().Emit( *this);
    }

    //! @brief Update all hydrogen positions
    void FragmentComplete::UpdateH()
    {
      for( auto itr( m_Atoms.Begin()), itr_end( m_Atoms.End()); itr != itr_end; ++itr)
      {
        if( itr->GetNumberCovalentlyBoundHydrogens())
        {
          HydrogensHandler::UpdateHCoordinates( *itr, m_Atoms);
        }
      }
    }

    //! @brief merge a given fragment with this fragment
    //! @param FRAGMENT the given fragment
    //! @param BONDS_BETWEEN interconnectivity between this and the given fragment where first member of triplet is
    //!        is index of atom of this fragment while third member is index of atom in the given fragment.
    void FragmentComplete::MergeWithFragment
    (
      const FragmentComplete &FRAGMENT,
      const storage::Vector< sdf::BondInfo> &BONDS_BETWEEN
    )
    {
      ResetCache();
      m_Atoms.AddAtomsWithConnectivity( FRAGMENT.m_Atoms, BONDS_BETWEEN);
    }

    //! @brief get changable atoms conformational interface
    //! @return iterator to changable atoms conformational interface, which allows the base class to call SetPosition
    iterate::Generic< AtomConformationalInterface> FragmentComplete::GetAtomsIteratorNonConst()
    {
      return iterate::Generic< AtomConformationalInterface>( m_Atoms.Begin(), m_Atoms.End());
    }

  } // namespace chemistry
} // namespace bcl
