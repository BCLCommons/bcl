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
#include "chemistry/bcl_chemistry_fragment_grow.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_merge_fragment_complete.h"
#include "find/bcl_find_collector_interface.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FragmentGrow::FragmentGrow() :
      m_FragmentPool(),
      m_AtomCollector(),
      m_AtomPicker(),
      m_FragmentPicker()
    {
    }

    //! @brief constructor from provided arguments
    //! @param FRAGMENT_POOL
    //! @param ATOM_COLLECTOR
    //! @param ATOM_PICKER
    //! @param FRAGMENT_PICKER
    FragmentGrow::FragmentGrow
    (
      const util::ShPtr< FragmentEnsemble> &FRAGMENT_POOL,
      const find::CollectorInterface< util::SiPtrList< const AtomConformationalInterface>, FragmentComplete> &ATOM_COLLECTOR,
      const find::PickInterface< util::SiPtr< const AtomConformationalInterface>, util::SiPtrList< const AtomConformationalInterface> > &ATOM_PICKER,
      const find::PickInterface< const FragmentComplete &, FragmentEnsemble> &FRAGMENT_PICKER
    ) :
      m_FragmentPool( FRAGMENT_POOL),
      m_AtomCollector( ATOM_COLLECTOR.Clone()),
      m_AtomPicker( ATOM_PICKER.Clone()),
      m_FragmentPicker( FRAGMENT_PICKER.Clone())
    {
    }

    //! @brief clone constructor
    FragmentGrow *FragmentGrow::Clone() const
    {
      return new FragmentGrow( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FragmentGrow::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief retrieves the fragment pool that this class uses by default
    //! @return a ShPtr to the fragment list
    util::ShPtr< FragmentEnsemble> FragmentGrow::GetFragmentPool() const
    {
      return m_FragmentPool;
    }

    //! @brief Set the fragment pool
    //! @param FRAGMENT_POOL the fragment pool to use
    void FragmentGrow::SetFragmentPool( const util::ShPtr< FragmentEnsemble> &FRAGMENT_POOL)
    {
      m_FragmentPool = FRAGMENT_POOL;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an fragment and generating a new fragment by growing on a valence
    //! @param FRAGMENT small molecule of interest
    //! @return MutateResult with Constitution after the mutate
    math::MutateResult< FragmentComplete> FragmentGrow::operator()( const FragmentComplete &FRAGMENT) const
    {

      // Call AddFragment
      util::ShPtr< FragmentComplete> new_molecule( GrowFragment( FRAGMENT));

      // return the new constitution
      return math::MutateResult< FragmentComplete>( new_molecule, *this);

    }

  ////////////////
  // operations //
  ////////////////

    //! @brief a function that adds a FragmentComplete to another FragmentComplete
    //! @param FRAGMENT the small molecule of interest
    //! @return FragmentComplete with the newly added fragment
    util::ShPtr< FragmentComplete> FragmentGrow::GrowFragment( const FragmentComplete &FRAGMENT) const
    {
      if( !m_FragmentPool.IsDefined())
      {
        return util::ShPtr< FragmentComplete>();
      }
      return AddFragmentFromList( FRAGMENT, *m_FragmentPool);
    }

    //! @brief create a copy of a molecule with a single random H removed
    //! @param MOLECULE the molecule to modify
    //! @return a new molecule with a single hydrogen removed
    FragmentComplete FragmentGrow::RemoveSingleH( const FragmentComplete &MOLECULE) const
    {
      AtomVector< AtomComplete> atom_v( MOLECULE.GetAtomVector());
      storage::Vector< size_t> h_atoms;
      h_atoms.AllocateMemory( atom_v.GetSize());

      for( size_t i( 0), end_i( atom_v.GetSize()); i < end_i; ++i)
      {
        if( atom_v( i).GetElementType() == GetElementTypes().e_Hydrogen)
        {
          h_atoms.PushBack( i);
        }
      }

      if( h_atoms.IsEmpty())
      {
        return MOLECULE;
      }

      h_atoms.Shuffle();
      storage::Vector< size_t> keep_atoms;
      keep_atoms.AllocateMemory( atom_v.GetSize() - 1);
      for( size_t i( 0), end_i( atom_v.GetSize()); i < end_i; ++i)
      {
        if( i != h_atoms( 0))
        {
          keep_atoms.PushBack( i);
        }
      }
      atom_v.Reorder( keep_atoms);
      return FragmentComplete( atom_v, MOLECULE.GetName());
    }

    //! @brief a function that adds a FragmentComplete to another FragmentComplete
    //! @param FRAGMENT the small molecule of interest
    //! @param FRAGMENT_LIST a list of fragments to use (overrides the internal fragment list)
    //! @return FragmentComplete with the newly added fragment
    util::ShPtr< FragmentComplete> FragmentGrow::AddFragmentFromList
    (
      const FragmentComplete &BASE_FRAGMENT,
      const FragmentEnsemble &FRAGMENT_LIST
    ) const
    {
      if( FRAGMENT_LIST.GetSize() == 0)
      {
        BCL_MessageStd( "fragment list is empty: cannot add a fragment");
        return util::ShPtr< FragmentComplete>();
      }

      BCL_MessageDbg( "starting fragment grow");

    ////////////////////
    // pick base atom //
    ////////////////////

      FragmentComplete base_fragment( BASE_FRAGMENT);

      util::SiPtrList< const AtomConformationalInterface> available_base_atoms( m_AtomCollector->Collect( base_fragment));
      if( available_base_atoms.IsEmpty())
      {
        base_fragment = RemoveSingleH( base_fragment);
        available_base_atoms = m_AtomCollector->Collect( base_fragment);
      }

      if( available_base_atoms.IsEmpty())
      {
        BCL_MessageStd( this->GetClassIdentifier() + " cannot grow fragment: no available sites on the parent molecule");
        return util::ShPtr< FragmentComplete>();
      }

      // pick the base valence atom
      util::SiPtr< const AtomConformationalInterface> picked_base_atom( m_AtomPicker->Pick( available_base_atoms));

    /////////////////////////////////////
    // pick fragment and fragment atom //
    /////////////////////////////////////

      // now pick the fragment to add to the picked_atom
      FragmentComplete picked_fragment( m_FragmentPicker->Pick( FRAGMENT_LIST));

      // if no valence atoms found, exit
      if( picked_fragment.GetNumberAtoms() == 0)
      {
        // warn user, return empty result
        BCL_MessageStd( "number of atoms in picked fragment is 0");
        return util::ShPtr< FragmentComplete>();
      }

      // collect possible grow fragment valence atoms
      util::SiPtrList< const AtomConformationalInterface> available_frag_atoms( m_AtomCollector->Collect( picked_fragment));
      if( available_frag_atoms.IsEmpty())
      {
        picked_fragment = RemoveSingleH( picked_fragment);
        available_frag_atoms = m_AtomCollector->Collect( picked_fragment);
      }

      if( available_frag_atoms.IsEmpty())
      {
        BCL_MessageStd( this->GetClassIdentifier() + " cannot grow fragment: no available sites on the fragment");
        return util::ShPtr< FragmentComplete>();
      }

      // picked valence atom in picked grow fragment
      util::SiPtr< const AtomConformationalInterface> picked_frag_atom( m_AtomPicker->Pick( available_frag_atoms));

    //////////
    // join //
    //////////

      return AddFragment( base_fragment, picked_fragment, picked_base_atom, picked_frag_atom);
    }

    //! @brief adds a fragment to a base molecule at a specific atom
    //! @param BASE the base molecule to add a fragment to
    //! @param FRAGMENT the piece to add to BASE
    //! @param BASE_ATOM the atom on BASE to connect to FRAGMENT
    //! @param FRAGMENT_ATOM the atom on FRAGMENT to connect to BASE
    //! @param MAX_BOND_ORDER the maximum multiplicity of the formed bond (0=max allowed)
    util::ShPtr< FragmentComplete> FragmentGrow::AddFragment
    (
      const FragmentComplete &BASE,
      const FragmentComplete &FRAGMENT,
      const util::SiPtr< const AtomConformationalInterface> &BASE_ATOM,
      const util::SiPtr< const AtomConformationalInterface> &FRAGMENT_ATOM,
      const size_t &MAX_BOND_ORDER
    )
    {
      size_t n_base_atoms( BASE.GetNumberAtoms());
      size_t n_frag_atoms( FRAGMENT.GetNumberAtoms());

      if( !n_base_atoms || !n_frag_atoms)
      {
        BCL_MessageStd( "Empty base molecule or fragment: cannot add fragment");
        return util::ShPtr< FragmentComplete>();
      }

      // Find BASE_ATOM in BASE; if it's not there then there is a problem
      size_t base_atom_index( BASE.GetAtomIndex( *BASE_ATOM));
      if( !util::IsDefined( base_atom_index))
      {
        BCL_MessageStd( "Could not find the specified base atom when adding a fragment");
        return util::ShPtr< FragmentComplete>();
      }
      // This should give the total available bond orders...?
      size_t base_valences( BASE_ATOM->GetNumberElectronsInValenceBonds());

      // Find FRAGMENT_ATOM in FRAGMENT; if it's not there then there is a problem
      size_t frag_atom_index( FRAGMENT.GetAtomIndex( *FRAGMENT_ATOM));
      if( !util::IsDefined( frag_atom_index))
      {
        BCL_MessageStd( "Could not find the specified fragment atom when adding a fragment");
        return util::ShPtr< FragmentComplete>();
      }
      size_t frag_valences( FRAGMENT_ATOM->GetNumberElectronsInValenceBonds());

      // Max bond order; either MAX_BOND_ORDER or the lowest common valences available
      size_t max_bond_order 
      (
        MAX_BOND_ORDER ? std::min( MAX_BOND_ORDER, std::min( base_valences, frag_valences))
                       : std::min( base_valences, frag_valences)
      ); 

      // what bond order should be formed, randomly chosen
      size_t desired_bond_order
      ( 
        max_bond_order ? random::GetGlobalRandom().Random< size_t>( max_bond_order)
                       : 1
      );

      while( !desired_bond_order)
      {
        desired_bond_order = random::GetGlobalRandom().Random< size_t>( max_bond_order);
      }
      
    //////////
    // join //
    //////////

      ConfigurationalBondType bond_type;
      switch( desired_bond_order)
      {
        case 1:
          {
            // Check if the new bond should be conjugated
            const bool is_conjugated_bond
            (
              FRAGMENT_ATOM->GetAtomType()->IsConjugated()
              && BASE_ATOM->GetAtomType()->IsConjugated()
            );
            bond_type = is_conjugated_bond
              ? GetConfigurationalBondTypes().e_ConjugatedSingleBond
              : GetConfigurationalBondTypes().e_NonConjugatedSingleBond;
          }
          break;
        case 2:
          bond_type = GetConfigurationalBondTypes().e_ConjugatedDoubleBond;
          break;
        case 3:
          bond_type = GetConfigurationalBondTypes().e_ConjugatedTripleBond;
          break;
        default:
          BCL_MessageStd( "Impossible bond order " + util::Format()( desired_bond_order));
          return util::ShPtr< FragmentComplete>();
          break;
      }

      // Add FRAGMENT to BASE
      storage::Pair< bool, FragmentComplete> new_fragment
      (
        MergeFragmentComplete::MergeFragments
        (
          BASE,
          FRAGMENT,
          bond_type,
          storage::Pair< size_t, size_t>( base_atom_index, frag_atom_index)
        )
      );

      // Standardize because things may be different now
      util::ShPtr< FragmentComplete> new_molecule;
      if( new_fragment.First())
      {
        AtomVector< AtomComplete> complete_vector( new_fragment.Second().GetAtomVector());
        const std::string name;
        AtomsCompleteStandardizer standardizer( complete_vector, name, true);
        new_molecule = util::ShPtr< FragmentComplete>( new FragmentComplete( complete_vector, ""));
      }
      return new_molecule;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FragmentGrow::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_FragmentPool, ISTREAM);
      io::Serialize::Read( m_AtomCollector, ISTREAM);
      io::Serialize::Read( m_AtomPicker, ISTREAM);
      io::Serialize::Read( m_FragmentPicker, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT number of indentations
    //! @return ostream which was written to
    std::ostream &FragmentGrow::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_FragmentPool, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AtomCollector, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AtomPicker, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_FragmentPicker, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace chemistry
} // namespace bcl
