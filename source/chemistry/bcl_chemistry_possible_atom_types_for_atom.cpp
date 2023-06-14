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
#include "chemistry/bcl_chemistry_possible_atom_types_for_atom.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_bond_conformational.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PossibleAtomTypesForAtom::s_Instance( GetObjectInstances().AddInstance( new PossibleAtomTypesForAtom()));

    //! @brief Create the map from atom environment string to possible atom types
    //! @param IN_AROMATIC_RING whether to only include types that could be in an aromatic ring
    //! @param EXPLICIT_CHARGE true iff when set to 0 the expected charge must be interpreted literally as a neutral atom;
    //! default is false, which preserves backwards compatibility with old behavior where 0 allows searching of all atom types
    storage::Map< std::string, PossibleAtomTypesForAtom> PossibleAtomTypesForAtom::CreateAtomicEnvironmentToTypesMap
    (
      const bool IN_AROMATIC_RING,
      const bool EXPLICIT_CHARGE
    )
    {
      storage::Map< std::string, PossibleAtomTypesForAtom> atomic_environment_to_possible_types;
      // store the atomic environment for each atom type in a map for fast lookup when
      // determining atom types
      for
      (
        AtomTypes::const_iterator itr( GetAtomTypes().Begin()), itr_end( GetAtomTypes().End());
        itr != itr_end;
        ++itr
      )
      {
        if( !( *itr)->IsGasteigerAtomType()) // skip non-gasteiger types
        {
          continue;
        }

        AtomType type( *itr);

        if( !IN_AROMATIC_RING && EXPLICIT_CHARGE)
        {
          const std::string search_str
          (
            type->GetElementType()->GetChemicalSymbol()
            + util::Format()( type->GetNumberElectronsInBonds())
            + util::Format()( type->GetNumberBonds())
            + util::Format()( type->GetFormalCharge())
          );

          BCL_Debug( search_str);

          // add the type to the map, ignoring charge, since it is redundant with the # bonds and e- in bonds
          atomic_environment_to_possible_types[ search_str].AddAtomType( type);
        }
        else if( !IN_AROMATIC_RING)
        {
          const std::string search_str
          (
            type->GetElementType()->GetChemicalSymbol()
            + util::Format()( type->GetNumberElectronsInBonds())
            + util::Format()( type->GetNumberBonds())
          );

          BCL_Debug( search_str);

          // add the type to the map, ignoring charge, since it is redundant with the # bonds and e- in bonds
          atomic_environment_to_possible_types[ search_str].AddAtomType( type);
        }
        else if( type->GetNumberBonds() < size_t( 2))
        {
          // must be at least two bonds to be inside a ring
          continue;
        }
        else if( type->IsConjugated() && type->GetHybridOrbitalType() != GetHybridOrbitalTypes().e_Unhybridized)
        {
          // depending on what the atom in the aromatic ring is connected to, it may or may not have an extra
          // e- based considering aromatic bond types to have 3 e- in bonds
          const std::string search_str
          (
            type->GetElementType()->GetChemicalSymbol() + util::Format()( type->GetNumberBonds())
          );

          // add search strings for type in an aromatic ring with charges -1 - 1  to the map
          for
          (
            size_t min_double_bonds( 0), max_double_bonds( type->GetNumberElectronsInBonds() - type->GetNumberBonds());
            min_double_bonds <= max_double_bonds;
            ++min_double_bonds
          )
          {
            // the type may be chosen when an anion was requested if its charge is <= 0
            // an atom type with < 0 charge will be preferred, however
            if( type->GetFormalCharge() <= 0)
            {
              atomic_environment_to_possible_types[ search_str + util::Format()( min_double_bonds) + "N"].AddAromaticAtomType( type, short( -1));
            }

            // since charges are often omitted, all atom types are considered if a neutral atom type was requested
            // in a ring with declared aromatic bonds, unless EXPLICIT_CHARGE is true
            if( type->GetFormalCharge() == 0)
            {
              atomic_environment_to_possible_types[ search_str + util::Format()( min_double_bonds) + "O"].AddAromaticAtomType( type, short( 0));
            }

            // the type may be chosen when a cation was requestion if its charge is > 0
            // an atom type with > 0 charge will be preferred, however
            if( type->GetFormalCharge() >= 0)
            {
              atomic_environment_to_possible_types[ search_str + util::Format()( min_double_bonds) + "P"].AddAromaticAtomType( type, short( 1));
            }
          }
        }
      }

      // Remove unhybridized types wherever possible; the unhybridized types are only used when there is no
      // hybridized alternative
      for
      (
        storage::Map< std::string, PossibleAtomTypesForAtom>::iterator
          itr( atomic_environment_to_possible_types.Begin()),
          itr_end( atomic_environment_to_possible_types.End());
        itr != itr_end;
        ++itr
      )
      {
        if( IN_AROMATIC_RING)
        {
          // get the last character, which indicates the desired charge
          const char last_char( itr->first[ itr->first.size() - 1]);
          const short charge( last_char - 'O');
          itr->second.FinalizeAromatic( charge);
        }
        else
        {
          // finalize; only keep the best type
          itr->second.Finalize();
        }
      }

      return atomic_environment_to_possible_types;
    }

    //! @brief write out the atom typing scheme
    //! @param OSTREAM stream to write the atom typing scheme to
    std::ostream &PossibleAtomTypesForAtom::WriteDetailedScheme( std::ostream &OSTREAM)
    {
      const storage::Map< std::string, PossibleAtomTypesForAtom>
        atomic_environment_to_possible_types( CreateAtomicEnvironmentToTypesMap( false));
      OSTREAM << "Atom typing scheme except in rings with explicit aromatic bond types (sdf id=4)\n";
      OSTREAM << "Element\t# bonds\t# e- in bonds\tAtom Type(s)\n";
      for
      (
        storage::Map< std::string, PossibleAtomTypesForAtom>::const_iterator
          itr( atomic_environment_to_possible_types.Begin()), itr_end( atomic_environment_to_possible_types.End());
        itr != itr_end;
        ++itr
      )
      {
        const std::string element( itr->first.substr( 0, itr->first.size() - 2));
        const std::string n_e_in_bonds( 1, itr->first[ itr->first.size() - 2]);
        const std::string n_bonds( 1, itr->first[ itr->first.size() - 1]);
        OSTREAM << element << '\t' << n_bonds << '\t' << n_e_in_bonds;
        storage::Vector< AtomType> alternate_types( itr->second.GetAlternateTypes());
        for
        (
          storage::List< AtomType>::const_iterator
            itr_type( itr->second.m_AtomTypesByDecreasingStability.Begin()),
            itr_type_end( itr->second.m_AtomTypesByDecreasingStability.End());
          itr_type != itr_type_end;
          ++itr_type
        )
        {
          OSTREAM << '\t' << itr_type->GetName();
        }
        OSTREAM << '\n';
      }

      const storage::Map< std::string, PossibleAtomTypesForAtom>
        element_bonds_in_arom_ring_to_types_map( CreateAtomicEnvironmentToTypesMap( true));
      OSTREAM << "\n\nAtom Types chosen in aromatic rings with explicit aromatic bond types (sdf id=4)\n";
      OSTREAM << "# extra pi e- refers to the number of electrons in pi bonds except those in the aromatic bonds, eg\n"
              << "a double bond outside a ring would count as 1\n"
              << "Charge is one of the following -1,0,+1.  During atom typing, the charge may be changed to\n"
              << "preserve aromaticity\n"
              << "Element\t# bonds\t# extra pi e-\tCharge Sign\tAtom Type(s)\n";
      for
      (
        storage::Map< std::string, PossibleAtomTypesForAtom>::const_iterator
          itr( element_bonds_in_arom_ring_to_types_map.Begin()), itr_end( element_bonds_in_arom_ring_to_types_map.End());
        itr != itr_end;
        ++itr
      )
      {
        const std::string element( itr->first.substr( 0, itr->first.size() - 3));
        const std::string n_bonds( 1, itr->first[ itr->first.size() - 3]);
        const std::string n_extra_pi_electrons( 1, itr->first[ itr->first.size() - 2]);
        const int charge( int( itr->first[ itr->first.size() - 1]) - int( 'O'));
        OSTREAM << element << '\t' << n_bonds << '\t' << n_extra_pi_electrons << '\t' << charge;
        for
        (
          storage::List< AtomType>::const_iterator
            itr_type( itr->second.m_AtomTypesByDecreasingStability.Begin()),
            itr_type_end( itr->second.m_AtomTypesByDecreasingStability.End());
          itr_type != itr_type_end;
          ++itr_type
        )
        {
          OSTREAM << '\t' << itr_type->GetName();
        }
        OSTREAM << '\n';
      }
      return OSTREAM;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    PossibleAtomTypesForAtom::PossibleAtomTypesForAtom() :
      m_NumberAtomTypesWithHybridization( GetHybridOrbitalTypes().GetEnumCount(), size_t( 0)),
      m_AtomTypesByDecreasingStability(),
      m_NumberConjugatedTypes( 0),
      m_FinalizeFunction( NULL)
    {
    }

    //! @brief constructor from the known information about the atom
    //! @param ELEMENT element type,
    //! @param NUMBER_ELECTRONS_IN_BONDS number of electrons in bonds for the atom type; for atom in declared aromatic environment, number of exocyclic bonds
    //! @param NUMBER_BONDS number of bonds for the atom
    //! @param SUSPECTED_CHARGE; expected charge, ignored if no atom type matching the other criteria if found
    //! @param IN_AROMATIC_RING true iff the atom has bonds of the aromatic unspecified type
    //! @param EXPLICIT_CHARGE true iff when set to 0 the SUSPECTED_CHARGE must be interpreted literally as a neutral atom;
    //! default is false, which preserves backwards compatibility with old behavior where 0 allows searching of all atom types
    PossibleAtomTypesForAtom::PossibleAtomTypesForAtom
    (
      const ElementType &ELEMENT,
      const size_t NUMBER_ELECTRONS_IN_BONDS,
      const size_t NUMBER_BONDS,
      const short SUSPECTED_CHARGE,
      const bool IN_AROMATIC_RING,
      const bool EXPLICIT_CHARGE
    ) :
      m_NumberAtomTypesWithHybridization( GetHybridOrbitalTypes().GetEnumCount(), size_t( 0)),
      m_AtomTypesByDecreasingStability(),
      m_NumberConjugatedTypes( 0),
      m_FinalizeFunction( NULL)
    {
      // create maps from atomic environment to possible types
      static const storage::Map< std::string, PossibleAtomTypesForAtom>
        s_atomic_env_outside_arom_ring_to_types_map( CreateAtomicEnvironmentToTypesMap( false, EXPLICIT_CHARGE));
      static const storage::Map< std::string, PossibleAtomTypesForAtom>
        s_element_bonds_in_arom_ring_to_types_map( CreateAtomicEnvironmentToTypesMap( true, EXPLICIT_CHARGE));

      storage::Map< std::string, PossibleAtomTypesForAtom>::const_iterator itr;

      if( IN_AROMATIC_RING)
      {
        // when aromatic bonds (sdf id = 4) are given, 1 of 2 different numbers of electrons in bonds are possible,
        // depending on the charge.  Thus, first search for the type with the specified charge
        const std::string primary_search_string
        (
          ELEMENT->GetChemicalSymbol()
          + util::Format()( NUMBER_BONDS)
          + util::Format()( NUMBER_ELECTRONS_IN_BONDS - NUMBER_BONDS - 1)
          + std::string( 1, char( 'O' + SUSPECTED_CHARGE))
        );

        itr = s_element_bonds_in_arom_ring_to_types_map.Find( primary_search_string);
        if( itr != s_element_bonds_in_arom_ring_to_types_map.End())
        {
          // use the type with known charge
          *this = itr->second;
        }
        return;
      }


      const std::string primary_search_string
      (
        EXPLICIT_CHARGE ?
        ELEMENT->GetChemicalSymbol()
        + util::Format()( NUMBER_ELECTRONS_IN_BONDS)
        + util::Format()( NUMBER_BONDS) :
        ELEMENT->GetChemicalSymbol()
        + util::Format()( NUMBER_ELECTRONS_IN_BONDS)
        + util::Format()( NUMBER_BONDS)
//        + util::Format()( SUSPECTED_CHARGE)
      );
      // look for the type in the type-outside-ring map
      itr = s_atomic_env_outside_arom_ring_to_types_map.Find( primary_search_string);
      if( itr != s_atomic_env_outside_arom_ring_to_types_map.End())
      {
        *this = itr->second;
      }
    }

    PossibleAtomTypesForAtom *PossibleAtomTypesForAtom::Clone() const
    {
      return new PossibleAtomTypesForAtom( *this);
    }

  /////////////////
  // data access //
  /////////////////

    const std::string &PossibleAtomTypesForAtom::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    bool PossibleAtomTypesForAtom::CouldHaveHybridization( const HybridOrbitalType &HYBRID) const
    {
      return m_NumberAtomTypesWithHybridization( HYBRID.GetIndex());
    }

    size_t PossibleAtomTypesForAtom::GetNumberPossibleTypes() const
    {
      return m_AtomTypesByDecreasingStability.GetSize();
    }

    AtomType PossibleAtomTypesForAtom::GetMostStableType() const
    {
      return !m_AtomTypesByDecreasingStability.IsEmpty()
             ? m_AtomTypesByDecreasingStability.FirstElement() : GetAtomTypes().e_Undefined;
    }

    //! @brief get the alternate atom types
    //! @return the alternative atom types
    storage::Vector< AtomType> PossibleAtomTypesForAtom::GetAlternateTypes() const
    {
      return
        storage::Vector< AtomType>( ++m_AtomTypesByDecreasingStability.Begin(), m_AtomTypesByDecreasingStability.End());
    }

    //! @brief get the alternate atom type with the given charge
    //! @param CHARGE the charge desired
    //! @return an alternative atom type
    AtomType PossibleAtomTypesForAtom::GetAlternateTypeWithCharge( const short &CHARGE) const
    {
      for
      (
        storage::List< AtomType>::const_iterator
          itr( m_AtomTypesByDecreasingStability.Begin()), itr_end( m_AtomTypesByDecreasingStability.End());
        itr != itr_end;
        ++itr
      )
      {
        if( ( *itr)->GetFormalCharge() == CHARGE)
        {
          return *itr;
        }
      }
      return GetAtomTypes().e_Undefined;
    }

    bool PossibleAtomTypesForAtom::CouldBeConjugated() const
    {
      // return true if any conjugated types are available
      return m_NumberConjugatedTypes;
    }

    bool PossibleAtomTypesForAtom::MustBeConjugated() const
    {
      // return true if there is at least one conjugated type and the number of conjugated types
      // is the same as the number of types
      return m_NumberConjugatedTypes && m_NumberConjugatedTypes == m_AtomTypesByDecreasingStability.GetSize();
    }

    //! @brief determine the maximal # of pi-electrons in the pi-electron system
    //! @return the maximal # of pi-electrons in the pi-electron system
    size_t PossibleAtomTypesForAtom::GetMaxElectronsParticipatingInPiSystem() const
    {
      storage::List< AtomType>::const_iterator itr( m_AtomTypesByDecreasingStability.Begin()),
                                               itr_end( m_AtomTypesByDecreasingStability.End());

      while( itr != itr_end && !( *itr)->IsConjugated())
      {
        ++itr;
      }

      // if we found a non-conjugated type before the end, then return false, otherwise, return true
      return itr == itr_end ? 0 : ( *itr)->GetMaxEContributionToPiSystem();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief add an atom type to be considered
    //! @param ATOM_TYPE the type of atom to consider
    void PossibleAtomTypesForAtom::AddAtomType( const AtomType &ATOM_TYPE)
    {
      const HybridOrbitalType hybrid( ATOM_TYPE->GetHybridOrbitalType());

      // handle the common case where there are currently no atom types
      // also always add atom types for unhybridized types and atom types with only 1 bond
      // which are handled separately when Finalize is called
      if
      (
        m_AtomTypesByDecreasingStability.IsEmpty()
        || hybrid == GetHybridOrbitalTypes().e_Unhybridized
        || ATOM_TYPE->GetNumberBonds() == size_t( 1)
      )
      {
        ++m_NumberAtomTypesWithHybridization( hybrid.GetIndex());
        if( ATOM_TYPE->IsConjugated())
        {
          ++m_NumberConjugatedTypes;
        }
        m_AtomTypesByDecreasingStability.PushBack( ATOM_TYPE);
        return;
      }

      // handle the complicated case
      const double stability( ATOM_TYPE->GetStabilityMetric());

      if( CouldHaveHybridization( hybrid)) // find the type with that hybridization, test its stability
      {
        // test stability criterion to decide whether we should insert this hybrid
        storage::List< AtomType>::iterator itr( m_AtomTypesByDecreasingStability.Begin());

        // find the hybrid orbital type in the most stable type vector
        // we know this type is there, otherwise CouldHaveHybridization would have returned false
        while( ( *itr)->GetHybridOrbitalType() != hybrid)
        {
          ++itr;
        }

        // test the stability
        if( stability > ( *itr)->GetStabilityMetric())
        {
          // ATOM_TYPE is more stable than the existing type with that hybridization, so replace it
          if( ( *itr)->IsConjugated())
          {
            --m_NumberConjugatedTypes;
          }
          m_AtomTypesByDecreasingStability.RemoveElement( itr);
          --m_NumberAtomTypesWithHybridization( hybrid.GetIndex());
        }
      }

      // only insert if there are no types remaining with that hybridization
      if( !CouldHaveHybridization( hybrid))
      {
        ++m_NumberAtomTypesWithHybridization( hybrid.GetIndex());
        storage::List< AtomType>::iterator itr( m_AtomTypesByDecreasingStability.Begin()),
                                           itr_end( m_AtomTypesByDecreasingStability.End());

        while( itr != itr_end && ( *itr)->GetStabilityMetric() > stability)
        {
          ++itr;
        }

        m_AtomTypesByDecreasingStability.InsertElement( itr, ATOM_TYPE);
        if( ATOM_TYPE->IsConjugated())
        {
          ++m_NumberConjugatedTypes;
        }
      }
    }

    //! @brief set this object to only consider the given atom type
    //! @param ATOM_TYPE the atom type desired
    void PossibleAtomTypesForAtom::SetToType( const AtomType &ATOM_TYPE)
    {
      if( GetNumberPossibleTypes() > size_t( 0))
      {
        BCL_Assert
        (
          std::find
          (
            m_AtomTypesByDecreasingStability.Begin(),
            m_AtomTypesByDecreasingStability.End(),
            ATOM_TYPE
          ) != m_AtomTypesByDecreasingStability.End(),
          "Tried to set atom type to " + ATOM_TYPE.GetName()
          + " but only had these types: " + util::Format()( m_AtomTypesByDecreasingStability)
        );
      }
      // reset this object
      m_NumberAtomTypesWithHybridization.SetAllElements( 0);
      m_AtomTypesByDecreasingStability.Reset();
      m_NumberConjugatedTypes = 0;

      // add the atom type
      ++m_NumberAtomTypesWithHybridization( ATOM_TYPE->GetHybridOrbitalType().GetIndex());
      if( ATOM_TYPE->IsConjugated())
      {
        ++m_NumberConjugatedTypes;
      }
      m_AtomTypesByDecreasingStability.PushBack( ATOM_TYPE);
      m_FinalizeFunction = NULL;
    }

    //! @brief set the final type based on the given atom and smallest ring size
    //! @param ATOM the atom of interest
    //! @param SMALLEST_RING_SIZE the size of the smallest ring this atom is part of
    void PossibleAtomTypesForAtom::Finalize( const AtomConformationalInterface &ATOM, const size_t &SMALLEST_RING_SIZE)
    {
      // check whether there is any need for finalization
      if( m_FinalizeFunction)
      {
        // use the defined finalize function
        ( this->*m_FinalizeFunction)( ATOM, SMALLEST_RING_SIZE);
      }
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    std::istream &PossibleAtomTypesForAtom::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_AtomTypesByDecreasingStability, ISTREAM);
      m_NumberAtomTypesWithHybridization.SetAllElements( 0);
      m_NumberConjugatedTypes = 0;
      for
      (
        storage::List< AtomType>::const_iterator
          itr( m_AtomTypesByDecreasingStability.Begin()), itr_end( m_AtomTypesByDecreasingStability.End());
        itr != itr_end;
        ++itr
      )
      {
        ++m_NumberAtomTypesWithHybridization( ( *itr)->GetHybridOrbitalType().GetIndex());
        if( ( *itr)->IsConjugated())
        {
          ++m_NumberConjugatedTypes;
        }
      }

      return ISTREAM;
    }

    std::ostream &PossibleAtomTypesForAtom::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_AtomTypesByDecreasingStability, OSTREAM, INDENT);

      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief remove a particular hybrid orbital type from the possible types, unless that would remove all possibilities
    //! @param HYBRID the type of hybrid orbital to remove
    void PossibleAtomTypesForAtom::RemoveHybridization( const HybridOrbitalType &HYBRID)
    {
      // if there are no alternatives to HYBRID for this type, just return
      if( m_NumberAtomTypesWithHybridization( HYBRID.GetIndex()) == GetNumberPossibleTypes())
      {
        return;
      }

      // find the hybrid orbital type in the most stable type vector
      // we know this type is there, otherwise CouldHaveHybridization would have returned false
      for
      (
        storage::List< AtomType>::iterator
          itr( m_AtomTypesByDecreasingStability.Begin()), itr_end( m_AtomTypesByDecreasingStability.End());
        itr != itr_end;
        // iteration in loop
      )
      {
        if( ( *itr)->GetHybridOrbitalType() == HYBRID)
        {
          if( ( *itr)->IsConjugated())
          {
            --m_NumberConjugatedTypes;
          }
          itr = m_AtomTypesByDecreasingStability.Remove( itr);
          --m_NumberAtomTypesWithHybridization( HYBRID.GetIndex());
        }
        else // incorrect hybridization,
        {
          ++itr;
        }
      }
    }

    void PossibleAtomTypesForAtom::AddAromaticAtomType( const AtomType &ATOM_TYPE, const short &DESIRED_CHARGE)
    {
      const HybridOrbitalType hybrid( ATOM_TYPE->GetHybridOrbitalType());
      const short charge_diff( std::abs( ATOM_TYPE->GetFormalCharge() - DESIRED_CHARGE));

      // only add types with smaller than 2 difference in charge; changing charge by 2 would mean that the atom
      // was given a -1 charge and that we would replace it with a type that would have a +1 charge, or vice versa
      if( charge_diff > short( 1))
      {
        return;
      }

      ++m_NumberConjugatedTypes;
      ++m_NumberAtomTypesWithHybridization( hybrid.GetIndex());

      // handle the complicated case
      const double stability( ATOM_TYPE->GetStabilityMetric());

      // find out where this atom type should be placed
      storage::List< AtomType>::iterator itr( m_AtomTypesByDecreasingStability.Begin());

      static const storage::Vector< size_t> hybridization_rank
      (
        storage::Vector< size_t>::Create( 3, 2, 0, 1)
      );

      const size_t hybrid_rank( hybridization_rank( hybrid.GetIndex()));

      // place itr at the first atom type that ATOM_TYPE should be considered before in search order
      for
      (
        storage::List< AtomType>::const_iterator itr_end( m_AtomTypesByDecreasingStability.End());
        itr != itr_end;
        ++itr
      )
      {
        const short itr_charge_diff( std::abs( ( *itr)->GetFormalCharge() - DESIRED_CHARGE));

        if( charge_diff < itr_charge_diff)
        {
          break;
        }
        if( charge_diff > itr_charge_diff)
        {
          continue;
        }

        const size_t itr_hybrid_rank( hybridization_rank( ( *itr)->GetHybridOrbitalType().GetIndex()));
        if( hybrid_rank < itr_hybrid_rank)
        {
          break;
        }
        if( hybrid_rank > itr_hybrid_rank)
        {
          continue;
        }

        // equal hybridization ranks and charge difference; go for stability
        if( stability > ( *itr)->GetStabilityMetric())
        {
          break;
        }
      }
      m_AtomTypesByDecreasingStability.InsertElement( itr, ATOM_TYPE);
    }

    //! @brief Select the best choice for the atom type wherever possible
    //! @see @link https://structbio.vanderbilt.edu:8443/display/MeilerLab/RethinkingAtomTypeDetection @endlink
    //! @details the link above contains the statistics and models used to select the current set of rules
    void PossibleAtomTypesForAtom::Finalize()
    {
      if( m_AtomTypesByDecreasingStability.GetSize() == 1)
      {
        return;
      }

      // extract # bonds, element type, main group
      const int n_bonds( GetMostStableType()->GetNumberBonds());
      const ElementType element_type( GetMostStableType()->GetElementType());
      const int main_group( element_type->GetMainGroup());

      // handle group I and II elements; prefer unhybridized types with bonding S orbitals
      // handle group VII elements; prefer unhybridized types with bonding P orbitals
      if( n_bonds <= 2 && ( main_group == 1 || main_group == 2 || main_group == 7))
      {
        // Hybridization generally does not benefit these elements, which generally only form ionic bonds
        FinalizeUnhybridized();
        return;
      }

      // handle remaining elements.  These will all prefer to be hybridized, so remove unhybridized choices first
      RemoveHybridization( GetHybridOrbitalTypes().e_Unhybridized);
      if( m_AtomTypesByDecreasingStability.GetSize() == 1)
      {
        return;
      }

      // get # of electrons in bonds
      const int n_e_in_bonds( GetMostStableType()->GetNumberElectronsInBonds());

      // determine preferred hybridization
      HybridOrbitalType preferred_hybridization;

      if( n_bonds == 1) // handle remaining cases with a single bond
      {
        // Formally charged atoms with a single bond
        // choose hybrid orbital type according to VSEPR number (bonds + lone pairs)
        preferred_hybridization = HybridOrbitalType( 4 - n_e_in_bonds);
      }
      else if( n_bonds >= 4 || ( element_type->GetPeriod() >= size_t( 3) && n_e_in_bonds == n_bonds))
      {
        // with 4+ bonds, only Te types are available under the standard gasteiger types
        // period 3+ atom types (except group 1, 2, and 7) engage in d-orbital bonding, which is not represented
        // in gasteiger atom types.  Choose the most similar type by bond angles (SP3)
        preferred_hybridization = GetHybridOrbitalTypes().e_SP3;
      }
      else if( n_e_in_bonds - n_bonds >= 2)
      {
        // at least one triple bond or two double bonds
        // prefer digonal types if there are two bonds, trigonal if there are three
        if( n_bonds == 2)
        {
          preferred_hybridization = GetHybridOrbitalTypes().e_SP;
        }
        else // if( n_bonds == 3)
        {
          preferred_hybridization = GetHybridOrbitalTypes().e_SP2;
        }
      }
      else if( element_type->GetPeriod() >= size_t( 3))
      {
        // period 3+, group 3-6 elements with 2-3 bonds (with the current atom types, this can only be P, S, and Si),
        // exactly one double bond
        // period 3+ atom types (except group 1, 2, and 7) engage in d-orbital bonding, which is not represented
        // in gasteiger atom types.

        // under the assumption that similar bond angles reflect the most similar types, all types in this category
        // are closest to tetrahedral except Si with 3 bonds (1 double), which is closest to trigonal
        // On the other hand, none of the candidate atom types are tetrahedral since Te cannot form unsaturated bonds,
        // it is only through d-orbitals that this is achieved for period 3+ elements.
        // Likewise, prefer to eliminate all types except Tr and Te; the BCL will assert if any new types are added
        // that violate this rule.
        preferred_hybridization = GetHybridOrbitalTypes().e_SP3;
      }
      else if( n_bonds == 3 && n_e_in_bonds == 4)
      {
        // 1 unsaturated bond; two single bonds.  The pi-orbital assures planarity
        preferred_hybridization = GetHybridOrbitalTypes().e_SP2;
      }
      // down to B, C, N, O, with either 2 single bonds, 1 single & 1 double bond, or 3 single bonds
      else if( element_type == GetElementTypes().e_Boron || element_type == GetElementTypes().e_Carbon)
      {
        // O, and N atom types depend heavily on their environment
        // boron and carbon with 2 - 3 bonds are always trigonal
        preferred_hybridization = GetHybridOrbitalTypes().e_SP2;
      }
      // down to N and O, with either 2 single bonds, 1 single & 1 double bond, or 3 single bonds
      else if( element_type == GetElementTypes().e_Nitrogen)
      {
        if( n_bonds == 3 && n_e_in_bonds == 3)
        {
          m_FinalizeFunction = &PossibleAtomTypesForAtom::FinalizeNitrogenThreeSingle;
        }
        else if( n_bonds == 2 && n_e_in_bonds == 3)
        {
          m_FinalizeFunction = &PossibleAtomTypesForAtom::FinalizeNitrogenSingleDouble;
        }
        else if( n_bonds == 2 && n_e_in_bonds == 2)
        {
          m_FinalizeFunction = &PossibleAtomTypesForAtom::FinalizeNitrogenTwoSingle;
        }
      }
      else if( element_type == GetElementTypes().e_Oxygen)
      {
        if( n_bonds == 3 && n_e_in_bonds == 3)
        {
          m_FinalizeFunction = &PossibleAtomTypesForAtom::FinalizeOxygenThreeSingle;
        }
        else if( n_bonds == 2 && n_e_in_bonds == 3)
        {
          m_FinalizeFunction = &PossibleAtomTypesForAtom::FinalizeOxygenSingleDouble;
        }
        else if( n_bonds == 2 && n_e_in_bonds == 2)
        {
          m_FinalizeFunction = &PossibleAtomTypesForAtom::FinalizeOxygenTwoSingle;
        }
      }

      if( preferred_hybridization == GetHybridOrbitalTypes().e_SP)
      {
        RemoveHybridization( GetHybridOrbitalTypes().e_SP3);
        RemoveHybridization( GetHybridOrbitalTypes().e_SP2);
      }
      else if( preferred_hybridization == GetHybridOrbitalTypes().e_SP2)
      {
        RemoveHybridization( GetHybridOrbitalTypes().e_SP);
        RemoveHybridization( GetHybridOrbitalTypes().e_SP3);
      }
      else if( preferred_hybridization == GetHybridOrbitalTypes().e_SP3)
      {
        RemoveHybridization( GetHybridOrbitalTypes().e_SP);
        RemoveHybridization( GetHybridOrbitalTypes().e_SP2);
      }

      if( !m_FinalizeFunction && m_AtomTypesByDecreasingStability.GetSize() > size_t( 1))
      {
        // write ambiguous choices out
        std::stringstream equiv_choices;
        for
        (
          storage::List< AtomType>::const_iterator
            itr( m_AtomTypesByDecreasingStability.Begin()), itr_end( m_AtomTypesByDecreasingStability.End());
          itr != itr_end;
          ++itr
        )
        {
          equiv_choices << ' ' << itr->GetName();
        }
        BCL_Exit( "Unhandled atom types! " + equiv_choices.str(), -1);
      }
    }

    //! @brief choose the preferred atom type (using VSEPR theory) assuming that the orbitals do not hybridize
    //! @details This is used for elements in group 1, 2, & 7, which do hybridize in the gasteiger scheme
    void PossibleAtomTypesForAtom::FinalizeUnhybridized()
    {
      RemoveHybridization( GetHybridOrbitalTypes().e_SP3);
      RemoveHybridization( GetHybridOrbitalTypes().e_SP2);
      RemoveHybridization( GetHybridOrbitalTypes().e_SP);
      if( m_AtomTypesByDecreasingStability.GetSize() == 1)
      {
        return;
      }

      // ensure that only two types remain that differ only based on whether they bond with a sigma or just p orbitals
      BCL_Assert
      (
        m_AtomTypesByDecreasingStability.GetSize() == 2,
        "Should have been exactly two choices remaining for unhybridized types"
      );

      // choose the type with bonding sigma orbitals, based on normal valence electron / Lewis shell model
      AtomType atom_type_0_bonding_s_orbitals( m_AtomTypesByDecreasingStability.FirstElement());
      AtomType atom_type_1_bonding_s_orbital( m_AtomTypesByDecreasingStability.LastElement());
      if( atom_type_0_bonding_s_orbitals->GetNumberUnhybridizedSigmaOrbitals())
      {
        BCL_Assert
        (
          !atom_type_1_bonding_s_orbital->GetNumberUnhybridizedSigmaOrbitals(),
          "One of the two unhybridized types should have not had a bonding sigma orbital"
        );
        std::swap( atom_type_0_bonding_s_orbitals, atom_type_1_bonding_s_orbital);
      }
      else
      {
        BCL_Assert
        (
          atom_type_1_bonding_s_orbital->GetNumberUnhybridizedSigmaOrbitals(),
          "One of the two unhybridized types should have had a bonding sigma orbital"
        );
      }
      if( atom_type_0_bonding_s_orbitals->GetElementType()->GetMainGroup() == 7)
      {
        // the two hybridizations available are SP2P2P2 and S2P2P2P; prefer the latter,
        // based on normal valence electron / Lewis shell model
        SetToType( atom_type_0_bonding_s_orbitals);
      }
      else
      {
        // the hybridizations available are S/P or SP/PP.  Prefer the type with that bonds with s-orbitals (S/SP), which
        // are lower energy
        SetToType( atom_type_1_bonding_s_orbital);
      }
    }

    //! @brief only keep the most stable types for the atom that span the set of desired pi orbital electrons (0-2)
    //! @param DESIRED_CHARGE the preferred charge
    //! used during construction of the maps when there is no part of standardization that
    //! should edit this class
    void PossibleAtomTypesForAtom::FinalizeAromatic( const short &DESIRED_CHARGE)
    {
      // hybridization order vector maps hybridization to the order for that hybridization
      static const storage::Vector< size_t> hybridization_rank
      (
        storage::Vector< size_t>::Create( 3, 2, 0, 1)
      );

      // to work in any aromatic system, one of several different types may be necessary to satisfy hueckels rule
      // thus, always choose the first type in the list with the given number of pi electrons that has the best
      // hybridization order of any of the atom types with that hybridization
      storage::Map< size_t, AtomType> rank_to_best_type;

      // foreach x electrons in pi system
      //   - Prefer getting something closer to desired charge
      //   - If distance from desired charge is the same, prefer getting something with better hybridization rank
      //   - All else being equal, prefer the more stable compound
      for( size_t pi_electrons( 0), limit_pi_electrons( 3); pi_electrons < limit_pi_electrons; ++pi_electrons)
      {
        size_t best_hybrid_orbital_rank( 10);
        short best_charge_diff( 10);
        AtomType best_type( GetAtomTypes().e_Undefined);
        size_t index( 0), best_index( util::GetUndefined< size_t>());
        storage::List< AtomType> equivalent_choices;
        for
        (
          storage::List< AtomType>::const_iterator
            itr( m_AtomTypesByDecreasingStability.Begin()), itr_end( m_AtomTypesByDecreasingStability.End());
          itr != itr_end;
          ++itr, ++index
        )
        {
          if
          (
            ( *itr)->GetMaxEContributionToPiSystem() != pi_electrons
            &&
            ( pi_electrons != size_t( 0) || ( *itr)->GetPiElectronContributionType() != AtomTypeData::e_ZeroOrTwo)
          )
          {
            continue;
          }

          const short itr_charge_diff( std::abs( ( *itr)->GetFormalCharge() - DESIRED_CHARGE));
          const size_t hybrid_orbital_rank( hybridization_rank( ( *itr)->GetHybridOrbitalType().GetIndex()));
          if( hybrid_orbital_rank < best_hybrid_orbital_rank && itr_charge_diff <= best_charge_diff)
          {
            best_hybrid_orbital_rank = hybrid_orbital_rank;
            best_charge_diff = itr_charge_diff;
            best_type = *itr;
            best_index = index;
            equivalent_choices.Reset();
          }
          else if( hybrid_orbital_rank == best_hybrid_orbital_rank && itr_charge_diff == best_charge_diff)
          {
            equivalent_choices.PushBack( *itr);
          }
        }

        if( best_type.IsDefined())
        {
          rank_to_best_type[ best_index] = best_type;
        }
      }

      m_AtomTypesByDecreasingStability.Reset();
      for
      (
        storage::Map< size_t, AtomType>::const_iterator
          itr( rank_to_best_type.Begin()), itr_end( rank_to_best_type.End());
        itr != itr_end;
        ++itr
      )
      {
        m_AtomTypesByDecreasingStability.PushBack( itr->second);
      }

      m_NumberConjugatedTypes = m_AtomTypesByDecreasingStability.GetSize();
    }

    //! @brief choose the final atom type for Nitrogen with two single bonds
    //! @param ATOM the atom of interest
    //! @param SMALLEST_RING_SIZE the size of the smallest ring this atom is part of
    void PossibleAtomTypesForAtom::FinalizeNitrogenTwoSingle
    (
      const AtomConformationalInterface &ATOM,
      const size_t &SMALLEST_RING_SIZE
    )
    {
      // in ring of size 3 - 5 -> Tetrahedral (51 cases)
      // Bound to N or S -> Tetrahedral (90 cases)
      // Else -> Trigonal (53 cases)
      if( SMALLEST_RING_SIZE >= size_t( 3) && SMALLEST_RING_SIZE <= size_t( 5))
      {
        SetToType( GetAtomTypes().N_Te2Te2TeTe);
      }
      else
      {
        // check for bonds to N or S
        const storage::Set< ElementType> element_types( GetConnectedElementTypes( ATOM));
        if
        (
          element_types.Contains( GetElementTypes().e_Sulfur)
          || element_types.Contains( GetElementTypes().e_Nitrogen)
        )
        {
          SetToType( GetAtomTypes().N_Te2Te2TeTe);
        }
        else
        {
          SetToType( GetAtomTypes().N_Tr2TrTrPi2);
        }
      }
    }

    //! @brief choose the final atom type for Nitrogen with three single bonds
    //! @param ATOM the atom of interest
    //! @param SMALLEST_RING_SIZE the size of the smallest ring this atom is part of
    void PossibleAtomTypesForAtom::FinalizeNitrogenThreeSingle
    (
      const AtomConformationalInterface &ATOM,
      const size_t &SMALLEST_RING_SIZE
    )
    {
      // Node 1 In a 3 membered ring Tetrahedral (Accuracy: 100%, 293 cases)
      // Bond orbital overlap likely precludes trigonal
      if( SMALLEST_RING_SIZE == size_t( 3))
      {
        SetToType( GetAtomTypes().N_Te2TeTeTe);
        return;
      }
      const size_t unsaturated_neighbor_count( CountUnsaturatedNeighbors( ATOM));
      if( unsaturated_neighbor_count >= size_t( 2))
      {
        // Node 2 2+ unsaturated neighbors Trigonal (Accuracy: 99.5%, 11070 cases)
        // 2+ unsaturated neighbors guarantees that maximal orbital overlap will be achieved with trigonal geometry
        SetToType( GetAtomTypes().N_TrTrTrPi2);
        return;
      }
      else if( IsBondedToAHalogen( ATOM))
      {
        // Node 3 Connected to any halogen Tetrahedral (Accuracy: 88.9%, 45 cases)
        // Lone pair repulsion and strong electronegativity induces extra P character on N
        SetToType( GetAtomTypes().N_Te2TeTeTe);
        return;
      }
      // Node 4 1 unsaturated neighbor Trigonal (Accuracy: 94.1%, 14440 cases)
      // A single unsaturated neighbor is otherwise sufficient to grant trigonality in most cases
      // Exceptions are generally small (4-5 membered), fused or bridged rings, in which case bond angle strain can be severe, but presumably does not change the hybridization of N
      else if( unsaturated_neighbor_count == size_t( 1))
      {
        SetToType( GetAtomTypes().N_TrTrTrPi2);
        return;
      }

      storage::Set< ElementType> element_types( GetConnectedElementTypes( ATOM));

      if( element_types.Contains( GetElementTypes().e_Hydrogen))
      {
        element_types.Erase( GetElementTypes().e_Hydrogen);
      }
      // Node 6a All C neighbors Tetrahedral (Accuracy: 96.4%, 5981 cases) No impetus to become planar
      if( element_types.GetSize() == size_t( 1) && *element_types.Begin() == GetElementTypes().e_Carbon)
      {
        SetToType( GetAtomTypes().N_Te2TeTeTe);
        return;
      }

      // Node 5 Connected to any element period > 2 except those in group 6 Trigonal (Accuracy: 98.1%, 327 cases)
      // This is likely due to the relative similarity of d & p orbitals, since elements in period > 2 appear to engage
      // heavily in d-orbital bonding
      for
      (
        storage::Set< ElementType>::const_iterator
          itr( element_types.Begin()), itr_end( element_types.End());
        itr != itr_end;
        ++itr
      )
      {
        if( ( *itr)->GetPeriod() > size_t( 2) && ( *itr)->GetMainGroup() != size_t( 6))
        {
          // as a rule, H-O-H is tetrahedral; also holds for any halogens, N, and O (Node 2)
          SetToType( GetAtomTypes().N_TrTrTrPi2);
          return;
        }
      }

      // check for being in multiple, non-aromatic rings (aromatic rings were already checked
      // for by looking for unsaturated neighbors)
      if( ATOM.CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, size_t( 1)) == size_t( 3))
      {
        // Node 6b in multiple rings (Accuracy: 96.4%, 5981 cases) No impetus to become planar
        SetToType( GetAtomTypes().N_Te2TeTeTe);
      }
      else if( SMALLEST_RING_SIZE == size_t( 6))
      {
        // Node 7: In a 6 membered ring (Accuracy: 70.4%, 252 cases)
        if( element_types.Contains( GetElementTypes().e_Oxygen))
        {
          // Node 7a at least 1 O Trigonal (181 cases)
          SetToType( GetAtomTypes().N_TrTrTrPi2);
        }
        else
        {
          // Node 7b no O Tetrahedral (71 cases)
          SetToType( GetAtomTypes().N_Te2TeTeTe);
        }
      }
      else
      {
        // Node 8 Everything else (Accuracy: 78.6%, 564 cases)
        size_t heteroatom_count( 0);
        for
        (
          storage::Vector< BondConformational>::const_iterator
            itr( ATOM.GetBonds().Begin()), itr_end( ATOM.GetBonds().End());
          itr != itr_end;
          ++itr
        )
        {
          if( itr->GetTargetAtom().GetElementType() != GetElementTypes().e_Carbon && itr->GetTargetAtom().GetElementType() != GetElementTypes().e_Hydrogen)
          {
            ++heteroatom_count;
          }
        }
        if( heteroatom_count >= size_t( 2))
        {
          // Node 8a: Connected to 2+ heteroatoms Trigonal (71 cases)
          // Heteroatoms are generally quite electronegative and can form conjugated systems
          SetToType( GetAtomTypes().N_TrTrTrPi2);
        }
        else
        {
          // Node 8b: Connected to 1 heteroatom (493 cases)
          // Very little impetus to become planar
          SetToType( GetAtomTypes().N_Te2TeTeTe);
        }
      }
    }

    //! @brief choose the final atom type for a single and a double bond
    //! @param ATOM the atom of interest
    //! @param SMALLEST_RING_SIZE the size of the smallest ring this atom is part of
    void PossibleAtomTypesForAtom::FinalizeNitrogenSingleDouble
    (
      const AtomConformationalInterface &ATOM,
      const size_t &SMALLEST_RING_SIZE
    )
    {
      // use the digonal type only if bound to an atom w/ period 3+ and
      // whose main group (5 for N, 6 for O) is <= the bonded atom, and is at least 4
      // For example, Al should not induce digonality
      for
      (
        storage::Vector< BondConformational>::const_iterator
          itr( ATOM.GetBonds().Begin()), itr_end( ATOM.GetBonds().End());
        itr != itr_end;
        ++itr
      )
      {
        const ElementType &element( itr->GetTargetAtom().GetElementType());
        if( element->GetPeriod() > size_t( 2) && ( element->GetMainGroup() == 4 || element->GetMainGroup() == 5))
        {
          SetToType( GetAtomTypes().N_DiDiPi2Pi);
          return;
        }
      }
      SetToType( GetAtomTypes().N_Tr2TrTrPi);
    }

    //! @brief choose the final atom type for Oxygen with two single bonds
    //! @param ATOM the atom of interest
    //! @param SMALLEST_RING_SIZE the size of the smallest ring this atom is part of
    void PossibleAtomTypesForAtom::FinalizeOxygenTwoSingle
    (
      const AtomConformationalInterface &ATOM,
      const size_t &SMALLEST_RING_SIZE
    )
    {
      // any valences imply hydrogens; hydrogens imply tetrahedral
      if( ATOM.GetBonds().GetSize() < size_t( 2))
      {
        SetToType( GetAtomTypes().O_Te2Te2TeTe);
      }
      // check for aromatic ring membership
      // only the first bond needs to be checked
      else if( ATOM.GetBonds()( 0).GetBondType()->GetBondData( ConfigurationalBondTypeData::e_IsAromatic) == 1)
      {
        SetToType( GetAtomTypes().O_Tr2TrTrPi2);
      }
      // Node 1 In a 3-5 membered ring Tetrahedral (Accuracy: 99.5%, 10210 cases)
      // Bond orbital overlap likely precludes trigonal
      else if( SMALLEST_RING_SIZE >= size_t( 3) && SMALLEST_RING_SIZE <= size_t( 5))
      {
        SetToType( GetAtomTypes().O_Te2Te2TeTe);
      }
      else
      {
        storage::Set< ElementType> element_types( GetConnectedElementTypes( ATOM));
        if
        (
          element_types.Contains( GetElementTypes().e_Hydrogen)
          || element_types.Contains( GetElementTypes().e_Nitrogen)
          || element_types.Contains( GetElementTypes().e_Oxygen)
          || IsBondedToAHalogen( ATOM)
        )
        {
          // as a rule, H-O-H is tetrahedral; also holds for any halogens, N, and O (Node 2)
          SetToType( GetAtomTypes().O_Te2Te2TeTe);
        }
        // special case for sulfur, which needs to be checked for saturation
        else if( element_types.Contains( GetElementTypes().e_Sulfur))
        {
          bool had_saturated_s( false);
          for
          (
            storage::Vector< BondConformational>::const_iterator
              itr( ATOM.GetBonds().Begin()), itr_end( ATOM.GetBonds().End());
            itr != itr_end;
            ++itr
          )
          {
            const ElementType &element( itr->GetTargetAtom().GetElementType());
            if( element == GetElementTypes().e_Sulfur && !IsUnsaturated( itr->GetTargetAtom()))
            {
              // part of node 2
              SetToType( GetAtomTypes().O_Te2Te2TeTe);
              had_saturated_s = true;
            }
          }
          if( !had_saturated_s)
          {
            // other hetero atoms, choose the trigonal type (Node 3)
            // Node 3
            SetToType( GetAtomTypes().O_Tr2TrTrPi2);
          }
        }
        else if( element_types.GetSize() == size_t( 1) && *element_types.Begin() == GetElementTypes().e_Carbon)
        {
          // only bonded to carbon, have to count unsaturated neighbors
          if( CountUnsaturatedNeighbors( ATOM))
          {
            SetToType( GetAtomTypes().O_Tr2TrTrPi2);
          }
          else
          {
            SetToType( GetAtomTypes().O_Te2Te2TeTe);
          }
        }
        else if( element_types.GetSize() == size_t( 1) && *element_types.Begin() == GetElementTypes().e_Silicon)
        {
          // only connected to Si; Si - O - Si has strong digonal character (typically > 130 degrees)
          SetToType( GetAtomTypes().O_DiDiPi2Pi2);
        }
        else
        {
          // other hetero atoms, choose the trigonal type (Node 3)
          // Node 3
          SetToType( GetAtomTypes().O_Tr2TrTrPi2);
        }
      }
    }

    //! @brief choose the final atom type for Oxygen with a single and a double bond
    //! @param ATOM the atom of interest
    //! @param SMALLEST_RING_SIZE the size of the smallest ring this atom is part of
    void PossibleAtomTypesForAtom::FinalizeOxygenSingleDouble
    (
      const AtomConformationalInterface &ATOM,
      const size_t &SMALLEST_RING_SIZE
    )
    {
      // use the digonal type only if bound to an atom w/ period 3+ and
      // whose main group (5 for N, 6 for O) is <= the bonded atom, and is at least 4
      // For example, Al should not induce digonality
      for
      (
        storage::Vector< BondConformational>::const_iterator
          itr( ATOM.GetBonds().Begin()), itr_end( ATOM.GetBonds().End());
        itr != itr_end;
        ++itr
      )
      {
        const ElementType &element( itr->GetTargetAtom().GetElementType());
        const int main_group( element->GetMainGroup());
        if( element->GetPeriod() > size_t( 2) && ( main_group >= 4 && main_group <= 6))
        {
          SetToType( GetAtomTypes().O_DiDiPi2Pi);
          return;
        }
      }
      SetToType( GetAtomTypes().O_Tr2TrTrPi);
    }

    //! @brief choose the final atom type for Oxygen with three single bonds
    //! @param ATOM the atom of interest
    //! @param SMALLEST_RING_SIZE the size of the smallest ring this atom is part of
    void PossibleAtomTypesForAtom::FinalizeOxygenThreeSingle
    (
      const AtomConformationalInterface &ATOM,
      const size_t &SMALLEST_RING_SIZE
    )
    {
      // Trigonal if bound to anything group 1 - 3, else tetrahedral
      for
      (
        storage::Vector< BondConformational>::const_iterator
          itr( ATOM.GetBonds().Begin()), itr_end( ATOM.GetBonds().End());
        itr != itr_end;
        ++itr
      )
      {
        const ElementType &element( itr->GetTargetAtom().GetElementType());
        const int main_group( element->GetMainGroup());
        if( main_group > 0 && main_group < 4)
        {
          SetToType( GetAtomTypes().O_TrTrTrPi2);
          return;
        }
      }
      SetToType( GetAtomTypes().O_Te2TeTeTe);
    }

    //! @brief get connected element types
    //! @param ATOM the atom of interest
    //! @return a set of the connected element types
    storage::Set< ElementType> PossibleAtomTypesForAtom::GetConnectedElementTypes
    (
      const AtomConformationalInterface &ATOM
    )
    {
      storage::Set< ElementType> element_types;
      for
      (
        storage::Vector< BondConformational>::const_iterator
          itr( ATOM.GetBonds().Begin()), itr_end( ATOM.GetBonds().End());
        itr != itr_end;
        ++itr
      )
      {
        element_types.Insert( itr->GetTargetAtom().GetElementType());
      }
      if( ATOM.GetNumberofValenceBondsWithOrder( 1))
      {
        element_types.Insert( GetElementTypes().e_Hydrogen);
      }
      return element_types;
    }

    //! @brief test whether a particular atom is unsaturated
    //! @param ATOM the atom of interest
    //! @return true if atom has no A. unsaturated bonds or B. is part of an aromatic ring or C. has empty orbitals
    bool PossibleAtomTypesForAtom::IsUnsaturated( const AtomConformationalInterface &ATOM)
    {
      const size_t number_explicit_bonds( ATOM.GetBonds().GetSize());
      const size_t number_single_bonds
      (
        ATOM.CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_BondOrder, size_t( 1))
      );
      if( number_single_bonds == number_explicit_bonds)
      {
        // check for aromatic bonds
        const size_t number_aromatic_bonds
        (
          ATOM.CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsAromatic, size_t( 1))
        );
        if( number_aromatic_bonds == size_t( 0))
        {
          // no aromatic or unsaturated bonds; last resort, check for empty orbitals
          if( ATOM.GetElementType()->GetMainGroup() - ATOM.GetCharge() != 3)
          {
            return false;
          }
        }
      }
      return true;
    }

    //! @brief count unsaturated neighbors
    //! @param ATOM the atom of interest
    //! @return the number of unsaturated neighbors around ATOM
    size_t PossibleAtomTypesForAtom::CountUnsaturatedNeighbors( const AtomConformationalInterface &ATOM)
    {
      size_t unsaturated_counts( 0);
      for
      (
        storage::Vector< BondConformational>::const_iterator
          itr( ATOM.GetBonds().Begin()), itr_end( ATOM.GetBonds().End());
        itr != itr_end;
        ++itr
      )
      {
        if( IsUnsaturated( itr->GetTargetAtom()))
        {
          ++unsaturated_counts;
        }
      }
      return unsaturated_counts;
    }

    //! @brief test whether atom is bonded to any halogens
    //! @param ATOM the atom of interest
    //! @return true if the atom is bonded to any halogens
    bool PossibleAtomTypesForAtom::IsBondedToAHalogen( const AtomConformationalInterface &ATOM)
    {
      for
      (
        storage::Vector< BondConformational>::const_iterator
          itr( ATOM.GetBonds().Begin()), itr_end( ATOM.GetBonds().End());
        itr != itr_end;
        ++itr
      )
      {
        if( itr->GetTargetAtom().GetElementType()->GetMainGroup() == size_t( 7))
        {
          return true;
        }
      }
      return false;
    }

  } // namespace chemistry
} // namespace bcl

