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
#include "chemistry/bcl_chemistry_cs_code.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_base.h"
#include "descriptor/bcl_descriptor_iterator.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief virtual copy constructor
    CsCode *CsCode::Clone() const
    {
      return new CsCode( *this);
    }

    //! @brief returns class name
    //! the class name as const ref std::string
    const std::string &CsCode::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CsCode::Read( std::istream &ISTREAM)
    {
      // read member
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &CsCode::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write
      return OSTREAM;
    }

  ///////////////
  // operators //
  ///////////////

  ////////////////
  // operations //
  ////////////////

    //! @brief compares two elements of the atom environment
    //! @return true if left-hand side has higher bond type or atomic number (for same bond type)
    bool CsCode::CompareAtomsAndBonds
    (
      const storage::Pair< ConfigurationalBondType, util::SiPtr< const AtomConformationalInterface> > &LHS,
      const storage::Pair< ConfigurationalBondType, util::SiPtr< const AtomConformationalInterface> > &RHS
    )
    {
      return
           LHS.First()->GetBondData( ConfigurationalBondTypeData::e_BondOrderOrAromatic) > RHS.First()->GetBondData( ConfigurationalBondTypeData::e_BondOrderOrAromatic)
        || (
                LHS.First()->GetBondData( ConfigurationalBondTypeData::e_BondOrderOrAromatic) == RHS.First()->GetBondData( ConfigurationalBondTypeData::e_BondOrderOrAromatic)
             && LHS.Second()->GetElementType() > RHS.Second()->GetElementType()
           );
    }

    //! @brief generate the code vector for each element of the atom environment
    //! @return the code vector for each element of the atom environment
    storage::Vector< float> CsCode::GenerateCodeVector
    (
      const ConformationInterface &MOLECULE,
      const AtomConformationalInterface &ATOM,
      const ConfigurationalBondType &BOND_TYPE,
      const size_t NR_HYDROGENS,
      const float  RING_CLOSURE,
      storage::Vector< descriptor::CheminfoProperty> &PROPERTIES
    )
    {
      storage::Vector< float> code;

      code.PushBack( ATOM.GetElementType()->GetElectronConfiguration().ValenceElectronsSP());
      code.PushBack( ATOM.GetElementType()->GetPeriod());
      code.PushBack( ATOM.GetAtomType()->GetHybridOrbitalType().GetIndex());
      code.PushBack( BOND_TYPE->GetBondData( ConfigurationalBondTypeData::e_BondOrderOrAromatic));
      code.PushBack( NR_HYDROGENS);
      code.PushBack( RING_CLOSURE);

      const size_t index( MOLECULE.GetAtomIndex( ATOM));

      // append all the properties
      // construct a descriptor iterator to the atom
      descriptor::Iterator< AtomConformationalInterface> itr_atom( MOLECULE.GetAtomsIterator());
      itr_atom.GotoPosition( index);
      for
      (
        storage::Vector< descriptor::CheminfoProperty>::iterator
          itr_property( PROPERTIES.Begin()), itr_property_end( PROPERTIES.End());
        itr_property != itr_property_end;
        ++itr_property
      )
      {
        code.PushBack( ( **itr_property)( itr_atom)( 0));
      }
      code.PushBack( linal::Distance( ATOM.GetPosition(), ATOM.GetPosition()));

      return code;
    }

    //! @brief parse the miscellaneous property string containing the carbon chemical shift (NMRshiftDB format)
    //! @return map between the atom index and the according chemical shift
    storage::Map< size_t, float> CsCode::ParseChemicalShifts( const ConformationInterface &MOLECULE)
    {
      storage::Map< size_t, float> shifts;

      const std::string spectrum( MOLECULE.GetMDLProperty( "Spectrum 13C 0"));
      const storage::Vector< std::string> shift_strings( util::SplitString( spectrum, "|"));

      for
      (
        storage::Vector< std::string>::const_iterator
          itr_shifts( shift_strings.Begin()), itr_shifts_end( --shift_strings.End());
        itr_shifts != itr_shifts_end;
        ++itr_shifts
      )
      {
        const storage::Vector< std::string> shift_elements( util::SplitString( *itr_shifts, ";"));
        shifts
        [
          util::ConvertStringToNumericalValue< size_t>( shift_elements( 2))
        ] = util::ConvertStringToNumericalValue< float>( shift_elements( 0));
      }

      return shifts;
    }

    //! @brief implement a limited depth search from the atom of interest through the environment
    void CsCode::LimitedDepthFirstSearch
    (
      const ConformationInterface &MOLECULE,
      const AtomConformationalInterface &ATOM,
      const AtomEnvironment &ATOM_ENVIRONMENT,
      const size_t DEPTH,
      storage::Vector< float> &CODE,
      storage::Vector< descriptor::CheminfoProperty> &PROPERTIES,
      const bool &ONLY_FOLLOW_CONJUGATED_BONDS
    )
    {
      // which depth?
      BCL_MessageDbg( "Depth: " + util::Format()( DEPTH));
      BCL_MessageDbg( "++++++++++++++++++++++++++++++++++++++++++");
      BCL_MessageDbg
      (
        "Atom of interest #: " + util::Format()( MOLECULE.GetAtomIndex( ATOM))
      );

      // get sorted neighbors of atom with their according bonds
      storage::Vector< storage::Pair< ConfigurationalBondType, util::SiPtr< const AtomConformationalInterface> > >
        sorted_neighbors;

      for
      (
        storage::Vector< BondConformational>::const_iterator
          itr_bonds( ATOM.GetBonds().Begin()), itr_bonds_end( ATOM.GetBonds().End());
        itr_bonds != itr_bonds_end;
        ++itr_bonds
      )
      {
        sorted_neighbors.PushBack
        (
          storage::Pair< ConfigurationalBondType, util::SiPtr< const AtomConformationalInterface> >
          (
            itr_bonds->GetBondType(),
            itr_bonds->GetTargetAtom()
          )
        );
      }

      if( ONLY_FOLLOW_CONJUGATED_BONDS)
      {
        storage::Vector< storage::Pair< ConfigurationalBondType, util::SiPtr< const AtomConformationalInterface> > > conjugated_neighbors;
        for
        (
          storage::Vector< storage::Pair< ConfigurationalBondType, util::SiPtr< const AtomConformationalInterface> > >::const_iterator
            itr( sorted_neighbors.Begin()), itr_end( sorted_neighbors.End());
          itr != itr_end;
          ++itr
        )
        {
          if( itr->Second()->GetAtomType()->IsConjugated())
          {
            conjugated_neighbors.PushBack( *itr);
          }
        }
        conjugated_neighbors.InternalData().swap( sorted_neighbors.InternalData());
      }
      sorted_neighbors.Sort( &CompareAtomsAndBonds);

      size_t ring_closure( 0);

      // output necessary properties iterating through atoms
      for
      (
        storage::Vector< storage::Pair< ConfigurationalBondType, util::SiPtr< const AtomConformationalInterface> > >::const_iterator
          itr_pairs( sorted_neighbors.Begin()), itr_pairs_end( sorted_neighbors.End());
        itr_pairs != itr_pairs_end;
        ++itr_pairs
      )
      {
        // skip hydrogens in first sphere
        if
        (
             ATOM_ENVIRONMENT.GetEnvironmentAtoms().GetValue( ATOM).First() == 0
          && itr_pairs->Second()->GetElementType() == GetElementTypes().e_Hydrogen
        )
        {
          continue;
        }

        if
        (
          ATOM_ENVIRONMENT.GetEnvironmentAtoms().GetValue( itr_pairs->Second()).First() >
          ATOM_ENVIRONMENT.GetEnvironmentAtoms().GetValue( ATOM).First()
        )
        {
          BCL_MessageDbg
          (
            "now on atom #: " + util::Format()( MOLECULE.GetAtomIndex( *itr_pairs->Second()))
            + " atom of interest #: " + util::Format()( MOLECULE.GetAtomIndex( ATOM))
          );
          if( MOLECULE.GetAtomIndex( ATOM) == 10)
          {
            BCL_MessageDbg
            (
              "sorted_neighbors: " + util::Format()( sorted_neighbors)
            );
          }

          CODE.Append
          (
            GenerateCodeVector
            (
              MOLECULE,
              *itr_pairs->Second(),
              itr_pairs->First(),
              itr_pairs->Second()->GetNumberCovalentlyBoundHydrogens(),
              ATOM_ENVIRONMENT.GetEnvironmentAtoms().GetValue( itr_pairs->Second()).Second(),
              PROPERTIES
            )
          );

          // go one level deeper iterating through atoms
          if( DEPTH > 1)
          {
            LimitedDepthFirstSearch
            (
              MOLECULE,
              *itr_pairs->Second(),
              ATOM_ENVIRONMENT,
              DEPTH - 1,
              CODE,
              PROPERTIES,
              ONLY_FOLLOW_CONJUGATED_BONDS
            );
          }
        }

        // padding seems necessary if ring closure occurs
        else
        {
          ++ring_closure;
        }
      }

      // padding for missing substituents
      // only pad ATOM_ENVIRONMENT layer 2 or higher
      if( ATOM_ENVIRONMENT.GetEnvironmentAtoms().GetValue( ATOM).First() > 0)
      {
        BCL_MessageDbg
        (
          "padding neighbors on atom #: " + util::Format()( MOLECULE.GetAtomIndex( ATOM))
            + " for " + util::Format()( ( 4 - sorted_neighbors.GetSize())) + " * " + util::Format()( ( math::Pow< size_t>( 3, DEPTH) - 1) / 2
          ));

        // DEPTH determines how many atoms have to be padded beyond the actual atom

        // special case for conjugated systems and sphere 1 (atoms connected to atom of interest)
        // if the atom of interest is not conjugated itself, it's removed from sorted_neighbors
        // but still considered for the code generation, i.e., one less substituent needs to be padded
        if
        (
              ATOM_ENVIRONMENT.GetEnvironmentAtoms().GetValue( ATOM).First() == 1
          &&  !ATOM.GetAtomType()->IsConjugated()
          &&  ONLY_FOLLOW_CONJUGATED_BONDS
        )
        {
          CODE.Append
          (
            storage::Vector< float>
            (
              ( 3 - sorted_neighbors.GetSize()) * ( math::Pow< size_t>( 3, DEPTH) - 1) / 2 * 13, 0.0
            )
          );
        }
        else
        {
          CODE.Append
          (
            storage::Vector< float>
            (
              ( 4 - sorted_neighbors.GetSize()) * ( math::Pow< size_t>( 3, DEPTH) - 1) / 2 * 13, 0.0
            )
          );
        }
      }

      // padding for ring closure
//        if( ATOM_ENVIRONMENT.GetEnvironmentAtoms().GetValue( ATOM).Second() > 0)
//        {
//          CODE.Append( linal::Vector< float>( ( math::Pow< size_t>( 3, DEPTH) - 1) / 2 * 13, 0.0));
//        }
      for( size_t counter( 1); counter < ring_closure; ++counter)
      {
        BCL_MessageDbg( "padding ring closure: " + util::Format()( ( math::Pow< size_t>( 3, DEPTH) - 1) / 2));

        CODE.Append( storage::Vector< float>( ( math::Pow< size_t>( 3, DEPTH) - 1) / 2 * 13, 0.0));
      }
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace chemistry
} // namespace bcl

