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
#include "chemistry/bcl_chemistry_atom_environment.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_bond_conformational.h"
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "util/bcl_util_si_ptr_list.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    AtomEnvironment::AtomEnvironment()
    {
    }

    //! construct from number of spheres
    AtomEnvironment::AtomEnvironment
    (
      const AtomConformationalInterface &ATOM_OF_INTEREST,
      const size_t NUMBER_SPHERES
    ) :
      m_AtomOfInterest( ATOM_OF_INTEREST)
    {
      // construct sphere 0 of the atom environment containing only the 'atom of interest'
      m_EnvironmentAtoms[ m_AtomOfInterest] = storage::Pair< size_t, double>( 0, 0);

      // keep track of the last sphere to iterate through constructing the next sphere
      util::SiPtrList< const AtomConformationalInterface> last_sphere;
      last_sphere.PushBack( m_AtomOfInterest);

      // construct spheres 1..n by iterating through the last sphere and try to add all new atoms
      // if atoms are already in the environment: ring closure
      for( size_t count_spheres( 1); count_spheres < NUMBER_SPHERES; ++count_spheres)
      {
        util::SiPtrList< const AtomConformationalInterface> next_sphere;

        // construct each sphere of the substituent cone
        for
        (
          util::SiPtrList< const AtomConformationalInterface>::const_iterator
            itr_atom( last_sphere.Begin()), itr_atom_end( last_sphere.End());
          itr_atom != itr_atom_end;
          ++itr_atom
        )
        {
          for
          (
            storage::Vector< BondConformational>::const_iterator
              itr_bonds( ( *itr_atom)->GetBonds().Begin()),
              itr_bonds_end( ( *itr_atom)->GetBonds().End());
            itr_bonds != itr_bonds_end;
            ++itr_bonds
          )
          {
            if( !m_EnvironmentAtoms.Has( itr_bonds->GetTargetAtom()))
            {
              m_EnvironmentAtoms[ itr_bonds->GetTargetAtom()] = storage::Pair< size_t, double>( count_spheres, 0);
              next_sphere.PushBack( itr_bonds->GetTargetAtom());
            }
            else if
            (
              m_EnvironmentAtoms[ itr_bonds->GetTargetAtom()].First() == count_spheres - 1
            )
            {
              m_EnvironmentAtoms[ itr_bonds->GetTargetAtom()].Second() = 1.0;
              m_EnvironmentAtoms[ *itr_atom].Second() = 1.0;
            }
            else if
            (
              m_EnvironmentAtoms[ itr_bonds->GetTargetAtom()].First() == count_spheres
            )
            {
              m_EnvironmentAtoms[ itr_bonds->GetTargetAtom()].Second() = 0.5;
            }
          }
        }

        last_sphere = next_sphere;
      }
    }

    //! virtual copy constructor
    AtomEnvironment *AtomEnvironment::Clone() const
    {
      return new AtomEnvironment( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AtomEnvironment::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief compute atom environment two bonds from reference atom
    //! @param A atom of interest
    //! @return return string indicating atom environment
    std::string AtomEnvironment::MakeAtomEnvironmentStringTwo( const AtomConformationalInterface &A)
    {
      // first layer
      storage::Vector< std::string> subs;
      subs.AllocateMemory( A.GetBonds().GetSize());
      for
      (
          auto itr( A.GetBonds().Begin()), itr_end( A.GetBonds().End());
          itr != itr_end;
          ++itr
      )
      {
        std::string c;
        c += '0' + itr->GetBondType()->GetBondData( ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness);
        c += itr->GetTargetAtom().GetElementType()->GetChemicalSymbol();

        // second layer
        storage::Vector< std::string> subt;
        subt.AllocateMemory( itr->GetTargetAtom().GetBonds().GetSize());
        for
        (
            auto itrb( itr->GetTargetAtom().GetBonds().Begin()), itrb_end( itr->GetTargetAtom().GetBonds().End());
            itrb != itrb_end;
            ++itrb
        )
        {
          std::string d;
          d += '0' + itrb->GetBondType()->GetBondData( ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness);
          d += itrb->GetTargetAtom().GetElementType()->GetChemicalSymbol();
          subt.PushBack( d);
        }
        subt.Sort( std::less< std::string>());
        subs.PushBack( util::Join( "", subt));
      }
      subs.Sort( std::less< std::string>());
      return A.GetElementType()->GetChemicalSymbol() + util::Join( "", subs);
    }

    //! @brief compute atom environment three bonds from reference atom
    //! @param A atom of interest
    //! @return return string indicating atom environment
    std::string AtomEnvironment::MakeAtomEnvironmentStringThree( const AtomConformationalInterface &A)
    {
      // first layer
      storage::Vector< std::string> subs;
      subs.AllocateMemory( A.GetBonds().GetSize());
      for
      (
          auto itr( A.GetBonds().Begin()), itr_end( A.GetBonds().End());
          itr != itr_end;
          ++itr
      )
      {
        std::string c;
        c += '0' + itr->GetBondType()->GetBondData( ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness);
        c += itr->GetTargetAtom().GetElementType()->GetChemicalSymbol();

        // second layer
        storage::Vector< std::string> subt;
        subt.AllocateMemory( itr->GetTargetAtom().GetBonds().GetSize());
        for
        (
            auto itrb( itr->GetTargetAtom().GetBonds().Begin()), itrb_end( itr->GetTargetAtom().GetBonds().End());
            itrb != itrb_end;
            ++itrb
        )
        {
          std::string d;
          d += '0' + itrb->GetBondType()->GetBondData( ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness);
          d += itrb->GetTargetAtom().GetElementType()->GetChemicalSymbol();

          // add chemical symbol from layer 2
          subt.PushBack( d);

          // third layer
          storage::Vector< std::string> subu;
          subu.AllocateMemory( itrb->GetTargetAtom().GetBonds().GetSize());
          for
          (
              auto itrc( itrb->GetTargetAtom().GetBonds().Begin()), itrc_end( itrb->GetTargetAtom().GetBonds().End());
              itrc != itrc_end;
              ++itrc
          )
          {
            std::string e;
            e += '0' + itrc->GetBondType()->GetBondData( ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness);
            e += itrc->GetTargetAtom().GetElementType()->GetChemicalSymbol();

            // add chemical symbol from layer 3
            subu.PushBack( e);
          }
          subu.Sort( std::less< std::string>());
          subt.PushBack( util::Join( "", subu));
        }
        subt.Sort( std::less< std::string>());
        subs.PushBack( util::Join( "", subt));
      }
      subs.Sort( std::less< std::string>());
      return A.GetElementType()->GetChemicalSymbol() + util::Join( "", subs);
    }

    //! @brief compute atom environment four bonds from reference atom
    //! @param A atom of interest
    //! @return return string indicating atom environment
    std::string AtomEnvironment::MakeAtomEnvironmentStringFour( const AtomConformationalInterface &A)
    {
      // first layer
      storage::Vector< std::string> subs;
      subs.AllocateMemory( A.GetBonds().GetSize());
      for
      (
          auto itr( A.GetBonds().Begin()), itr_end( A.GetBonds().End());
          itr != itr_end;
          ++itr
      )
      {
        std::string c;
        c += '0' + itr->GetBondType()->GetBondData( ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness);
        c += itr->GetTargetAtom().GetElementType()->GetChemicalSymbol();

        // second layer
        storage::Vector< std::string> subt;
        subt.AllocateMemory( itr->GetTargetAtom().GetBonds().GetSize());
        for
        (
            auto itrb( itr->GetTargetAtom().GetBonds().Begin()), itrb_end( itr->GetTargetAtom().GetBonds().End());
            itrb != itrb_end;
            ++itrb
        )
        {
          std::string d;
          d += '0' + itrb->GetBondType()->GetBondData( ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness);
          d += itrb->GetTargetAtom().GetElementType()->GetChemicalSymbol();

          // add chemical symbol from layer 2
          subt.PushBack( d);

          // third layer
          storage::Vector< std::string> subu;
          subu.AllocateMemory( itrb->GetTargetAtom().GetBonds().GetSize());
          for
          (
              auto itrc( itrb->GetTargetAtom().GetBonds().Begin()), itrc_end( itrb->GetTargetAtom().GetBonds().End());
              itrc != itrc_end;
              ++itrc
          )
          {
            std::string e;
            e += '0' + itrc->GetBondType()->GetBondData( ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness);
            e += itrc->GetTargetAtom().GetElementType()->GetChemicalSymbol();

            // add chemical symbol from layer 3
            subu.PushBack( e);

            // fourth layer
            storage::Vector< std::string> subz;
            subz.AllocateMemory( itrc->GetTargetAtom().GetBonds().GetSize());
            for
            (
                auto itrz( itrc->GetTargetAtom().GetBonds().Begin()), itrz_end( itrc->GetTargetAtom().GetBonds().End());
                itrz != itrz_end;
                ++itrz
            )
            {
              std::string e;
              e += '0' + itrz->GetBondType()->GetBondData( ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness);
              e += itrz->GetTargetAtom().GetElementType()->GetChemicalSymbol();

              // add chemical symbol from layer 4
              subz.PushBack( e);
            }
            subz.Sort( std::less< std::string>());
            subu.PushBack( util::Join( "", subz));
          }
          subu.Sort( std::less< std::string>());
          subt.PushBack( util::Join( "", subu));
        }
        subt.Sort( std::less< std::string>());
        subs.PushBack( util::Join( "", subt));
      }
      subs.Sort( std::less< std::string>());
      return A.GetElementType()->GetChemicalSymbol() + util::Join( "", subs);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write AtomEnvironment to std::ostream
    std::ostream &AtomEnvironment::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write base class
      io::Serialize::Write( m_AtomOfInterest, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_EnvironmentAtoms, OSTREAM, INDENT);

      // return
      return OSTREAM;
    } // Write

    //! read AtomEnvironment from io::IFStream
    std::istream &AtomEnvironment::Read( std::istream &ISTREAM)
    {
      // read base class
      io::Serialize::Read( m_AtomOfInterest, ISTREAM);
      io::Serialize::Read( m_EnvironmentAtoms, ISTREAM);

      // return
      return ISTREAM;
    } // Read

  } // namespace chemistry
} // namespace bcl
