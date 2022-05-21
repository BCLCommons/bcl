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

// include the namespace headers
#include "chemistry/bcl_chemistry.h"
// include the header of this class
#include "chemistry/bcl_chemistry_molecule_environment.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_complete.h"
#include "chemistry/bcl_chemistry_atom_vector.h"
#include "chemistry/bcl_chemistry_fragment_complete.h" //!> input is a fragment complete
#include "io/bcl_io_ofstream.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_triplet.h"
#include "util/bcl_util_function_wrapper.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> MoleculeEnvironment::s_Elem
    (
      util::Enumerated< ConformationComparisonInterface>::AddInstance( new MoleculeEnvironment( AtomEnvironmentBender::e_Element))
    );
    //const util::SiPtr< const util::ObjectInterface> MoleculeEnvironment::s_ElemRC
    //(
    //  util::Enumerated< ConformationComparisonInterface>::AddInstance( new MoleculeEnvironment( AtomEnvironmentBender::e_ElemRC))
    //);
    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> MoleculeEnvironment::s_Atom
    (
      util::Enumerated< ConformationComparisonInterface>::AddInstance( new MoleculeEnvironment( AtomEnvironmentBender::e_Atom))
    );
    //const util::SiPtr< const util::ObjectInterface> MoleculeEnvironment::s_AtomRC
    //(
    //  util::Enumerated< ConformationComparisonInterface>::AddInstance( new MoleculeEnvironment( AtomEnvironmentBender::e_AtomRC))
    //);

    //! @brief atom type name as string
    //! @param ATOM_TYPE the atom type
    //! @return the string for the measure
    const std::string &MoleculeEnvironment::GetAtomTypeName( const t_AtomType &ATOM_TYPE)
    {
      return AtomEnvironmentBender::GetAtomTypeInfo( ATOM_TYPE).First();
    }

  //////////////////
  // Constructions//
  //////////////////

    //! constructor with Atom_type specification
    MoleculeEnvironment::MoleculeEnvironment( const t_AtomType &ATOM_TYPE) : m_AtomType( ATOM_TYPE)
    {
    }

    //! @brief constructor for building an atom environment from bond distance limit, index of atom, atom type, and fragment complete
    MoleculeEnvironment::MoleculeEnvironment( const t_AtomType &ATOM_TYPE, const ConformationInterface &FRAGMENT) :
      m_AtomType( ATOM_TYPE)
    {
      size_t index( 0);
      m_MoleculeEnv.AllocateMemory( FRAGMENT.GetSize());
      for
      (
          iterate::Generic< const AtomConformationalInterface> itr( FRAGMENT.GetAtomsIterator());
          itr.NotAtEnd();
          ++itr, ++index
      )
      {
        if( itr->GetElementType() != GetElementTypes().e_Hydrogen)
        {
          m_MoleculeEnv.PushBack( AtomEnvironmentBender( index, m_AtomType.GetEnum(), FRAGMENT));
        }
      }
      m_MoleculeEnv.Sort( std::less< AtomEnvironmentBender>());
    }

    MoleculeEnvironment::MoleculeEnvironment( const MoleculeEnvironment &OTHER) :
      m_MoleculeEnv( OTHER.m_MoleculeEnv),
      m_AtomType( OTHER.m_AtomType)
    {
    }

    //! copy constructor
    MoleculeEnvironment *MoleculeEnvironment::Clone() const
    {
      return new MoleculeEnvironment( *this);
    }

    //! @brief compares two atom environments
    bool MoleculeEnvironment::operator ==( const MoleculeEnvironment &MOLECULE) const
    {
      return m_MoleculeEnv == MOLECULE.m_MoleculeEnv;
    }

    //! @brief compares two atom environments
    bool MoleculeEnvironment::operator !=( const MoleculeEnvironment &MOLECULE) const
    {
      return m_MoleculeEnv != MOLECULE.m_MoleculeEnv;
    }

    //! @brief compute the tanimoto score between molecule environments of THIS and OTHER
    double MoleculeEnvironment::operator()
    (
      const ConformationInterface& THIS,
      const ConformationInterface& OTHER
    ) const
    {
      MoleculeEnvironment this_molecule( m_AtomType, FragmentComplete( THIS));
      MoleculeEnvironment other_molecule( m_AtomType, FragmentComplete( OTHER));
      return this_molecule.BuserScore( other_molecule);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeEnvironment::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( AtomEnvironmentBender::GetAtomTypeInfo( m_AtomType).Second());
      return serializer;
    }
    //////////////////
    // data access ///
    //////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MoleculeEnvironment::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the atom type that this atom enviroment calculator calculates
    const MoleculeEnvironment::t_AtomType &MoleculeEnvironment::GetAtomType() const
    {
      return m_AtomType;
    }

    //! @brief returns the data label
    //! @return data label as string
    const std::string &MoleculeEnvironment::GetAlias() const
    {
      return ( AtomEnvironmentBender::GetAtomTypeInfo( m_AtomType).First());
    }

    //! @brief returns the atom environment
    const MoleculeEnvironment::t_MoleculeEnv &MoleculeEnvironment::GetMoleculeEnvironment() const
    {
      return m_MoleculeEnv;
    }

    // @brief converts vector of sorted objects into a map of objects and their counts
    storage::Map< AtomEnvironmentBender, size_t> MoleculeEnvironment::ConvertVectorToMap() const
    {
      storage::Map< AtomEnvironmentBender, size_t> map;
      for( t_MoleculeEnv::const_iterator iter = m_MoleculeEnv.Begin(), end_iter = m_MoleculeEnv.End(); iter != end_iter; ++iter)
      {
        AddToMap( *iter, map);
      }
      return map;
    }

    //! @brief adds the Key into MAP or increment its count( the value assiated with that key)
    void MoleculeEnvironment::AddToMap( const AtomEnvironmentBender &KEY, storage::Map< AtomEnvironmentBender, size_t> &MAP)
    {
      std::pair< typename storage::Map< AtomEnvironmentBender, size_t>::iterator, bool> insert_pair
      (
        MAP.Insert( storage::Pair< AtomEnvironmentBender, size_t>( KEY, size_t( 1)))
      );
      if( !( insert_pair.second))
      {
        insert_pair.first->second += 1;
      }
    }

    // @brief Compute the TanimotoScore between this and OTHER
    double MoleculeEnvironment::TanimotoScore ( const MoleculeEnvironment &OTHER) const
    {
      double sum_a( GetMoleculeEnvironment().GetSize());
      double sum_b( OTHER.GetMoleculeEnvironment().GetSize());
      storage::Map< AtomEnvironmentBender, size_t> lhs( ConvertVectorToMap());
      t_MoleculeEnv key_lhs( lhs.GetKeysAsVector());
      storage::Map< AtomEnvironmentBender, size_t> rhs( OTHER.ConvertVectorToMap());
      t_MoleculeEnv key_rhs( rhs.GetKeysAsVector());
      std::vector< AtomEnvironmentBender> shared_atom;
      std::set_intersection
      (
        key_lhs.Begin(), key_lhs.End(), key_rhs.Begin(), key_rhs.End(), std::back_inserter( shared_atom)
      );
      double sum_ab( 0);
      for( std::vector< AtomEnvironmentBender>::const_iterator iter = shared_atom.begin(), end_iter = shared_atom.end(); iter != end_iter; ++iter)
      {
        size_t a( lhs.Find( *iter)->second);
        size_t b( rhs.Find( *iter)->second);
        sum_ab += std::min( a, b);
      }
      if( !sum_ab)
      {
        return 0;                                          //!> If sum_ab == 0
      }
      return double( sum_ab / ( sum_a + sum_b - sum_ab));
    }

    // @brief Computes the Buser score between this and OTHER
    double MoleculeEnvironment::BuserScore ( const MoleculeEnvironment &OTHER) const
    {
      //! A: # attributes in commom between the two molecule
      //! B: # attributes present in the lsh but not the rhs
      //! C: # attributes present in the rhs but not the lhs
      //! D: # average length of rhs and lhs
      double a( 0.0);
      double b( GetMoleculeEnvironment().GetSize());
      double c( OTHER.GetMoleculeEnvironment().GetSize());
      double d( ( b + c) / 2.0);
      storage::Map< AtomEnvironmentBender, size_t> lhs( ConvertVectorToMap());
      storage::Map< AtomEnvironmentBender, size_t> rhs( OTHER.ConvertVectorToMap());
      auto l_iter( lhs.Begin());
      auto r_iter( rhs.Begin());
      while( l_iter != lhs.End() && r_iter != rhs.End())
      {
        if( l_iter->first < r_iter->first)
        {
          ++l_iter;
        }
        else if( r_iter->first < l_iter->first)
        {
          ++r_iter;
        }
        else
        {
          a += std::min( l_iter->second, r_iter->second);
          ++r_iter;
          ++l_iter;
        }
      }
      b -= a;
      c -= a;
      return ( math::Sqrt( a * d) + a) / ( math::Sqrt( a * d) + a + b + c);
    }
  } // namespace chemistry
} // namespace bcl

