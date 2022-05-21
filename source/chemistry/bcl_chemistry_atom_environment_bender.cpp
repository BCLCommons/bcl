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
#include "chemistry/bcl_chemistry_atom_environment_bender.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_complete.h"
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "chemistry/bcl_chemistry_fragment_complete.h" //!> input is a fragment complete
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_bond_info.h"
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    // separate, anonymous namespace to prevent exporting symbols
    namespace
    {
      // add each of the possible instances to the enumerated instances
      util::ObjectInterface *AddInstances()
      {
        // keep a pointer to the last created instance
        util::ObjectInterface *last_instance( NULL);
        for( size_t atom_type( 0); atom_type < AtomEnvironmentBender::s_NumberAtomTypes; ++atom_type)
        {
          last_instance =
              util::Enumerated< util::SerializableInterface>::AddInstance
              (
                new AtomEnvironmentBender( static_cast< AtomEnvironmentBender::Atom_type>( atom_type))
              );
        }
        return last_instance;
      }
    }

    //single instance of that class
    const util::SiPtr< const util::ObjectInterface> AtomEnvironmentBender::s_Instances( AddInstances());

    //! @brief Get the name, description, and function of the given atom type
    //! @param ATOM_TYPE the atom type of interest
    //! @return the short name or abbreviation for the class
    const storage::Triplet< std::string, std::string, AtomEnvironmentBender::t_functions>
    &AtomEnvironmentBender::GetAtomTypeInfo( const Atom_type &ATOM_TYPE)
    {
      typedef storage::Triplet< std::string, std::string, AtomEnvironmentBender::t_functions> t_Info;
      AtomEnvironmentBender::t_functions  elem_function( &AtomEnvironmentBender::ElementHash, &AtomEnvironmentBender::UnElemHash );
      //AtomEnvironmentBender::t_functions  elem_functionRC( &AtomEnvironmentBender::ElemRCHash, &AtomEnvironmentBender::UnElemRCHash );
      AtomEnvironmentBender::t_functions  atom_function( &AtomEnvironmentBender::AtomHash, &AtomEnvironmentBender::UnAtomHash );
      //AtomEnvironmentBender::t_functions  atom_functionRC( &AtomEnvironmentBender::AtomRCHash, &AtomEnvironmentBender::UnAtomRCHash );
      static const t_Info s_info[ s_NumberAtomTypes] =
      {
          t_Info( "Element", "Atomic Number and bond orders", elem_function),
          t_Info( "Atom", "Atomic Number and bond orders", atom_function),
      };
      return s_info[ ATOM_TYPE];
    };

    //! @brief atom type name as string
    //! @param ATOM_TYPE the atom type
    //! @return the string for the measure
    const std::string &AtomEnvironmentBender::GetAtomTypeName( const Atom_type &ATOM_TYPE)
    {
      return GetAtomTypeInfo( ATOM_TYPE).First();
    }

  //////////////////
  // Constructions//
  //////////////////

    //! default constructor
    AtomEnvironmentBender::AtomEnvironmentBender() : m_BondNumLimit( 1)
    {
    }

    //! constructor with Atom_type specification
    AtomEnvironmentBender::AtomEnvironmentBender( const Atom_type &ATOM_TYPE) : m_BondNumLimit( 1), m_AtomType( ATOM_TYPE)
    {
    }

    //! constructor from AE string representation
    AtomEnvironmentBender::AtomEnvironmentBender( const Atom_type &ATOM_TYPE, const std::string &STRING) :
        m_BondNumLimit( 1),
        m_AtomType( ATOM_TYPE),
        m_AtomEnvironment( StringToAE( STRING, ATOM_TYPE))
    {
    }

    //! @brief constructor for building an atom environment from bond distance limit, index of atom, atom type, and fragment complete
    AtomEnvironmentBender::AtomEnvironmentBender
    (
      size_t ATOM_INDEX, const Atom_type &ATOM_TYPE, const ConformationInterface &FRAGMENT
    ) :
      m_AtomIndex( ATOM_INDEX), m_BondNumLimit( 1), m_AtomType( ATOM_TYPE),
      m_AtomEnvironment( m_BondNumLimit + 1, storage::Map< size_t, size_t>())
    {
      t_output s_output( AtomEnvironmentBender::GetEncodedAtomEnvironment( FRAGMENT));
      storage::Vector< std::string> encoded_atom_environment = s_output.First();
      storage::Vector< size_t> bond_orders = s_output.Second();
      DecodeAtomEnvironment
      (
        encoded_atom_environment,
        HashMolecule( FRAGMENT, bond_orders),
        m_BondNumLimit
      );
    }

    //! copy constructor
    AtomEnvironmentBender *AtomEnvironmentBender::Clone() const
    {
      return new AtomEnvironmentBender( *this);
    }

    //! @brief compares two atom environments
    bool AtomEnvironmentBender::operator ==( const AtomEnvironmentBender &ATOM) const
    {
      return ( m_AtomEnvironment == ATOM.m_AtomEnvironment);
    }

    //! @brief compares two atom environments
    bool AtomEnvironmentBender::operator <( const AtomEnvironmentBender &ATOM) const
    {
      return ( m_AtomEnvironment < ATOM.m_AtomEnvironment);
    }

    //! @brief compares two atom environments
    bool AtomEnvironmentBender::operator !=( const AtomEnvironmentBender &ATOM) const
    {
      return ( m_AtomEnvironment != ATOM.m_AtomEnvironment);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AtomEnvironmentBender::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( GetAtomTypeInfo( m_AtomType).Second());
      return serializer;
    }
    //////////////////
    // data access ///
    //////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AtomEnvironmentBender::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief access the atom of interest
    std::size_t AtomEnvironmentBender::GetAtomOfInterestIndex() const
    {
      return m_AtomIndex;
    }

    //! @brief get the atom type that this atom enviroment calculator calculates
    const AtomEnvironmentBender::Atom_type &AtomEnvironmentBender::GetAtomType() const
    {
      return m_AtomType;
    }

    //! @brief returns the data label
    //! @return data label as string
    const std::string &AtomEnvironmentBender::GetAlias() const
    {
      return ( GetAtomTypeInfo( m_AtomType).First());
    }

    //! @brief returns the atom environment
    const AtomEnvironmentBender::t_AtomEnvironment &AtomEnvironmentBender::GetAtomEnvironment() const
    {
      return m_AtomEnvironment;
    }

    //! @brief read the fragment complete of the molecule and and the index of the atom of interest,
    //!        and then provide encoded version of the atom environment.
    //! @param FRAGMENT the fragment complete representing the molecule
    //! @param BOND_LIMIT the limit number of bond away from the atom of interest
    AtomEnvironmentBender::t_output AtomEnvironmentBender::GetEncodedAtomEnvironment
    (
      const ConformationInterface &FRAGMENT
    )
    {
      //!> get the number of atoms in the molecule
      size_t moleculeSize = FRAGMENT.GetSize();

      //!> initializes encoded environment with all '-1' values
      storage::Vector< std::string> encoded_atom_environment( m_BondNumLimit + 1, std::string( moleculeSize, '0'));
      storage::Vector< size_t> bond_orders( moleculeSize, size_t( 0));

      // the atom itself would be the only '1' in the first string
      // this function will recursively do the rest of work
      encoded_atom_environment( 0)[ m_AtomIndex] = '1';
      GetEncodedAtomEnvironmentHelper( FRAGMENT, m_AtomIndex, m_BondNumLimit, 1, encoded_atom_environment, bond_orders);
      AtomEnvironmentBender::t_output s_output( encoded_atom_environment, bond_orders);
      return s_output;
    }

    //! @brief recursively find the neighbor atoms in the atom environment of the atom of interest and then update
    //! @param ENCODED_ATOM_ENVIRONMENT the encoded atom enviroment, each string represents a phere,
    //          each character represent absence/presece of the corresponding atom
    void AtomEnvironmentBender::GetEncodedAtomEnvironmentHelper
    (
      const ConformationInterface &FRAGMENT, size_t ATOM_INDEX, size_t BOND_LIMIT,
      size_t BOND_NUM, storage::Vector< std::string> &ENCODED_ATOM_ENVIRONMENT,
      storage::Vector< size_t> &BOND_ORDERS
    )
    {
      if( BOND_NUM <= BOND_LIMIT)
      {
        iterate::Generic< const AtomConformationalInterface> atom( FRAGMENT.GetAtomsIterator());
        atom.GotoPosition( ATOM_INDEX);
        const storage::Vector< BondConformational> &bond_vector( atom->GetBonds());
        for
        (
            storage::Vector< BondConformational>::const_iterator itr = bond_vector.Begin(),
            itr_end = bond_vector.End(); itr != itr_end;
            ++itr
        )
        {
          if( itr->IsDefined() && ( itr->GetTargetAtom().GetElementType() != GetElementTypes().e_Hydrogen))
          {
            size_t bond_order( 4); //preset the bondorder to aromatic bond
            // if not aromatic bond, change to normal bond order scheme (1, 2, 3)
            if( !( itr->GetBondType()->GetBondData( ConfigurationalBondTypeData::e_IsAromatic)))
            {
              bond_order = itr->GetBondType()->GetBondData( ConfigurationalBondTypeData::e_BondOrder);
            }
            //size_t bond_order( itr->GetBondType()->GetBondData( ConfigurationalBondTypeData::e_BondOrder) );
            size_t target_atom_index = FRAGMENT.GetAtomIndex( itr->GetTargetAtom());
            UpdateEncodedAtomEnvironment( BOND_NUM, target_atom_index, ENCODED_ATOM_ENVIRONMENT, bond_order, BOND_ORDERS);
            GetEncodedAtomEnvironmentHelper
            (
              FRAGMENT, target_atom_index, BOND_LIMIT, BOND_NUM + 1, ENCODED_ATOM_ENVIRONMENT, BOND_ORDERS
            );
          }
        }
      }
    }
    //! @brief returns the corresponding character representing the bond order
    //! @param BOND_ORDER a integer representing the bond order
    std::string AtomEnvironmentBender::GetBondOrder( size_t BOND_ORDER)
    {
      switch( BOND_ORDER)
      {
        case 1: return "-";
        case 2: return "=";
        case 3: return "#";
        case 4: return "~";
        default: return "_";
      }
    }

    //! @brief updates the encoded atom environment
    //! @param TARGET_ATOM_INDEX the index of the target atom, which bond to the atom of interest
    //! @param ENCODED_ATOM_ENVIRONMENT a vector of string that represents the neighbor atoms in different spheres
    void AtomEnvironmentBender::UpdateEncodedAtomEnvironment
    (
      size_t BOND_NUM, size_t TARGET_ATOM_INDEX,
      storage::Vector< std::string> &ENCODED_ATOM_ENVIRONMENT,
      size_t BOND_ORDER, storage::Vector< size_t> &BOND_ORDERS
    )
    {
      if( ( ENCODED_ATOM_ENVIRONMENT( BOND_NUM)[ TARGET_ATOM_INDEX] != '0') && ( BOND_ORDERS( TARGET_ATOM_INDEX) < BOND_ORDER))
      {
        BOND_ORDERS( TARGET_ATOM_INDEX) = BOND_ORDER;
        return;
      }
      bool update = true;
      for( size_t i = 0; i <= BOND_NUM; i++)
      {
        if( ENCODED_ATOM_ENVIRONMENT( i)[ TARGET_ATOM_INDEX] != '0')
        {
          update = false;
          break;
        }
      }
      if( update)
      {
        ENCODED_ATOM_ENVIRONMENT( BOND_NUM)[ TARGET_ATOM_INDEX] = '1';
        BOND_ORDERS( TARGET_ATOM_INDEX) = BOND_ORDER;
      }
    }

    //! @brief converts from the encoded atom environment to the normal representation of the atom environment
    //! @param ENCODED_ATOM_ENVIRONMENT encoded atom environment in term of a vector of string of 0 and 1
    //! @param HASHED_MOLECULE is vector of hashed atoms of the molecule
    void AtomEnvironmentBender::DecodeAtomEnvironment
    (
      const storage::Vector< std::string> &ENCODED_ATOM_ENVIRONMENT,
      const storage::Vector< size_t> &HASHED_MOLECULE,
      size_t BOND_LIMIT
    )
    {
      size_t size( HASHED_MOLECULE.GetSize());
      for( size_t bond_num( 0); bond_num <= BOND_LIMIT; bond_num++)
      {
        for( size_t atom_index = 0; atom_index <= size; atom_index++)
        {
          if( ENCODED_ATOM_ENVIRONMENT( bond_num)[ atom_index] == '1')
          {
            this->AddAtom( HASHED_MOLECULE( atom_index), bond_num);   //!> adds the hashed atom into m_AtomEnvironment
          }
        }
      }
    }

    //! @brief adds the Key into MAP or increment its count( the value assiated with that key)
    template< typename K>
    void AtomEnvironmentBender::AddToMap( const K &KEY, storage::Map< K, size_t> &MAP)
    {
      std::pair< typename storage::Map< K, size_t>::iterator, bool> insert_pair
      (
        MAP.Insert( storage::Pair< K, size_t>( KEY, size_t( 1)))
      );
      if( !( insert_pair.second))
      {
        insert_pair.first->second += 1;
      }
    }

    //! @brief adds the hashed atom into m_AtomEnvironment
    //! @param HASHED_ATOM
    //! @param BOND_DISTANCE
    void AtomEnvironmentBender::AddAtom( size_t HASHED_ATOM, size_t BOND_DISTANCE)
    {
      AddToMap< size_t>( HASHED_ATOM, m_AtomEnvironment( BOND_DISTANCE));
    }

    ////////////////////////////////////////
    // helper functions for hash functions /
    ////////////////////////////////////////

    //! @brief returns the chemical symbol of an element from its atomic number ATOM_NUMBER
    const std::string &AtomEnvironmentBender::GetAtomicSymbol( size_t ATOM_NUMBER, bool ATOM_FLAG)
    {
      if( ATOM_FLAG) // if this fingerprint uses atom type
      {
        AtomType atom( ATOM_NUMBER);
        return atom.GetName();
      }
      else // if this fingerprint uses element type
      {
        ElementType element( ATOM_NUMBER);
        return element->GetChemicalSymbol();
      }
    }

    //! @brief return the enum index of the atom type
    const size_t AtomEnvironmentBender::GetAtomTypeIndex( const AtomConformationalInterface &ATOM)
    {
      return ATOM.GetAtomType()->GetIndex();
    }

    //! @brief return the enum index of the element type
    //! @param ATOM: the index of the atom of interest
    const size_t AtomEnvironmentBender::GetAtomicNumber( const AtomConformationalInterface &ATOM)
    {
      return ATOM.GetElementType().GetIndex();
    }

    ///////////////////////////////
    // Hash and unhash functions /
    //////////////////////////////

    //! @brief choose the hash function based on the choice of atom type, then hash every atoms of a fragment complete
    //! @return a number which is a compressed version of the atom type
    storage::Vector< size_t> AtomEnvironmentBender::HashMolecule
    (
      const ConformationInterface &FRAGMENT,
      const storage::Vector< size_t> &BOND_ORDERS
    ) const
    {
      size_t size = FRAGMENT.GetNumberAtoms();
      storage::Vector< size_t> hashed_atoms( size, size_t( 0));
      size_t index( 0);
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr( FRAGMENT.GetAtomsIterator());
        itr.NotAtEnd();
        ++itr, ++index
      )
      {
        function::MemberUnaryConst< AtomEnvironmentBender, const AtomEnvironmentBender::t_input &, size_t>
        function( GetAtomTypeInfo( m_AtomType).Third().First());
        t_input s_input( *itr, BOND_ORDERS( index));
        hashed_atoms( index) = function( *this, s_input); //!> add the hashed atom into the corresponding array element
      }
      return hashed_atoms;
    }

    //! @brief unhashes the atom environment into its string representation
    std::string AtomEnvironmentBender::UnHash() const
    {
      function::MemberUnaryConst< AtomEnvironmentBender, const size_t &, std::string> function( GetAtomTypeInfo( m_AtomType).Third().Second());
      std::ostringstream oss;
      for( size_t i = 0; i <= m_BondNumLimit; ++i)
      {
        oss << "[";
        for
        (
            storage::Map< size_t, size_t>::const_iterator itr = m_AtomEnvironment( i).Begin(),
            itr_end = m_AtomEnvironment( i).End();
            itr != itr_end; ++itr
        )
        {
          for( size_t j = 0; j < itr->second; ++j)
          {
            oss << function( *this, itr->first) + " ";
          }
        }
        oss << "]";
      }
      return oss.str();
    }

    //! @brief hashes the Element atom type info of the INPUT
    size_t AtomEnvironmentBender::ElementHash( const t_input &INPUT) const
    {
      return ( INPUT.Second() << 7) | GetAtomicNumber( INPUT.First());
    }

    //! @brief converts the hashed atom into its string representation
    std::string AtomEnvironmentBender::UnElemHash( const size_t &HASHEDATOM) const
    { // need to check the number of bits?
      std::string symbol( GetAtomicSymbol( HASHEDATOM & 127, false));
      std::string bond_order( GetBondOrder( ( HASHEDATOM >> 7) & 7));
      return ( bond_order + symbol);
    }

    //! @brief hashes the Element atom type info of the INPUT
    size_t AtomEnvironmentBender::AtomHash( const t_input &INPUT) const
    {
      return ( INPUT.Second() << 8) | GetAtomTypeIndex( INPUT.First());
    }

    //! @brief converts the hashed atom into its string representation
    std::string AtomEnvironmentBender::UnAtomHash( const size_t &HASHEDATOM) const
    { // need to check the number of bits?
      std::string symbol( GetAtomicSymbol( HASHEDATOM & 255, true));
      std::string bond_order( GetBondOrder( ( HASHEDATOM >> 8) & 7));
      return ( bond_order + symbol);
    }

    //! @brief converts the hashed string back to the AE
    AtomEnvironmentBender::t_AtomEnvironment AtomEnvironmentBender::StringToAE( const std::string &STRING, const Atom_type ATOM_TYPE)
    {
      t_AtomEnvironment AE;
      storage::Vector< std::string> tokens( util::SplitString( STRING, "[]"));
      size_t count( 0);
      for( storage::Vector< std::string>::const_iterator iter = tokens.Begin(); iter != tokens.End(); ++iter, count++)
      {
        storage::Vector< std::string> atoms( util::SplitString( *iter, " "));
        storage::Map< size_t, size_t> layer;
        for( storage::Vector< std::string>::const_iterator itr = atoms.Begin(); itr != atoms.End(); ++itr)
        {
          size_t hashed_atom;
          if( ATOM_TYPE == e_Element)
          {
            hashed_atom = ElementStringToAE(*itr);
          }
          else
          {
            hashed_atom = AtomStringToAE( *itr);
          }
          AddToMap< size_t>( hashed_atom, layer);
        }
        AE.PushBack( layer);
      }
      return AE;
    } // namespace chemistry

    //! @brief converts the hashed string back to the Element AE
    size_t AtomEnvironmentBender::ElementStringToAE( const std::string &STRING)
    {
      return CharToBond( STRING[ 0]) << 7 | StringToElement( STRING.substr( 1, std::string::npos));
    }

    //! @brief converts the hashed string back to the Element AE
    size_t AtomEnvironmentBender::AtomStringToAE( const std::string &STRING)
    {
      return CharToBond( STRING[ 0]) << 8 | StringToAtom( STRING.substr( 1, std::string::npos));
    }

    //! @brief converts the string back to bond info
    size_t AtomEnvironmentBender::CharToBond( const char &CHAR)
    {
      switch( CHAR)
      {
        case '-': return 1;
        case '=': return 2;
        case '#': return 3;
        case '~': return 4;
        default: return 0;
      }
      return 0;
    }

    //! @brief converts the string back to the Element info
    size_t AtomEnvironmentBender::StringToElement( const std::string &STRING)
    {
      return GetElementTypes().ElementTypeLookup( STRING).GetIndex();
    }

    //! @brief converts the string back to the Atom type info
    size_t AtomEnvironmentBender::StringToAtom( const std::string &STRING)
    {
      return AtomType( STRING).GetIndex();
    }
  } // namespace chemistry
} // namespace bcl
