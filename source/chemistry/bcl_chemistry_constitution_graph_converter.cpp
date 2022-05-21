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
#include "chemistry/bcl_chemistry_constitution_graph_converter.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  ///////////
  // Enums //
  ///////////

    //! @brief Data as string
    //! @param DATA the data whose name is desired
    //! @return the name as string
    const std::string &ConstitutionGraphConverter::GetAtomComparisonType( const AtomComparisonType &DATA)
    {
      static const std::string s_Names[ 1 + s_NumberAtomComparisonTypes] =
      {
        "Identity",
        "ElementType",
        "AtomType",
        GetStaticClassName< AtomComparisonType>()
      };
      return s_Names[ DATA];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor, compares by atom type and e_BondOrderAmideOrAromaticWithRingness
    ConstitutionGraphConverter::ConstitutionGraphConverter() :
      m_AtomRepresentation( e_AtomType),
      m_BondRepresentation( ConstitutionalBondTypeData::e_BondOrderAmideOrAromaticWithRingness)
    {
    }

    //! @brief constructor from coloring schemes
    //! @param ATOM_COMPARISON_TYPE means by which to color vertices of a molecule
    //! @param BOND_TYPE_INFO the desired information use for the bonds
    ConstitutionGraphConverter::ConstitutionGraphConverter
    (
      const AtomComparisonType &ATOM_COMPARISON_TYPE,
      const ConstitutionalBondTypeData::Data &BOND_TYPE_INFO
    ) :
      m_AtomRepresentation( ATOM_COMPARISON_TYPE),
      m_BondRepresentation( BOND_TYPE_INFO)
    {
    }

    //! @brief virtual copy constructor
    ConstitutionGraphConverter *ConstitutionGraphConverter::Clone() const
    {
      return new ConstitutionGraphConverter( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ConstitutionGraphConverter::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief Create a graph of a constitution
    //! @param CONSTITUTION a constitution
    //! @return The constitution converted into a graph with the given coloring schemes
    graph::ConstGraph< size_t, size_t>
      ConstitutionGraphConverter::operator()( const ConstitutionInterface &CONSTITUTION) const
    {
      storage::Vector< size_t> vertices;
      vertices.AllocateMemory( CONSTITUTION.GetNumberAtoms());

      // get data for each atom
      for
      (
        iterate::Generic< const AtomConstitutionalInterface> itr( CONSTITUTION.GetAtomsIterator());
        itr.NotAtEnd();
        ++itr
      )
      {
        vertices.PushBack( GetAtomData( *itr));
      }

      // construct the graph and return it
      return graph::ConstGraph< size_t, size_t>
             (
               vertices,
               CONSTITUTION.GetAdjacencyList( m_BondRepresentation),
               GetConstitutionalBondTypes().e_Undefined
             );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ConstitutionGraphConverter::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_AtomRepresentation, ISTREAM);
      io::Serialize::Read( m_BondRepresentation, ISTREAM);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &ConstitutionGraphConverter::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_AtomRepresentation, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_BondRepresentation, OSTREAM);
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief helper function to get AtomComparisonType from an atom
    //! @param ATOM atom to retrieve data from
    //! @param TYPE atom comparison type
    //! @return data for the atome atom
    size_t ConstitutionGraphConverter::GetAtomData( const AtomConstitutionalInterface &ATOM) const
    {
      switch( m_AtomRepresentation)
      {
        case e_Identity:
          return size_t( 1);
        case e_ElementType:
          return ATOM.GetElementType().GetIndex();
        case e_AtomType:
          return ATOM.GetAtomType().GetIndex();
        case s_NumberAtomComparisonTypes:
        default:
          return util::GetUndefined< size_t>();
      }
      return util::GetUndefined< size_t>();
    }

  } // namespace chemistry
} // namespace bcl

