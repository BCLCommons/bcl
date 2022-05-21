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
#include "biol/bcl_biol_atom_type_data.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom_types.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AtomTypeData::s_Instance( GetObjectInstances().AddInstance( new AtomTypeData()));

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct undefined atom type
    AtomTypeData::AtomTypeData() :
      m_AtomName( ""),
      m_ElementType( util::GetUndefined< chemistry::ElementType>()),
      m_InBackBone( false),
      m_InSideChain( false),
      m_HelixCoordinates( util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>()),
      m_StrandCoordinates( util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>()),
      m_BondLengths(),
      m_Connections()
    {
    }

    //! @brief construct undefined atom type from given information
    //! @param ATOM_NAME name of the atom
    //! @param ELEMENT_TYPE element type of the atom
    //! @param IN_BACKBONE if atom type is in backbone
    //! @param IN_SIDE_CHAIN true if atom type is in side chain
    //! @param HELIX_CYLINDER_COORDINATES position on the cylinder for helix
    //! @param STRAND_CYLINDER_COORDINATES position on the cylinder for strand
    AtomTypeData::AtomTypeData
    (
      const std::string &ATOM_NAME,
      const chemistry::ElementType &ELEMENT_TYPE,
      const bool IN_BACKBONE,
      const bool IN_SIDE_CHAIN,
      const coord::CylinderCoordinates &HELIX_CYLINDER_COORDINATES,
      const coord::CylinderCoordinates &STRAND_CYLINDER_COORDINATES
    ) :
      m_AtomName( ATOM_NAME),
      m_ElementType( ELEMENT_TYPE),
      m_InBackBone( IN_BACKBONE),
      m_InSideChain( IN_SIDE_CHAIN),
      m_HelixCoordinates( HELIX_CYLINDER_COORDINATES),
      m_StrandCoordinates( STRAND_CYLINDER_COORDINATES),
      m_BondLengths(),
      m_Connections()
    {
    }

    //! @brief virtual copy constructor
    AtomTypeData *AtomTypeData::Clone() const
    {
      return new AtomTypeData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AtomTypeData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return Cylinder Coordinates for the specified SS_TYPE
    //! @param SS_TYPE SSType of interest
    //! @return Cylinder Coordinates for the specified SS_TYPE
    const coord::CylinderCoordinates &AtomTypeData::GetCoordinates( const SSType &SS_TYPE) const
    {
      static coord::CylinderCoordinates s_coordinates;

      if( SS_TYPE == GetSSTypes().HELIX)
      {
        return m_HelixCoordinates;
      }
      if( SS_TYPE == GetSSTypes().STRAND)
      {
        return m_StrandCoordinates;
      }
      BCL_Exit( "No cylindrical coordinates are stored for the provided sstype " + util::Format()( SS_TYPE), -1);

      return s_coordinates;
    }

    //! @brief return Cylinder Coordinates for the specified SS_TYPE
    //! @param SS_TYPE SSType of interest
    //! @return Cylinder Coordinates for the specified SS_TYPE
    const storage::Set< AtomType> &AtomTypeData::GetConnections() const
    {
      return m_Connections;
    }

    //! @brief return atom types that this atom type will connect to with a double bond, if present in the AA
    //! @return atom types that this atom type will connect to with a double bond, if present in the AA
    const storage::Set< AtomType> &AtomTypeData::GetDoubleBondConnections() const
    {
      return m_DoubleBondConnections;
    }

    //! @brief gets the bond length between this atom type and the passed atom type
    //! @param ATOM_TYPE atom type bound to
    //! @return the bond length between this atom type and the passed atom type
    double AtomTypeData::GetBondLength( const AtomType &ATOM_TYPE) const
    {
      // search the map for the atom type
      const storage::Map< AtomType, double>::const_iterator find_itr( m_BondLengths.Find( ATOM_TYPE));

      // return an undefined double if it wasn't found, otherwise return the length
      return find_itr == m_BondLengths.End() ? util::GetUndefined< double>() : find_itr->second;
    }

    //! @brief sets the bond length to another atom type
    //! @param ATOM_TYPE atom type bound to
    //! @param LENGTH average bond length
    void AtomTypeData::SetBondLength( const AtomType &ATOM_TYPE, const double LENGTH) const
    {
      m_BondLengths[ ATOM_TYPE] = LENGTH;
    }

    //! @brief sets the connections that the atom always makes, if available
    //! @param CONNECTIONS set of atom types that this atom type always connects to
    void AtomTypeData::SetConnections( const storage::Set< AtomType> &CONNECTIONS) const
    {
      m_Connections = CONNECTIONS;
    }

    //! @brief sets the double bond connections that the atom always makes, if available
    //! @param CONNECTIONS set of atom types that this atom type always makes double bonds to
    void AtomTypeData::SetDoubleBondConnections( const storage::Set< AtomType> &CONNECTIONS) const
    {
      m_DoubleBondConnections = CONNECTIONS;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AtomTypeData::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_AtomName, ISTREAM);
      io::Serialize::Read( m_ElementType, ISTREAM);
      io::Serialize::Read( m_InBackBone, ISTREAM);
      io::Serialize::Read( m_HelixCoordinates, ISTREAM);
      io::Serialize::Read( m_StrandCoordinates, ISTREAM);
      io::Serialize::Read( m_BondLengths, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @return ostream which was written to
    std::ostream &AtomTypeData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_AtomName, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ElementType, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_InBackBone, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_HelixCoordinates, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_StrandCoordinates, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_BondLengths, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
