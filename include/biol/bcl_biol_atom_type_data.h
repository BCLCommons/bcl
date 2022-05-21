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

#ifndef BCL_BIOL_ATOM_TYPE_DATA_H_
#define BCL_BIOL_ATOM_TYPE_DATA_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_biol_ss_types.h"
#include "chemistry/bcl_chemistry_element_types.h"
#include "coord/bcl_coord_cylinder_coordinates.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomTypeData
    //! @brief This is a low level helper class to store biological atom type properties
    //! @details This class acts as the storage class for the enumerator AtomTypes
    //!
    //! @see @link example_biol_atom_type_data.cpp @endlink
    //! @author woetzen, karakam
    //! @date 11/25/2007
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AtomTypeData :
      public util::ObjectInterface
    {
      friend class AtomTypes;

    private:

    //////////
    // data //
    //////////

      std::string                             m_AtomName;            //!< name of the atom defined by IUPAC
      chemistry::ElementType                  m_ElementType;         //!< element type of the atom
      bool                                    m_InBackBone;          //!< true if atom type is in backbone
      bool                                    m_InSideChain;         //!< true if atom type is in side chain
      coord::CylinderCoordinates              m_HelixCoordinates;    //!< Position on the cylinder for helix
      coord::CylinderCoordinates              m_StrandCoordinates;   //!< Position on the cylinder for strand
      mutable storage::Map< AtomType, double> m_BondLengths;         //!< Map of bond lengths to other atom types

      //! List of all possible connections
      //! If this atom type and the referenced atom type both exist in an AA, they are connected
      //! The sole exception is CD - N, which is only connected in proline, and so is not referenced by m_Connections
      mutable storage::Set< AtomType>         m_Connections;

      //! Connections (if any) which is always double bonded to this atom type, if they are present in both AAs
      //! AND neither of them is saturated with H
      //! This is either a member of m_Connections, or is undefined
      mutable storage::Set< AtomType>         m_DoubleBondConnections;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct undefined atom type
      AtomTypeData();

      //! @brief construct undefined atom type from given information
      //! @param ATOM_NAME name of the atom
      //! @param ELEMENT_TYPE element type of the atom
      //! @param IN_BACKBONE if atom type is in backbone
      //! @param IN_SIDE_CHAIN true if atom type is in side chain
      //! @param HELIX_CYLINDER_COORDINATES position on the cylinder for helix
      //! @param STRAND_CYLINDER_COORDINATES position on the cylinder for strand
      AtomTypeData
      (
        const std::string &ATOM_NAME,
        const chemistry::ElementType &ELEMENT_TYPE,
        const bool IN_BACKBONE,
        const bool IN_SIDE_CHAIN,
        const coord::CylinderCoordinates &HELIX_CYLINDER_COORDINATES,
        const coord::CylinderCoordinates &STRAND_CYLINDER_COORDINATES
      );

      //! @brief virtual copy constructor
      AtomTypeData *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the atom name
      //! @return the atom name by IUPAC
      const std::string &GetName() const
      {
        return m_AtomName;
      }

      //! @brief get the atom element type
      //! @return the atom element type
      const chemistry::ElementType &GetElementType() const
      {
        return m_ElementType;
      }

      //! @brief determines whether this atomtype is in backbone or not
      //! @return whether this atomtype is in backbone or not
      bool IsBackBone() const
      {
        return m_InBackBone;
      }

      //! @brief determines whether this atomtype is in side chain or not
      //! @return whether this atomtype is in side chain or not
      bool IsSideChain() const
      {
        return m_InSideChain;
      }

      //! @brief return Helix Cylinder Coordinates
      //! @return Helix Cylinder Coordinates
      const coord::CylinderCoordinates &GetHelixCoordinates() const
      {
        return m_HelixCoordinates;
      }

      //! @brief return Helix Cylinder Coordinates
      //! @return Helix Cylinder Coordinates
      const coord::CylinderCoordinates &GetStrandCoordinates() const
      {
        return m_StrandCoordinates;
      }

      //! @brief return Cylinder Coordinates for the specified SS_TYPE
      //! @param SS_TYPE SSType of interest
      //! @return Cylinder Coordinates for the specified SS_TYPE
      const coord::CylinderCoordinates &GetCoordinates( const SSType &SS_TYPE) const;

      //! @brief return atom types that this atom type will connect to, if present in the AA
      //! @return atom types that this atom type will connect to, if present in the AA
      const storage::Set< AtomType> &GetConnections() const;

      //! @brief return atom types that this atom type will connect to with a double bond, if present in the AA
      //! @return atom types that this atom type will connect to with a double bond, if present in the AA
      const storage::Set< AtomType> &GetDoubleBondConnections() const;

      //! @brief gets the bond length between this atom type and the passed atom type
      //! @param ATOM_TYPE atom type bound to
      //! @return the bond length between this atom type and the passed atom type
      double GetBondLength( const AtomType &ATOM_TYPE) const;

    private:

      //! @brief sets the bond length to another atom type
      //! @param ATOM_TYPE atom type bound to
      //! @param LENGTH average bond length
      void SetBondLength( const AtomType &ATOM_TYPE, const double LENGTH) const;

      //! @brief sets the connections that the atom always makes, if available
      //! @param CONNECTIONS set of atom types that this atom type always connects to
      void SetConnections( const storage::Set< AtomType> &CONNECTIONS) const;

      //! @brief sets the double bond connections that the atom always makes, if available
      //! @param CONNECTIONS set of atom types that this atom type always connects to with a double bond
      void SetDoubleBondConnections( const storage::Set< AtomType> &CONNECTIONS) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream
      //! @param INDENT number of indentations
      //! @return ostream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; //class AtomTypeData

  } // namespace biol
} // namespace bcl

#endif // BCL_BIOL_ATOM_TYPE_DATA_H_

