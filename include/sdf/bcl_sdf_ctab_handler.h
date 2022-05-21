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

#ifndef BCL_SDF_CTAB_HANDLER_H_
#define BCL_SDF_CTAB_HANDLER_H_

// include the namespace header
#include "bcl_sdf.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_sdf_atom_info.h"
#include "bcl_sdf_bond_info.h"
#include "bcl_sdf_mdl_header.h"
#include "bcl_sdf_mdl_line_types.h"
#include "bcl_sdf_mdl_property.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CTabHandler
    //! @brief handles reading and writing connection table (CTab) formatted files.
    //!        CTab stands for connection table, which is a common connectivity specification used across all MDL file formats
    //!        A ctab entry contains a few main sections, namely the counts line, the atom block, the bond block, MDL 
    //!        properties, and an ending delimiter.  For details see the MDL ctfile format specification (in this directory).
    //!        Some quick details are:
    //!        Counts line: first two fields (up to 3 digits each) indicate number of atoms and bonds respectively
    //!          the counts line also specifies the format, V2000 or V3000, of the ctab format.  V2000 is the most common
    //!        Atom block: Specifies x y and z coordinates, element type, charge, valence, and other atom info
    //!        Bond block: specifies atoms in the bond, the bond order, ring/chain topology, stereochemistry, etc
    //!        MDL Properties block: specifies additional information.  Lines in this block are in the format of 'M  XYZ'
    //!          where XYZ is a label to describe what property is described
    //!
    //! @see @link example_sdf_ctab_handler.cpp @endlink
    //! @author geanesar
    //! @date May 06, 2016
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CTabHandler :
      public util::ObjectInterface
    {
    private:

    ////////////////////
    // helper classes //
    ////////////////////
      
    //////////
    // data //
    //////////
      
      //! whether the CTab has been read and contains valid data
      bool m_IsValid;

      //! header information (atom and bond counts)
      MdlHeader                     m_Header;
      
      //! atom information
      storage::Vector< AtomInfo>    m_AtomInfos;
      
      //! bond information
      storage::Vector< BondInfo>    m_BondInfos;

      //! atom mappings for e.g. reactions
      storage::Vector< size_t>      m_AtomMapping;

      //! Store whether configurational info was read
      bool                          m_AtomTypesRead;

      //! Store whether configurational info was read
      bool                          m_BondTypesRead;

      //! Store whether configurational info was read
      bool                          m_DoubleBondIsometryWasRead;

      //! Store whether configurational info was read
      bool                          m_ChiralityWasRead;

      //! cached MDL properties
      storage::Map< std::string, storage::Vector< std::string> > m_MDLProperties;

    public:

      //! s_Instance enables reading/writing ShPtr to this object
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      CTabHandler();

      //! @brief constructor from input stream
      explicit CTabHandler( std::istream &ISTREAM);

      //! @brief constructor from a pre-read set of lines
      CTabHandler( const storage::List< std::string> &LINES);

      //! @brief constructor from a list of AtomInfo and BondInfos
      CTabHandler
      ( 
        const storage::Vector< AtomInfo> &ATOM_INFOS,
        const storage::Vector< BondInfo> &BOND_INFOS
      );

      //! @brief Clone function
      //! @return pointer to new CTabHandler
      CTabHandler *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get all atom info stored in handler
      //! @return list of atom info stored in handler
      const storage::Vector< AtomInfo> &GetAtomInfo() const
      {
        return m_AtomInfos;
      }

      //! @brief get all Bond info stored in handler
      //! @return list of Bond info stored in handler
      const storage::Vector< BondInfo> &GetBondInfo() const
      {
        return m_BondInfos;
      }

      //! @brief get the cached MDL properties stored in this class
      //! @return a map with the property name (key) and associated values (values)
      const storage::Map< std::string, storage::Vector< std::string> > &GetCachedProperties() const
      {
        return m_MDLProperties;
      }

      //! @brief get the atom mapping vector
      //! @return a vector with atom mapping info, one entry for every atom; 0 = not mapped
      const storage::Vector< size_t> &GetAtomMapping() const
      {
        return m_AtomMapping;
      }

      //! @brief set the atom mapping value for a specific atom
      //! @param ATOM_NO the index of the atom to map
      //! @param VALUE the value to map; value of 0 = not mapped
      void SetAtomMapping( const size_t &ATOM_NO, const size_t &VALUE)
      {
        BCL_Assert
        (
          ATOM_NO < m_AtomMapping.GetSize(),
          "CTabHandler: tried to set the mapping for an invalid atom index.  Tried " + 
          util::Format()( ATOM_NO) + " but max value is " + util::Format()( m_AtomMapping.GetSize() - 1)
        );
        m_AtomMapping( ATOM_NO) = VALUE;
      }

      //! @brief clear/reset all useful data in this class
      void Reset()
      {
        m_IsValid = false;
        m_Header = MdlHeader();
        m_AtomInfos.Reset();
        m_BondInfos.Reset();
        m_AtomMapping.Reset();
        m_AtomTypesRead = false;
        m_BondTypesRead = false;
        m_DoubleBondIsometryWasRead = false;
        m_ChiralityWasRead = false;
      }

      //! @brief check if the handler is valid (non-empty, valid data)
      //! @return true if the mdl handler is in a valid state
      bool IsValid() const
      {
        return m_IsValid;
      }

      //! @brief check whether configuration information was read
      //! @return true if configuration information was read
      bool WasChiralityRead() const
      {
        return m_ChiralityWasRead;
      }

      //! @brief check whether configuration information was read
      //! @return true if configuration information was read
      bool WereAtomTypesRead() const
      {
        return m_AtomTypesRead;
      }

      //! @brief check whether configuration information was read
      //! @return true if configuration information was read
      bool WereBondTypesRead() const
      {
        return m_BondTypesRead;
      }

      //! @brief check whether configuration information was read
      //! @return true if configuration information was read
      bool WasDoubleBondIsometryRead() const
      {
        return m_DoubleBondIsometryWasRead;
      }

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief Read the CTab using iterators to strings
      //! @param LINE_BEGIN a line that represents a header/counts line
      //! @param LINE_END one-past-end of possible lines
      //! @details if the CTab ends before LINE_END, not all lines will be read
      //! @return an iterator to the first line that was not read
      storage::List< std::string>::const_iterator ReadCTab
      ( 
        const storage::List< std::string>::const_iterator &LINE_BEGIN,
        const storage::List< std::string>::const_iterator &LINE_END
      );

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &ReadCTab( std::istream &ISTREAM);

      //! @brief write to std::ostream in mdl format
      //! @param OSTREAM the stream to write to
      //! @param FORCE_WRITE_ATOM_TYPES whether to force writing BCL atom types (and other info)
      std::ostream &WriteCTab( std::ostream &OSTREAM, const bool &FORCE_WRITE_ATOM_TYPES = false) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief check if a line is the terminating line
      //! @return true if the line just contains 'M  END'
      static bool IsTerminalLine( const std::string &LINE);

    protected:

      //! @brief test whether a line contains only spaces or is otherwise empty
      //! @param STRING the string to test
      static bool ContainsNonspaceCharacters( const std::string &STRING);

      //! @brief apply instructions from an MDL property to this class
      //! @param PROP the property to apply
      void ApplyMDLProperty( const MdlProperty &PROP);

      //! @brief determine the MDL properties that must be used when writing CTab information
      //! @return a vector of MDL PropertyEnums that denote which properties must be computed
      storage::Vector< MdlProperty::PropertyEnum> GetNecessaryMDLProperties( const bool &FORCE_ATOM_TYPES) const;

      //! @brief from a set of atoms and bonds, determine whether atom types must be written
      //! @param ATOMS vector of atom info
      //! @param BONDS vector of bond info
      //! @return true if atom types must be written to ensure that they will be the same upon rereading the molecule
      static bool TestMustWriteAtomTypes
      (
        const storage::Vector< AtomInfo> &ATOMS,
        const storage::Vector< BondInfo> &BONDS
      );

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class CTabHandler

  } // namespace sdf
} // namespace bcl

#endif //BCL_SDF_CTAB_HANDLER_H_

