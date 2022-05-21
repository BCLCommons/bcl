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

#ifndef BCL_SDF_MDL_HANDLER_H_
#define BCL_SDF_MDL_HANDLER_H_

// include the namespace header
#include "bcl_sdf.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_sdf_atom_info.h"
#include "bcl_sdf_bond_info.h"
#include "bcl_sdf_ctab_handler.h"
#include "bcl_sdf_mdl_line_types.h"
#include "bcl_sdf_molfile_handler.h"
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
    //! @class MdlHandler
    //! @brief manages all MDL lines from a MDL file/block in terms of reading and parsing mdl lines of
    //!        certain line types as well as writing the lines back to create a MDL file from a MDL Handler.
    //!
    //! @see @link example_sdf_mdl_handler.cpp @endlink
    //! @author butkiem1, loweew, woetzen, mendenjl
    //! @date May 7, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MdlHandler :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      //! determine whether the SDF file is valid
      mutable bool m_IsValid;

      //! Molfile handler for
      mutable MolfileHandler m_Molfile;

      //! map that stores associated misc properties with their identifiers
      mutable storage::Map< std::string, std::string> m_MiscProperties;

      //! whether the buffered lines were parsed yet
      mutable bool m_WasParsed;

      //! buffered lines from an SDF file
      storage::List< std::string> m_Lines;

    public:

      //! s_Instance enables reading/writing ShPtr to this object
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MdlHandler();

      //! @brief default constructor
      MdlHandler( std::istream &ISTREAM);

      //! @brief constructor from a pre-parsed set of lines
      MdlHandler( const storage::List< std::string> &LINES);

      //! @brief Clone function
      //! @return pointer to new MdlHandler
      MdlHandler *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief clear all information from this class
      void Reset()
      {
        m_IsValid = false;
        m_Molfile.Reset();
        m_MiscProperties.Reset();
        m_WasParsed = false;
        m_Lines.Reset();
      }

      //! @brief return if the class is valid or not
      bool IsValid() const
      {
        FinalizeParsing();
        return m_IsValid && m_Molfile.IsValid();
      }

      //! @brief get a description of the molecule
      //! @return a string containing the molecule's description
      const std::string &GetDescription() const
      {
        FinalizeParsing();
        return m_Molfile.GetDescription();
      }

      //! @brief get all atom info stored in handler
      //! @return list of atom info stored in handler
      const storage::Vector< AtomInfo> &GetAtomInfo() const
      {
        FinalizeParsing();
        return m_Molfile.GetAtomInfo();
      }

      //! @brief get all Bond info stored in handler
      //! @return list of Bond info stored in handler
      const storage::Vector< BondInfo> &GetBondInfo() const
      {
        FinalizeParsing();
        return m_Molfile.GetBondInfo();
      }

      //! @brief get the cached MDL properties from the contained molfile
      //! @return the cached MDL properties
      const storage::Map< std::string, storage::Vector< std::string> > &GetCachedProperties() const
      {
        FinalizeParsing();
        return m_Molfile.GetCachedProperties();
      }

      //! @brief get the atom mapping from the underlying molfile
      //! @return the atom mapping vector
      const storage::Vector< size_t> &GetAtomMapping() const
      {
        FinalizeParsing();
        return m_Molfile.GetAtomMapping();
      }

      //! @brief get all MiscProp lines stored in handler
      //! @return list of mdl MiscProp lines stored in handler
      const storage::Map< std::string, std::string> &GetMiscProperties() const
      {
        FinalizeParsing();
        return m_MiscProperties;
      }

      //! @brief check whether configuration information was read
      //! @return true if configuration information was read
      bool WasChiralityRead() const
      {
        FinalizeParsing();
        return m_Molfile.WasChiralityRead();
      }

      //! @brief check whether configuration information was read
      //! @return true if configuration information was read
      bool WereAtomTypesRead() const
      {
        FinalizeParsing();
        return m_Molfile.WereAtomTypesRead();
      }

      //! @brief check whether configuration information was read
      //! @return true if configuration information was read
      bool WereBondTypesRead() const
      {
        FinalizeParsing();
        return m_Molfile.WereBondTypesRead();
      }

      //! @brief check whether configuration information was read
      //! @return true if configuration information was read
      bool WasDoubleBondIsometryRead() const
      {
        FinalizeParsing();
        return m_Molfile.WasDoubleBondIsometryRead();
      }

      const MolfileHandler &GetMolfile() const
      {
        FinalizeParsing();
        return m_Molfile;
      }

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief iterate through buffered input and convert to useful data
      void FinalizeParsing() const;

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &ReadFromSDF( std::istream &ISTREAM);

      //! @brief write to std::ostream in mdl format
      std::ostream &WriteToSDF( std::ostream &OSTREAM) const;

      //! @brief write to std::ostream in mdl format
      //! @param OSTREAM the output stream to write to
      //! @param MOLFILE the molfile to write
      //! @param MISC_PROPERTIES the miscellaneous properties to write
      //! @return the output stream that was written
      static std::ostream &WriteToSDF
      (
        std::ostream &OSTREAM,
        const MolfileHandler &MOLFILE,
        const storage::Map< std::string, std::string> &MISC_PROPERTIES
      );

      //! @brief write to std::ostream in mdl format
      //! @param OSTREAM the output stream to write to
      //! @param DESCRIPTION the 3-line description to use
      //! @param ATOM_INFO the atom infos to write
      //! @param BOND_INFO the bond infos to write
      //! @param MISC_PROPERTIES the miscellaneous properties to write
      static std::ostream &WriteToSDF
      (
        std::ostream &OSTREAM,
        const std::string &DESCRIPTION,
        const storage::Vector< AtomInfo> &ATOM_INFO,
        const storage::Vector< BondInfo> &BOND_INFO,
        const storage::Map< std::string, std::string> &MISC_PROPERTIES
      );

      //! @brief write a simple string that should be unique for a constitution, independent of H
      //! @note hash is dependent on ordering of atoms; use accordingly
      static std::string CreateConstitutionalHashString
      (
        const storage::Vector< AtomInfo> &ATOM_INFO,
        const storage::Vector< BondInfo> &BOND_INFO,
        chemistry::ConfigurationalBondTypeData::Data BOND_TYPE = chemistry::ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness
      );

      //! @brief write a simple string that should be unique for this molecular configuration
      //! @note hash is dependent on ordering of atoms, etc.; use accordingly
      static std::string CreateConfigurationalHashString
      (
        const storage::Vector< AtomInfo> &ATOM_INFO,
        const storage::Vector< BondInfo> &BOND_INFO,
        chemistry::ConfigurationalBondTypeData::Data BOND_TYPE = chemistry::ConfigurationalBondTypeData::e_BondOrderAmideWithIsometryOrAromaticWithRingness
      );

      //! @brief write a simple string that should be unique for this molecular conformation
      //! @note hash is dependent on orientation of molecule, ordering of atoms, etc.; use accordingly
      static std::string CreateConformationalHashString
      (
        const storage::Vector< AtomInfo> &ATOM_INFO,
        const storage::Vector< BondInfo> &BOND_INFO
      );

      //! @brief get access to a global flag defining whether hydrogens should be added
      //! @return access to a global flag defining whether hydrogens should be added
      static util::ShPtr< command::FlagInterface> &GetAddAtomMdlLineFlag();

      //! @brief add hydrogen handling preferences to the command line flag
      //! @param CMD command to add the hydrogen handling preference flags to
      static void AddAtomMdlLineFlag( command::Command &CMD);

      //! @brief return the data label, if line is a data label line
      //! @param LINE line from mdl section
      //! @return string that constains datalable, string will be empty for non-data lable lines
      static std::string GetMDLDataLabel( const std::string &LINE);

      //! @brief write misc properties as mdl lines into std::ostream
      //! @param OSTREAM ostream to write MiscProperties to
      static void WriteMiscProperties
      (
        std::ostream &OSTREAM,
        const storage::Map< std::string, std::string> &MISC_PROPERTIES
      );

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief check if a line is the terminating line
      //! @return true if the line just contains $$$$
      static bool IsTerminalLine( const std::string &LINE);

    protected:

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

      //! @brief test whether a line contains only spaces or is otherwise empty
      //! @param STRING the string to test
      static bool ContainsNonspaceCharacters( const std::string &STRING);

    }; // class MdlHandler

  } // namespace sdf
} // namespace bcl

#endif //BCL_SDF_MDL_HANDLER_H_

