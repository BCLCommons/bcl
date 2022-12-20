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

#ifndef BCL_SMILES_RDKIT_SMILES_PARSER_H_
#define BCL_SMILES_RDKIT_SMILES_PARSER_H_

// include the namespace header
#include "bcl_smiles.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "sdf/bcl_sdf_atom_info.h"
#include "sdf/bcl_sdf_bond_info.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace smiles
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RdkitSmilesParser
    //! @brief This class handles reading and writing of SMILES/SMARTS formatted files using the RDKit external
    //! library and functionality to convert between RDKit and BCL data structures.
    //!
    //! @see @link example_smiles_rdkit_smiles_parser.cpp @endlink
    //! @author brownbp1
    //! @date Aug 29, 2022
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RdkitSmilesParser :
      public util::ObjectInterface
    {

    ////////////////////
    // helper classes //
    ////////////////////
      
    //////////
    // data //
    //////////
      
    public:

      // enum for filetypes
      enum FileType
      {
        e_SMILES = 0,
        e_SMARTS,
        s_TotalFileTypes
      };

    private:

      // header line from input file
      storage::Vector< std::string> m_Header;

      // all data read in from input file
      storage::Vector< storage::Map< std::string, std::string>> m_FileData;

      // molecules that have been read in from input file
      storage::Vector< chemistry::FragmentComplete> m_Molecules;

      // file type to read/write
      FileType m_FileType;

    public:

      //! s_Instance enables reading/writing ShPtr to this object
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      RdkitSmilesParser();

      //! @brief constructor from input stream
      explicit RdkitSmilesParser( std::istream &ISTREAM);

      //! @brief constructor from a pre-read set of lines
      RdkitSmilesParser( const storage::List< std::string> &LINES);

      //! @brief constructor from an input stream and file type specification
      RdkitSmilesParser( std::istream &ISTREAM, const FileType &TYPE);

      //! @brief constructor from a pre-read set of lines and file type specification
      RdkitSmilesParser( const storage::List< std::string> &LINES, const FileType &TYPE);

//      //! @brief write constructor from input molecules and file-type
//      RdkitSmilesParser
//      (
//        std::ostream &OSTREAM,
//        const storage::Vector< chemistry::FragmentComplete> &MOLECULES,
//        const FileType &TYPE
//      );

      //! @brief Clone function
      //! @return pointer to new RdkitSmilesParser
      RdkitSmilesParser *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief return the file type
      //! returns the file type as a string
      static const std::string &GetFileTypeAsString( const FileType &ENUM);

      //! DataEnum simplifies the usage of the Data enum of this class
      typedef util::WrapperEnum< FileType, &GetFileTypeAsString, s_TotalFileTypes> DataEnum;

      //! @brief returns the file data labels from the input file header
      //! @return a vector of strings where each string is a column header label
      const storage::Vector< std::string> &GetHeader() const
      {
        return m_Header;
      }

      //! @brief returns the file data label for the molecule string from the header
      //! @return the string corresponding to the molecule input type (e.g., Smiles or SMILES)
      const std::string &GetMolStrLabel() const
      {
        return m_Header( 0);
      }

      //! @brief returns the contents of the file
      //! @return a vector of maps where each map matches a column header label (key)
      //! to a data value (value); SMILES/SMARTs string must always be first label
      const storage::Vector< storage::Map< std::string, std::string>> &GetFileData() const
          {
            return m_FileData;
      }

      //! @brief return the molecules that were read in from the input file
      //! @return pointers to molecules
      storage::Vector< chemistry::FragmentComplete> &GetMolecules()
      {
        return m_Molecules;
      }

      //! @brief return the file type
      //! returns the file type as an enum
      const FileType &GetFileType() const
      {
        return m_FileType;
      }

      //! @brief set the file type
      void SetFileType( const FileType &ENUM)
      {
        m_FileType = ENUM;
      }

      //! @brief set the header
      void SetHeader( const storage::Vector< std::string> &HEADER)
      {
        m_Header = HEADER;
      }

      void SetMolecules( const storage::Vector< chemistry::FragmentComplete> &MOLECULES)
      {
        m_Molecules = MOLECULES;
      }

      //! @brief clear/reset all useful data in this class
      void ResetAll()
      {
        m_Header = storage::Vector< std::string>();
        m_FileData = storage::Vector< storage::Map< std::string, std::string > >();
        m_Molecules = storage::Vector< chemistry::FragmentComplete>();
        m_FileType = e_SMILES;
      }

      void ResetHeader()
      {
        m_Header = storage::Vector< std::string>();
      }

      void ResetFileData()
      {
        m_FileData = storage::Vector< storage::Map< std::string, std::string> >();
      }

      void ResetMolecules()
      {
        m_Molecules = storage::Vector< chemistry::FragmentComplete>();
      }

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief Read the SMILES/SMARTS file using iterators to strings
      //! @param LINE_BEGIN a line that represents a header/counts line
      //! @param LINE_END one-past-end of possible lines
      //! @param READ_TO_DATA_ONLY if true, do not store molecules in memory, but populate
      //! the file data member; useful if you are reading many molecules but working on
      //! one at a time.
      //! @details if the file ends before LINE_END, not all lines will be read
      //! @return an iterator to the first line that was not read
      storage::List< std::string>::const_iterator ReadFile
      ( 
        const storage::List< std::string>::const_iterator &LINE_BEGIN,
        const storage::List< std::string>::const_iterator &LINE_END,
        const bool READ_TO_DATA_ONLY = true
      );

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @param READ_TO_DATA_ONLY if true, do not store molecules in memory, but populate
      //! the file data member; useful if you are reading many molecules but working on
      //! one at a time.
      //! @return istream which was read from
      std::istream &ReadFile( std::istream &ISTREAM, const bool READ_TO_DATA_ONLY = true);

      //! @brief write to std::ostream in SMILES or SMARTS format
      //! @param OSTREAM the stream to write to
      std::ostream &WriteFile( std::ostream &OSTREAM) const;

      //! @brief write header lines to std::ostream to top the output SMILES or SMARTS file
      //! @param OSTREAM the stream to write to
      std::ostream &WriteHeaderLine( std::ostream &OSTREAM, std::string const &DELIMITER = std::string( ",")) const;

      //! @brief write to std::ostream to fill molecule and property data into the output SMILES or SMARTS file
      //! @param OSTREAM the stream to write to
      //! @param MOL molecule data to write to file; only properties in m_Headers will be written
      std::ostream &WriteNonHeaderLine
      (
        std::ostream &OSTREAM,
        const chemistry::FragmentComplete &MOL,
        const bool INCLUDE_PROPERTIES = false,
        const std::string &DELIMITER = std::string( ",")
      ) const;

      //! @brief build a molecule from the file data member
      //! @param MOL_INDEX the index of the file data that will be filled
      //! @returns molecule filled from data SMILES/SMARTS and any properties
      chemistry::FragmentComplete FillMolFromData( const size_t MOL_INDEX) const;

      //! @brief write to std::ostream directly from a molecule without object construction
      //! @param OSTREAM the stream to write to
      //! @param MOL the molecule whose data we will write
      //! @returns a stream with the new molecule SMILES; note that this write function
      //! does not write properties to file because no header line is produced
      static std::ostream &WriteSMILESFromMol
      (
        std::ostream &OSTREAM,
        const chemistry::FragmentComplete &MOL
      );

      //! @brief write to std::ostream directly from a molecule without object construction
      //! @param OSTREAM the stream to write to
      //! @param NAME the molecule name
      //! @param ATOM_INFO the atominfo describing the molecule
      //! @param BOND_INFO the bondingo describing the molecule
      //! @returns a stream with the new molecule SMILES; note that this write function
      //! does not write properties to file because no header line is produced
      static std::ostream &WriteSMILESFromMolInfo
      (
        std::ostream &OSTREAM,
        const std::string &NAME,
        const storage::Vector< sdf::AtomInfo> &ATOM_INFO,
        const storage::Vector< sdf::BondInfo> &BOND_INFO
      );

      //! @brief write to std::ostream directly from a molecule without object construction
      //! @param OSTREAM the stream to write to
      //! @param MOL the molecule whose data we will write
      //! @returns a stream with the new molecule SMARTS; note that this write function
      //! does not write properties to file because no header line is produced
      static std::ostream &WriteSMARTSFromMol
      (
        std::ostream &OSTREAM,
        const chemistry::FragmentComplete &MOL
      );

      //! @brief write to std::ostream directly from a molecule without object construction
      //! @param OSTREAM the stream to write to
      //! @param NAME the molecule name
      //! @param ATOM_INFO the atominfo describing the molecule
      //! @param BOND_INFO the bondingo describing the molecule
      //! @returns a stream with the new molecule SMARTS; note that this write function
      //! does not write properties to file because no header line is produced
      static std::ostream &WriteSMARTSFromMolInfo
      (
        std::ostream &OSTREAM,
        const std::string &NAME,
        const storage::Vector< sdf::AtomInfo> &ATOM_INFO,
        const storage::Vector< sdf::BondInfo> &BOND_INFO
      );

      //! @brief convert input molecule into SMILES string
      //! @param MOL molecule to be converted to SMILES string
      static std::string ConvertMolToSMILES( const chemistry::FragmentComplete &MOL);

      //! @brief convert input molecule into SMARTS string
      //! @param MOL molecule to be converted to SMARTS string
      static std::string ConvertMolToSMARTS( const chemistry::FragmentComplete &MOL);

      //! @brief convert SMILES string to molecule
      //! @param MOL molecule to be converted to SMILES string
      static chemistry::FragmentComplete ConvertSMILESToMOL
      (
        const std::string &SMILES,
        const bool IS_SMARTS = false,
        const bool ADD_H = false,
        const bool GEN3D = false,
        const size_t GEOOPT_ITERATIONS = size_t( 0),
        const std::string &GEOOPT_MMFF94S_VARIANT = "MMFF94s",
        const float GEOOPT_NONBONDED_THRESH = float( 10.0),
        const bool GEOOPT_IGNORE_INTERFRAG_INTERACTIONS = true
      );

      //! @brief convert SMARTS string to molecule
      //! @param MOL molecule to be converted to SMARTS string
      static chemistry::FragmentComplete ConvertSMARTSToMOL
      (
        const std::string &SMILES,
        const bool ADD_H = false,
        const bool GEN3D = false,
        const size_t GEOOPT_ITERATIONS = size_t( 0),
        const std::string &GEOOPT_MMFF94S_VARIANT = "MMFF94s",
        const float GEOOPT_NONBONDED_THRESH = float( 10.0),
        const bool GEOOPT_IGNORE_INTERFRAG_INTERACTIONS = true
      );

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief test whether a line contains only spaces or is otherwise empty
      //! @param STRING the string to test
      static bool ContainsNonspaceCharacters( const std::string &STRING);

      //! @brief check validity of first column header name; this column MUST
      //! indicate the output file type (SMILES or SMARTS)
      //! @returns true if valid; false otherwise
      bool ContainsValidHeader() const;

      //! @brief check if there is a header; assume that if the first line
      //! of the input file contains only one column that that column is intended
      //! for SMILES or SMARTS IDs only and that the first row should not be
      //! interpreted as a header row
      bool ContainsHeaderLine() const;

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

    }; // class RdkitSmilesParser

  } // namespace smiles
} // namespace bcl

#endif //BCL_SMILES_RDKIT_SMILES_PARSER_H_

