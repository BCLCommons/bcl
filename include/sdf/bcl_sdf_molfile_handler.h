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

#ifndef BCL_SDF_MOLFILE_HANDLER_H_
#define BCL_SDF_MOLFILE_HANDLER_H_

// include the namespace header
#include "bcl_sdf.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_sdf_atom_info.h"
#include "bcl_sdf_bond_info.h"
#include "bcl_sdf_ctab_handler.h"
#include "bcl_sdf_mdl_line_types.h"
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
    //! @class MolfileHandler
    //! @brief manages all MDL lines from a MDL file/block in terms of reading and parsing mdl lines of
    //!        certain line types as well as writing the lines back to create a MDL file from a MDL Handler.
    //!
    //! @see @link example_sdf_molfile_handler.cpp @endlink
    //! @author geanesar
    //! @date May 10, 2016
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MolfileHandler :
      public CTabHandler
    {

    public:

      //! number molecule descriptor lines
      static const size_t s_NumberDescriptionLines;

    private:

    //////////
    // data //
    //////////

      //! mdl description lines
      std::string                    m_Description;

    public:

      //! s_Instance enables reading/writing ShPtr to this object
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MolfileHandler();

      //! @brief constructor from a description and atom/bond info
      MolfileHandler
      (
        const std::string &DESCRIPTION,
        const storage::Vector< AtomInfo> &ATOM_INFOS,
        const storage::Vector< BondInfo> &BOND_INFOS
      );

      //! @brief constructor from input stream
      explicit MolfileHandler( std::istream &ISTREAM);

      //! @brief constructor from a pre-read set of lines
      MolfileHandler( const storage::List< std::string> &LINES);

      //! @brief Clone function
      //! @return pointer to new MolfileHandler
      MolfileHandler *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief clear all information from this class
      void Reset()
      {
        m_Description.clear();
        CTabHandler::Reset();
      }

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get all description lines stored in handler
      //! @return list of mdl description lines stored in handler
      const std::string &GetDescription() const
      {
        return m_Description;
      }

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &ReadMolfile( std::istream &ISTREAM);

      //! @brief write to std::ostream in mdl format
      std::ostream &WriteMolfile( std::ostream &OSTREAM, const bool &FORCE_WRITE_ATOM_TYPES = false) const;

      //! @brief read a molfile from a set of iterators
      //! @param LINE_BEGIN where to begin reading the molfile from
      //! @param LINE_END one-past-last iterator of where reading should happen
      //! @return an iterator to the first unread line
      storage::List< std::string>::const_iterator ReadMolfile
      (
        const storage::List< std::string>::const_iterator &LINE_BEGIN,
        const storage::List< std::string>::const_iterator &LINE_END
      );

      //! @brief write to std::ostream in mdl format
      //! @param OSTREAM the stream to write to
      //! @param DESCRIPTION the three-line description to use
      //! @param CTAB the connection table data for this molfile
      //! @param FORCE_WRITE_ATOM_TYPES whether to force writing atom types and other BCL info
      //! @return the output stream that was written to
      static std::ostream &WriteMolfile
      (
        std::ostream &OSTREAM,
        const std::string &DESCRIPTION,
        const CTabHandler &CTAB,
        const bool &FORCE_WRITE_ATOM_TYPES = false
      );

      //! @brief write to std::ostream in mdl format
      static std::ostream &WriteMolfile
      (
        std::ostream &OSTREAM,
        const std::string &DESCRIPTION,
        const storage::Vector< AtomInfo> &ATOM_INFO,
        const storage::Vector< BondInfo> &BOND_INFO
      );

      //! @brief set the description; checks # of newlines, ensures that it is exactly s_NumberDescriptionLines
      //! If the # of newlines is < s_NumberDescriptionLines - 1, adds new lines as necessary
      //! Any newlines >= s_NumberDescriptionLines are replaced with spaces
      static std::string StandardizeDescription( const std::string &DESCRIPTION);

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

    }; // class MolfileHandler

  } // namespace sdf
} // namespace bcl

#endif //BCL_SDF_MOLFILE_HANDLER_H_

