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

#ifndef BCL_CHEMISTRY_ROTAMER_LIBRARY_FILE_H_
#define BCL_CHEMISTRY_ROTAMER_LIBRARY_FILE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_rotamer_library_interface.h"
#include "command/bcl_command_parameter_check_file_in_search_path.h"
#include "io/bcl_io_stream_buffer_classes.h"
#include "io/bcl_io_stream_buffer_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RotamerLibraryFile
    //! @brief is a handler for storing and retrieving smallmolecule files from a file system.
    //!
    //! @details SmallMolecules are stored in a file structure that follows that logic:
    //!        {path/{SOURCE}filename.sdf} where the filename can have a prefix component {SOURCE}
    //!
    //! SOURCE is a prefix which can be empty
    //!
    //! In an file the key corresponds to the index of molecules which is auto incremented and is only numeric
    //! when initialized as attached, the largest key is located
    //! when initialized as created or overwrite, the key starts with 0
    //!
    //! @see @link example_chemistry_rotamer_library_file.cpp @endlink
    //! @author kothiwsk
    //! @date Jul 01, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RotamerLibraryFile :
      public RotamerLibraryInterface
    {

    protected:

    //////////
    // data //
    //////////

      //! prefix of file that contains storage
      mutable std::string m_Filename;

      //!
      size_t      m_FileNumber;

      // split the file into a format so that it can be iterated
      mutable storage::Vector< storage::Vector< std::string> > m_Lines;

      //! compression used for storing files
      io::StreamBufferClass m_Compression;

      std::string                       m_Alias;

    public:

    //////////
    // data //
    //////////

      //! store and retrieve instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_StoreInstance;
      static const util::SiPtr< const util::ObjectInterface> s_CODInstance;

      //! @brief create the command line file search path object
      static command::ParameterCheckFileInSearchPath GetRotamerFinder();

      //! @brief a stopwatch for all reading
      static util::Stopwatch s_ReadingTimer;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor, try to find the rotamer library
      RotamerLibraryFile();

      //! @brief Clone function
      //! @return pointer to new RotamerLibraryFile
      RotamerLibraryFile *Clone() const;

      //! @brief constructor, try to find the rotamer library
      RotamerLibraryFile( const std::string &ALIAS);

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

      //! @brief get the filename
      //! @return the filename
      const std::string &GetFilename() const
      {
        return m_Filename;
      }

      //! @brief get the filename
      //! @return the filename
      const std::string &GetPrefix() const
      {
        return m_Alias.empty() ? GetDefault() : m_Alias;
      }

      //! @brief get the filename
      //! @return the filename
      std::string GetFileExtension() const;

      //! @brief detect whether the rotamer library is defined
      bool IsDefined() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief get the directed-graph whose vertices are constitution ids and child nodes are substructures of parent node
      graph::ConstGraph< size_t, size_t> RetrieveRotamerSubstructureTree() const;

      //! @brief get the directed-graph whose vertices are constitution ids and child nodes are substructures of parent node
      graph::ConstGraph< size_t, size_t> RetrieveRotamerRequirementsGraph() const;

      //! @brief retrive constitution associated with id.
      storage::Vector< graph::ConstGraph< size_t, size_t> > RetrieveConstitution( const storage::Vector< size_t> &IDS) const;

      //! @brief retrieve configurations associated with the given constitution id.
      storage::Vector< storage::Set< size_t> > RetrieveConfigurationMapping() const;

      //! @brief get constitutions that are root in the substructure tree i.e. they are not contained in any other fragments
      storage::Vector< graph::ConstGraph< size_t, size_t> > GetRootConstitutions() const;

      //! @brief retrieve configurations associated with a given set of constitution ids
      FragmentEnsemble RetrieveAssociatedConfigurations( const storage::Set< size_t> &IDS) const;

      //! @brief get the bond angle map
      t_BondAngleMap GetBondAngleMap() const;

    private:

      //! @brief get the directed-graph whose vertices are constitution ids and child nodes are substructures of parent node
      void CreateImpl
      (
        const size_t &NUMBER_OF_ROOT_NODES,
        const storage::Vector< storage::Set< size_t> > &ROTAMER_SUBSTRUCTURE_TREE,
        const util::ShPtrVector< FragmentConstitutionShared> &CONSTITUTIONS,
        const storage::Vector< storage::Set< size_t> > &CONSTITUTION_TO_CONFORMATION_MAPPING,
        const FragmentEnsemble &CONFORMATIONS,
        const storage::Vector< size_t> &PARENTS
      ) const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief Load in necessary information, if needed
      void LoadFiles() const;

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook
      (
        const util::ObjectDataLabel &LABEL,
        std::ostream &ERR_STREAM
      );

      //! @brief get license agreement text for rotamer library
      void GetRotlibLicenseText() const;

    }; // class RotamerLibraryFile

  } // namespace chemistry

} // namespace bcl

#endif // BCL_CHEMISTRY_ROTAMER_LIBRARY_FILE_H_
