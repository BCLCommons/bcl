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

#ifndef BCL_APP_DESCRIPTOR_CONVERT_CODE_OBJECT_FILE_H_
#define BCL_APP_DESCRIPTOR_CONVERT_CODE_OBJECT_FILE_H_

// include the namespace header
#include "app/bcl_app_apps.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DescriptorConvertCodeObjectFile
    //! @brief application DescriptorConvertCodeObjectFile converts a code object file from an old formating version into the new
    //!        version given a complete  set of old version descriptors and new version descriptors.
    //!
    //! @details
    //!
    //! @author butkiem1
    //!
    //! @date 03/09/2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DescriptorConvertCodeObjectFile :
      public InterfaceRelease
    {

    private:

    //////////
    // data //
    //////////

      //! descriptor code object file name
      util::ShPtr< command::ParameterInterface> m_CodeObjectFile;

      //! descriptor code object file name
      util::ShPtr< command::ParameterInterface> m_ConvertedCodeObjectFile;

      //! mapping file - 1st is from, 2nd is to
      util::ShPtr< command::FlagInterface>      m_CodeMapping;

      //! split flag - if set, split the original code file into individual features, using the appropriate dataset
      //! retriever (to obtain sizes for all features)
      util::ShPtr< command::FlagInterface>      m_Split;

      //! merge flag - if set, merge the given code file with this code file, considering partials
      util::ShPtr< command::FlagInterface>      m_Merge;

      //! update flag - if set, update the original code file using the appropriate dataset retriever
      //! Descriptor files are automatically updated when generating datasets, but manually filtered
      //! descriptor files need to be updated prior to model training
      util::ShPtr< command::FlagInterface>      m_Update;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      DescriptorConvertCodeObjectFile();

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      DescriptorConvertCodeObjectFile *Clone() const
      {
        return new DescriptorConvertCodeObjectFile( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      // instantiate enumerator for PrepareSmallMoleculeEnsemble class
      static const ApplicationType DescriptorConvertCodeObjectFile_Instance;

    }; // DescriptorConvertCodeObjectFile

  } // namespace app
} // namespace bcl

#endif // BCL_APP_DESCRIPTOR_CONVERT_CODE_OBJECT_FILE_H_
