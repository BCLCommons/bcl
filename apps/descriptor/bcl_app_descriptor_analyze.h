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

#ifndef BCL_APP_DESCRIPTOR_ANALYZE_H_
#define BCL_APP_DESCRIPTOR_ANALYZE_H_

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
    //! @class DescriptorAnalyze
    //! @details Object data labels will be converted into Partial representations for every descriptor column.
    //!          Then the percent overlap is calculated for every pair.
    //!
    //!          overlap matrix has following format:
    //!
    //!                                     cobj1 cobj2  cobj3
    //!                              cobj4  0.1    0.2    0.1
    //!                              cobj5  0.5    0.3    0.8
    //!                              cobj6  0.6    0.4    0.2
    //!
    //! @author butkiem1
    //! @see @link example_app_descriptor_analyze.cpp @endlink
    //! @date 08/03/2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DescriptorAnalyze :
      public InterfaceRelease
    {

    private:

    //////////
    // data //
    //////////

      //! descriptor code object file name
      util::ShPtr< command::FlagInterface> m_CodeObjectFilesRow;

      //! descriptor code object file name
      util::ShPtr< command::FlagInterface> m_CodeObjectFilesColumn;

      //! descriptor code object file name
      util::ShPtr< command::FlagInterface> m_OutputMatrixFileName;

      //! dataset retriever; used to determine the size of each feature
      util::ShPtr< command::FlagInterface> m_DatasetRetriever;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      DescriptorAnalyze();

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      DescriptorAnalyze *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      // This function is used to retrieve prior names for applications that are still valid (for backwards compatibility)
      // but which should not be displayed, e.g. for help
      storage::Vector< std::string> GetDeprecatedAppNames() const
      {
        return storage::Vector< std::string>( size_t( 1), "AnalyzeCodeObjectFile");
      }

    ////////////////
    // operations //
    ////////////////

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

    public:

      // instantiate enumerator for PrepareSmallMoleculeEnsemble class
      static const ApplicationType DescriptorAnalyze_Instance;

    }; // AnalyzeCodeObjectFile

  } // namespace app
} // namespace bcl

#endif // BCL_APP_DESCRIPTOR_ANALYZE_H_
