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

#ifndef BCL_APP_DESCRIPTOR_REFINE_BY_SCORE_H_
#define BCL_APP_DESCRIPTOR_REFINE_BY_SCORE_H_

// include the interface for all apps
#include "app/bcl_app_apps.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_command.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DescriptorRefineByScore
    //! @brief Application for refining a descriptor file using a descriptor score file
    //!
    //! @see @link example_app_descriptor_refine_by_score.cpp @endlink
    //! @author mendenjl
    //! @date Apr 18, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DescriptorRefineByScore :
      public InterfaceRelease
    {

    private:

    //////////
    // data //
    //////////

      //! flag for the data set score file
      util::ShPtr< command::FlagInterface> m_FlagDataSetScoreFilename;
      //! flag for the method of choosing descriptors
      util::ShPtr< command::FlagInterface> m_FlagSelectionMethod;
      //! flag for the output filename
      util::ShPtr< command::FlagInterface> m_FlagOutputFilename;
      //! flag to output information about the selection to m_FlagOutputFilename - last extension + info
      util::ShPtr< command::FlagInterface> m_FlagOutputInfo;

      //! flag for score file comparison
      util::ShPtr< command::FlagInterface> m_FlagCompareScoreFiles;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      DescriptorRefineByScore();

    public:

      //! @brief Clone function
      //! @return pointer to new RefineDescriptors
      DescriptorRefineByScore *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      // This function is used to retrieve prior names for applications that are still valid (for backwards compatibility)
      // but which should not be displayed, e.g. for help
      storage::Vector< std::string> GetDeprecatedAppNames() const;

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

      // instantiate enumerator for DescriptorRefineByScore class
      static const ApplicationType s_Instance;

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

    }; // DescriptorRefineByScore

  } // namespace app
} // namespace bcl

#endif // BCL_APP_DESCRIPTOR_REFINE_BY_SCORE_H_
