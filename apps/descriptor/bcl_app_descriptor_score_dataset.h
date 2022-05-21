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

#ifndef BCL_APP_DESCRIPTOR_SCORE_DATASET_H_
#define BCL_APP_DESCRIPTOR_SCORE_DATASET_H_

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
    //! @class DescriptorScoreDataset
    //! @brief Application for scoring individual feature columns. Used during descriptor development and selection
    //!
    //! @see @link example_app_descriptor_score_dataset.cpp @endlink
    //! @author mendenjl
    //! @date Apr 18, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DescriptorScoreDataset :
      public InterfaceRelease
    {

    private:

    //////////
    // data //
    //////////

      //! flag for the data set
      util::ShPtr< command::FlagInterface> m_FlagDataSet;
      //! flag for the list of feature names/types
      util::ShPtr< command::FlagInterface> m_FlagFeatureCode;
      //! flag for the list of result names/types
      util::ShPtr< command::FlagInterface> m_FlagResultCode;
      //! flag for the output filename
      util::ShPtr< command::FlagInterface> m_FlagOutputFilename;
      //! flag for the score filename
      util::ShPtr< command::FlagInterface> m_FlagScore;
      //! flag to produce terse output (smaller file, but no sorted descriptors)
      util::ShPtr< command::FlagInterface> m_FlagTerse;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      DescriptorScoreDataset();

    public:

      //! @brief Clone function
      //! @return pointer to new DescriptorScoreDataset
      DescriptorScoreDataset *Clone() const;

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
        return storage::Vector< std::string>( size_t( 1), "ScoreFeatureResultDataset");
      }

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

      // instantiate enumerator for this class
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

    }; // DescriptorScoreDataset

  } // namespace app
} // namespace bcl

#endif // BCL_APP_DESCRIPTOR_SCORE_DATASET_H_
