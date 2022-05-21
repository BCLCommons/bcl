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

#ifndef BCL_APP_DESCRIPTOR_SEQUENTIAL_FEATURE_SELECTION_H_
#define BCL_APP_DESCRIPTOR_SEQUENTIAL_FEATURE_SELECTION_H_

// include header of this class
#include "app/bcl_app_apps.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command.h"
#include "command/bcl_command_flag_interface.h"
#include "util/bcl_util.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DescriptorSequentialFeatureSelection
    //! @brief application DescriptorSequentialFeatureSelection creates specific descriptor code object files used
    //!        in descriptor selection process using app::TrainModel
    //!
    //! @details This application constructs a feature code representation for a particular step in the descriptor
    //!          optimization process based on the entire descriptor set, the number of iterations as step id.
    //!
    //! @author butkiem1
    //!
    //! @date 07/09/2010
    //!
    //! @see @link example_app_descriptor_sequential_feature_selection.cpp @endlink
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DescriptorSequentialFeatureSelection :
      public InterfaceRelease
    {

    private:

    //////////
    // data //
    //////////

      //! descriptor code object for features
      util::ShPtr< command::ParameterInterface> m_EntireDescriptorCodeObjectFileName;
      //! round number in descriptor optimization step
      //! output code object file for current iteration number
      util::ShPtr< command::ParameterInterface> m_FinalCodeObjectOutputFilename;

      //! round number in descriptor optimization step
      util::ShPtr< command::FlagInterface> m_FlagRoundNumber;

      //! write out initial descriptors for this iteration
      util::ShPtr< command::FlagInterface> m_FlagWriteInitialCodeObjectFile;

      //! flag for setting descriptor selection type
      util::ShPtr< command::FlagInterface> m_FlagDescriptorSelectionType;

      //! flag for storing meta data information for descriptor selection
      util::ShPtr< command::FlagInterface> m_FlagMetaDataStorage;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      DescriptorSequentialFeatureSelection();

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      DescriptorSequentialFeatureSelection *Clone() const
      {
        return new DescriptorSequentialFeatureSelection( *this);
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

      // This function is used to retrieve prior names for applications that are still valid (for backwards compatibility)
      // but which should not be displayed, e.g. for help
      storage::Vector< std::string> GetDeprecatedAppNames() const;

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

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
      static const ApplicationType DescriptorSequentialFeatureSelection_Instance;

    }; // DescriptorSequentialFeatureSelection

  } // namespace app
} // namespace bcl

#endif // BCL_APP_DESCRIPTOR_SEQUENTIAL_FEATURE_SELECTION_H_
