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

#ifndef BCL_APP_MOLECULE_PROPERTIES_H_
#define BCL_APP_MOLECULE_PROPERTIES_H_

// include the interface for all apps
#include "app/bcl_app_apps.h"

// include other forward headers - sorted alphabetically
#include "chemistry/bcl_chemistry.fwd.hh"
#include "descriptor/bcl_descriptor.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_command.h"
#include "io/bcl_io_ofstream.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_implementation.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeProperties
    //! @brief Application for performing operations with molecular properties
    //! @details Capabilities:
    //! MoleculeProperties performs operations on molecular properties,
    //! such as weight, HDon (# of hydrogen donors), etc (see -help for full list)
    //! Properties can be calculated, added to molecules and written out in file, analyzed with histograms or statistics,
    //! tabulated, removed, or renamed
    //! The properties can be numeric or string-based
    //!
    //! @see @link example_app_molecule_properties.cpp @endlink
    //! @author mendenjl
    //! @date 09/16/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeProperties :
      public InterfaceRelease
    {
    private:

    /////////////
    // typedef //
    /////////////

      // typedef for string properties, for keeping names reasonable length
      typedef util::Implementation< chemistry::StringPropertyInterface> StringProperty;

    //////////
    // data //
    //////////

      // histograms

      //! numeric property to take basic histogram of
      util::ShPtr< command::FlagInterface> m_HistogramNumericPropertyFlag;

      //! string property to take histogram of
      util::ShPtr< command::FlagInterface> m_HistogramStringPropertyFlag;

      //! filename to write histograms into
      util::ShPtr< command::FlagInterface> m_HistogramOutputFilenameFlag;

      //! a vector to hold all the string histograms desired
      mutable storage::Vector< storage::Map< std::string, size_t> > m_StringHistograms;

      //! make a vector to store the numeric histograms
      mutable storage::Vector< math::Histogram> m_NumericHistograms;

      //! a vector to store all the properties that will be looked up to make the table
      mutable storage::Vector< StringProperty> m_StringProperties;

      //! make a vector to hold the properties that will be used to calculate the numeric histograms
      mutable storage::Vector< descriptor::CheminfoProperty> m_NumericProperties;

      // statistics

      //! flag for properties to compute basic statistics on
      util::ShPtr< command::FlagInterface> m_AveMinMaxStdPropertyFlag;

      //! make a vector to store the mean/std statistics for each numeric property
      mutable storage::Vector< descriptor::CheminfoProperty> m_StatisticsProperties;
      mutable storage::Vector< math::RunningAverageSD< float> > m_MeanStds;
      mutable storage::Vector< math::RunningMinMax< float> > m_MinMaxs;

      // tabulation

      //! properties from which to make a table
      util::ShPtr< command::FlagInterface> m_TablePropertiesFlag;

      //! filename to write out property table to
      util::ShPtr< command::FlagInterface> m_TableOutputFilenameFlag;

      //! a vector to store all the properties that will be looked up to make the string histograms
      mutable storage::Vector< StringProperty> m_TableProperties;

      //! output file for table
      mutable io::OFStream m_TableOutputFile;

      // editing

      //! property names to remove (can be bcl-generated)
      util::ShPtr< command::FlagInterface> m_RemovePropertiesFlag;

      //! property names to remove (can be bcl-generated)
      util::ShPtr< command::FlagInterface> m_RemoveAllPropertiesFlag;

      //! any properties that should be renamed
      util::ShPtr< command::FlagInterface> m_RenamePropertiesFlag;

      //! any properties that should be added
      util::ShPtr< command::FlagInterface> m_AddPropertiesFlag;

      //! any strings that should be added as properties
      util::ShPtr< command::FlagInterface> m_AddPropertyStringsFlag;

      //! filename to write molecules into
      util::ShPtr< command::FlagInterface> m_OutputFilenameFlag;

      //! vector of properties to add
      mutable storage::List< storage::Pair< std::string, StringProperty> > m_PropertiesToAdd;

      //! properties to remove
      mutable storage::Vector< std::string> m_PropertiesToRemove;

      //! properties to rename
      mutable storage::Map< std::string, std::string> m_PropertiesToRename;

      //! whether to add an index property
      mutable bool m_AddIndexProperty;

      //! output file for molecules (if changing properties)
      mutable io::OFStream m_MoleculeOutputFile;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! @brief default constructor
      MoleculeProperties();

      //! @brief copy constructor; only copy flags for applications
      MoleculeProperties( const MoleculeProperties &PARENT);

    public:

      // instantiate enumerator for MoleculeProperties class
      static const ApplicationType MoleculeProperties_Instance;

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      MoleculeProperties *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief returns web text information
      //! @return text (html allowed but not required) that will be displayed on the website
      //! @note use confluence syntax, e.g. !my_image.png! to refer to an image under documentation/image
      const std::string &GetWebText() const;

      //! @brief check that all the parameter choices were valid
      //! @return true if all the parameter choices were valid
      bool CheckParametersAreAcceptable() const;

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

      //! @brief format a string for writing into a csv file
      //! @param STRING string to format
      //! @return the formatted string
      static std::string FormatCSVString( const std::string &STRING);

      //! @brief initialize the histograms and properties from the command line flags
      void Initialize() const;

      //! @brief add the data from a small molecule to the histograms
      //! @param MOLECULE the small molecule whose data to add
      //! @param MOLECULE_INDEX the index of the molecule currently operated on
      void AddData( chemistry::FragmentComplete &MOLECULE, const size_t &MOLECULE_INDEX) const;

    }; // MoleculeProperties

  } // namespace app
} // namespace bcl

#endif // BCL_APP_MOLECULE_PROPERTIES_H_
