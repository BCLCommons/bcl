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

#ifndef BCL_APP_MOLECULE_SPLIT_H_
#define BCL_APP_MOLECULE_SPLIT_H_

// include the interface for all apps
#include "app/bcl_app_apps.h"
#include "app/bcl_app_interface_release.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "chemistry/bcl_chemistry_fragment_split_interface.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "io/bcl_io_file.h"
#include "storage/bcl_storage_template_instantiations.h"
#include "util/bcl_util_implementation.h"

namespace bcl
{
  namespace app
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeSplit
    //! @brief Application for editing small molecule ensembles
    //! @details some of the capabilities:
    //!           1. Merging ensembles
    //!           2. Adding/Removing hydrogens
    //!           3. Writing the molecules in an ensembles that...
    //!              A. do (and/or do not) contain particular fragments (given in a file)
    //!              B. do not correspond to any molecule in another ensemble (at the constitutional level)
    //!              C. they contain a particular misc. property
    //!              D. a particular misc property value is exactly (non)equal a string (e.g. Inhibitor)
    //!              E. a particular misc property value contains a particular string (e.g. nan)
    //!              F. a particular small molecule misc property value is defined and satisfies a comparison (e.g. EC50 < 0.1)
    //!              G. are complexes
    //!
    //! @see @link example_app_molecule_split.cpp @endlink
    //! @author mendenjl, brownbp1
    //! @date May 15, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeSplit :
      public InterfaceRelease
    {
    public:

        //! A static instance of this class
        static const ApplicationType MoleculeSplit_Instance;

    private:

    //////////
    // data //
    //////////

      //! file indicating implementation of chemistry::FragmentSplitInterface
      util::ShPtr< command::FlagInterface> m_ImplementationFlag;

      //! flag indicating the minimum number of atoms a split fragment must have to be output
      util::ShPtr< command::FlagInterface> m_MinFragSizeFlag;

      //! whether to recenter the molecules
      util::ShPtr< command::FlagInterface> m_RecenterFlag;

      //! keep the MDL properties of original molecule on all split components
      util::ShPtr< command::FlagInterface> m_PreserveMDLPropertiesFlag;

      //! output filename
      util::ShPtr< command::FlagInterface> m_OutputFilenameFlag;

      mutable util::Implementation< chemistry::FragmentSplitInterface> m_Split; //!< obtains a implementation

      //! output file
      mutable io::OFStream m_OutputFile;

      //! molecule index currently being split
      mutable size_t m_MoleculeIndex;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      MoleculeSplit();

      //! copy constructor; skips i/o streams
      MoleculeSplit( const MoleculeSplit &PARENT);

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      MoleculeSplit *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

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
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief writes a single molecule out
      //! @param MOLECULE the molecule or fragment to write
      void Write( chemistry::FragmentComplete &MOLECULE) const;

    }; // MoleculeSplit

  } // namespace app
} // namespace bcl

#endif // BCL_APP_MOLECULE_SPLIT_H_
