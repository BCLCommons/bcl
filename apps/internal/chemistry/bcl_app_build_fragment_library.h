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

#ifndef BCL_APP_BUILD_FRAGMENT_LIBRARY_H_
#define BCL_APP_BUILD_FRAGMENT_LIBRARY_H_

// include the interface for all apps
#include "app/bcl_app_apps.h"

// include other forward headers - sorted alphabetically
#include "chemistry/bcl_chemistry.fwd.hh"
#include "descriptor/bcl_descriptor.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_interface.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BuildFragmentLibrary
    //! @brief Application fragments molecules of a library into fragments
    //!
    //! @see @link example_app_build_fragment_library.cpp @endlink
    //! @author kothiwsk
    //! @date January 17, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BuildFragmentLibrary :
      public InterfaceRelease
    {

    private:

    //////////
    // data //
    //////////

      // commands

      //! output filename for fragments obtained by fragmentation
      util::ShPtr< command::FlagInterface> m_OutputFileFlag;

      //! set the maximum number of breakable bonds that a molecule can have for it to be fragmented
      util::ShPtr< command::FlagInterface> m_MaxNumberofBonds;

      //! set the maximum number of rotatable bonds that a fragment should have to be written out
      util::ShPtr< command::FlagInterface> m_MaxRotatableBonds;

      //! flag to find constitutions rather than configurations (default)
      util::ShPtr< command::FlagInterface> m_FindConstitutions;

      //! flag to ignore uncommmon rings (generally should be on unless generating rotamers specifically for those rings)
      util::ShPtr< command::FlagInterface> m_IgnoreUncommonRings;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      BuildFragmentLibrary();

      //! copy constructor; ignores everything but flags
      BuildFragmentLibrary( const BuildFragmentLibrary &PARENT);

    public:

      // instantiate enumerator for PrepareSmallMoleculeEnsemble class
      static const ApplicationType BuildFragmentLibrary_Instance;

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      BuildFragmentLibrary *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get the full description
      //! @return the full description
      std::string GetDescription() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetWebText() const;

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

    }; // BuildFragmentationLibrary

  } // namespace app
} // namespace bcl
#endif // BCL_APP_BUILD_FRAGMENT_LIBRARY_H_
