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

#ifndef BCL_APP_MAKE_CHIMERIC_MOLECULE_H_
#define BCL_APP_MAKE_CHIMERIC_MOLECULE_H_

// include the interface for all apps
#include "app/bcl_app_apps.h"

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
    //! @class MakeChimericMolecule
    //! @brief Application for combining small molecules into a disconnected ensemble supermolecule
    //!
    //! @see @link example_app_make_chimeric_molecule.cpp @endlink
    //! @author brownbp1
    //! @date Sep 19, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MakeChimericMolecule :
      public Interface
    {
    public:

        //! A static instance of this class
        static const ApplicationType MakeChimericMolecule_Instance;

    private:

    //////////
    // data //
    //////////

      //! whether to recenter the molecules
      util::ShPtr< command::FlagInterface> m_RecenterFlag;

      //! output filename
      util::ShPtr< command::FlagInterface> m_OutputFilenameFlag;

      //! output file
      mutable io::OFStream m_OutputFile;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      MakeChimericMolecule();

      //! copy constructor; skips i/o streams
      MakeChimericMolecule( const MakeChimericMolecule &PARENT);

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      MakeChimericMolecule *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

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

    }; // MakeChimericMolecule

  } // namespace app
} // namespace bcl

#endif // BCL_APP_MAKE_CHIMERIC_MOLECULE_H_
