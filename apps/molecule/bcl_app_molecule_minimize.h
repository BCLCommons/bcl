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

#ifndef BCL_APP_MOLECULE_MINIMIZE_H_
#define BCL_APP_MOLECULE_MINIMIZE_H_

// include the interface for all apps
#include "app/bcl_app_apps.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_interface.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "io/bcl_io_file.h"

// includes from rdkit
//#include "GraphMol/RWMol.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeMinimize
    //! @brief Application for processing performing energy minimization of molecular geometries.
    //! @details This application makes use of the RDKit external library. It converts a BCL molecule into
    //! an RDKit molecule, utilizes the RDkit ForceFields framework to perform the energy minimization, and
    //! then converts the resultant molecule back to a BCL molecule before output.
    //!
    //! @see @link example_app_molecule_minimize.cpp @endlink
    //! @author brownbp1
    //! @date Aug 27, 2022
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeMinimize :
      public InterfaceRelease
    {
    private:

    //////////
    // data //
    //////////

      //! Output filename base
      util::ShPtr< command::FlagInterface> m_Output;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      MoleculeMinimize();

    public:

      // instantiate enumerator for MoleculeMinimize class
      static const ApplicationType MoleculeMinimize_Instance;

      //! @brief Clone function
      //! @return pointer to new MoleculeMinimize
      MoleculeMinimize *Clone() const;

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

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

//      //! @brief minimize with constraints
//      void MinimizeWithConstraints( std::shared_ptr< ::RDKit::RWMol> &MOL) const;

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

    }; // MoleculeMinimize

  } // namespace app
} // namespace bcl

#endif // BCL_APP_MOLECULE_MINIMIZE_H_
