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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#ifndef BCL_APP_GENERATE_HIERARCHICAL_TREE_H_
#define BCL_APP_GENERATE_HIERARCHICAL_TREE_H_

// include the interface for all apps
#include "app/bcl_app_apps.h"

// include other forward headers - sorted alphabetically
#include "chemistry/bcl_chemistry.fwd.hh"
#include "descriptor/bcl_descriptor.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "storage/bcl_storage_template_instantiations.h"
#include "util/bcl_util_implementation.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class GenerateHierarchicalTree
    //! @brief application for loading the rotamer library into a database or into a file tree structure
    //! @details Takes a rotamer library and loads it into csd_rotamer_library database or creates a flat file version of
    //!          rotamer library that can be searched in a tree like manner.
    //! @author kothiwsk
    //! @date 03/07/2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API GenerateHierarchicalTree :
      public InterfaceRelease
    {

    private:

    //////////
    // data //
    //////////

      //! scaffold sdf
      util::ShPtr< command::FlagInterface> m_LibraryFlag;

      //! output sdf file
      util::ShPtr< command::FlagInterface> m_FormatFlag;

      //! prune flag
      util::ShPtr< command::FlagInterface> m_PruneFlag;

      //! rotamer library interface format
      mutable util::Implementation< chemistry::RotamerLibraryInterface> m_Format;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      GenerateHierarchicalTree();

      //! copy constructor
      GenerateHierarchicalTree( const GenerateHierarchicalTree &APP);

    public:

      // instantiate enumerator for PrepareSmallMoleculeEnsemble class
      static const ApplicationType GenerateHierarchicalTree_Instance;

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      GenerateHierarchicalTree *Clone() const;

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

    }; // RotamerLibraryDb

  } // namespace app
} // namespace bcl
#endif // BCL_APP_GENERATE_HIERARCHICAL_TREE_H_
