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

#ifndef BCL_APP_MOLECULE_UNIQUE_H_
#define BCL_APP_MOLECULE_UNIQUE_H_

// include the interface for all apps
#include "app/bcl_app_apps.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "command/bcl_command_flag_interface.h"
#include "io/bcl_io_ofstream.h"
#include "util/bcl_util_stopwatch.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeUnique
    //! @brief Application for removing molecules duplicated in an ensemble
    //! @details Can remove repeated molecules, with a repeat defined in any of the following ways:
    //!          1. Same Atoms, Connectivity, Bond orders, and 3D coordinates (default)
    //!          2. Same Atoms, Connectivity, and Stereochemistry (e.g. same configuration)
    //!          3. Same Atoms & Connectivity (e.g. same constitution)
    //!
    //! @see @link example_app_molecule_unique.cpp @endlink
    //! @author mendenjl
    //! @date May 14, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeUnique :
      public InterfaceRelease
    {
    private:

    //////////
    // data //
    //////////

      //! whether to remove duplicate molecules (same atom types, bond types, and positions) from files
      util::ShPtr< command::FlagInterface> m_CompareFlag;

      //! Where to place the molecules
      util::ShPtr< command::FlagInterface> m_OutputFilenameFlag;

      //! Where to place the molecules
      util::ShPtr< command::FlagInterface> m_OutputDuplicatesFilenameFlag;

      //! Bin size to determine conformation
      util::ShPtr< command::FlagInterface> m_ConformerComparerFlag;

      //! Merge descriptors for duplicated molecules - in case of duplicates, the first instances values are retained
      util::ShPtr< command::FlagInterface> m_MergeDescriptorsFlag;

      //! Overwrite descriptors for duplicated molecules, only the last instances values are retained
      util::ShPtr< command::FlagInterface> m_OverwriteDescriptorsFlag;

      //! Require the name to be identical for merging. This is appropriate for building rotamer libraries from sources where
      //! multiple of the same molecule are often present in a given unit cell, but all come from the same diffraction
      //! experiment, which is indicated by the name
      util::ShPtr< command::FlagInterface> m_SameMoleculeAndSameNameFlag;

      //! Count of duplicates
      mutable size_t m_DuplicatesCount;

      //! Count of unique molecules
      mutable size_t m_UniqueCount;

      //! output for unique
      mutable io::OFStream m_OutputUnique;

      //! output for duplicates
      mutable io::OFStream m_OutputDuplicates;

      //! input timer
      mutable util::Stopwatch m_InputTimer;

      //! output timer
      mutable util::Stopwatch m_OutputTimer;

      //! removal timer
      mutable util::Stopwatch m_RemovalTimer;

      //! fragment feed
      mutable util::ShPtr< chemistry::FragmentFeed> m_SpFeed;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      MoleculeUnique();

    public:

      //! copy constructor, only copy flags
      MoleculeUnique( const MoleculeUnique &APP) :
        m_CompareFlag( APP.m_CompareFlag),
        m_OutputFilenameFlag( APP.m_OutputFilenameFlag),
        m_OutputDuplicatesFilenameFlag( APP.m_OutputDuplicatesFilenameFlag),
        m_ConformerComparerFlag( APP.m_ConformerComparerFlag),
        m_MergeDescriptorsFlag( APP.m_MergeDescriptorsFlag),
        m_OverwriteDescriptorsFlag( APP.m_OverwriteDescriptorsFlag),
        m_SameMoleculeAndSameNameFlag( APP.m_SameMoleculeAndSameNameFlag),
        m_InputTimer( false),
        m_OutputTimer( false),
        m_RemovalTimer( false)
      {
      }

      // instantiate enumerator for application instance
      static const ApplicationType MoleculesUnique_Instance;

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      MoleculeUnique *Clone() const
      {
        return new MoleculeUnique( *this);
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

      //! @brief get the full description
      //! @return the full description
      std::string GetDescription() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief check that all the parameter choices were valid
      //! @return true if all the parameter choices were valid
      bool CheckParametersAreAcceptable() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetWebText() const;

      //! @brief update the duplicates/unique count and write the molecules out
      //! @param WAS_UNIQUE whether the molecule was unique
      //! @param FRAGMENT the fragment to write
      void PostProcess( const bool &WAS_UNIQUE) const;
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

    }; // MoleculesUnique

  } // namespace app
} // namespace bcl
#endif // BCL_APP_MOLECULE_UNIQUE_H_
