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

#ifndef BCL_APP_ALIGN_TO_SCAFFOLD_H_
#define BCL_APP_ALIGN_TO_SCAFFOLD_H_

// includes from bcl - sorted alphabetically
#include "app/bcl_app_apps.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_parameter_check_file_existence.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AlignToScaffold
    //! @brief Aligns ensemble of molecules to a common scaffold
    //! @details Takes an ensemble of molecules and a scaffold sdf and aligns all of the molecules in the ensemble
    //!          onto the scaffold. Outputs ensemble in specified filename.
    //!
    //! @see @link example_app_align_to_scaffold.cpp @endlink
    //! @author kothiwsk, sliwosgr
    //! @date 03/07/2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AlignToScaffold :
      public InterfaceRelease
    {
    public:

      // static member of this class
      static const ApplicationType AlignToScaffold_Instance;

    private:

    //////////
    // data //
    //////////

      //! scaffold sdf
      util::ShPtr< command::ParameterInterface> m_InputScaffold;

      //! small molecule ensemble sdf
      util::ShPtr< command::ParameterInterface> m_EnsembleFileName;

      //! output sdf file
      util::ShPtr< command::ParameterInterface> m_OutputFilename;

      //! flag to take list of atoms that need to be aligned
      util::ShPtr< command::FlagInterface> m_AlignScaffoldAtoms;

      //! flag to take list of atoms that need to be aligned
      util::ShPtr< command::FlagInterface> m_AlignEnsembleAtoms;

      //! flag to take list of atoms that need to be aligned
      util::ShPtr< command::FlagInterface> m_AlignRigid;

      //! Scheme used for comparing whether two bonds are equivalent
      util::ShPtr< command::FlagInterface> m_BondComparisonType;

      //! Scheme used for comparing whether two atoms are equivalent
      util::ShPtr< command::FlagInterface> m_AtomComparisonType;

      //! Graph solution type for substructure comparison
      util::ShPtr< command::FlagInterface> m_SolutionTypeFlag;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      AlignToScaffold();

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      AlignToScaffold *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

      // Find the isomorphism (overlap) between the scaffold and small molecule. Will be called for each molecule in the
      // input ensemble.
      storage::Map< size_t, size_t> FindIsomorphism
      (
        const chemistry::FragmentComplete &MOLECULE,
        graph::ConstGraph< size_t, size_t> &SCAFFOLD_GRAPH
      ) const;

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

    }; // AlignToScaffold

  } // namespace app
} // namespace bcl

#endif // BCL_APP_ALIGN_TO_SCAFFOLD_H_
