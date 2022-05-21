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

#ifndef BCL_APP_ALIGN_BINDING_POSES_H_
#define BCL_APP_ALIGN_BINDING_POSES_H_

// includes from bcl - sorted alphabetically
#include "app/bcl_app_apps.h"
#include "assemble/bcl_assemble_protein_model.h"
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
    //! @class AlignBindingPoses
    //! @brief Aligns the binding poses of molecules in PDBs by superimposing reference PDB structures.
    //! @details because PDB format does not support proper ligand information (bond order; connectivity)
    //!          we often need to use ligand that were curated by PDB-Bind. However, the input PDBs are rarely aligned to
    //!          one another. This application allows superimposition of PDBs and their associated ligands in separate files
    //!
    //! @author mendenjl
    //! @date Oct 17, 2017
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AlignBindingPoses :
      public Interface
    {
    public:

      // static member of this class
      static const ApplicationType AlignBindingPoses_Instance;

    private:

    //////////
    // data //
    //////////

      //! scaffold sdf
      util::ShPtr< command::ParameterInterface> m_InputMolInBindingPos;

      //! PDB structure, may have ligand in the binding pose
      util::ShPtr< command::ParameterInterface> m_TemplatePDB;

      //! output sdf file
      util::ShPtr< command::ParameterInterface> m_OutputFilename;

      //! flag to take list of atoms that need to be aligned
      util::ShPtr< command::FlagInterface> m_Molecules;

      //! flag to take list of atoms that need to be aligned
      util::ShPtr< command::FlagInterface> m_Pdbs;

      //! flag for defining the method of superimposition
      util::ShPtr< command::FlagStatic> m_SuperimposeMethodFlag;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      AlignBindingPoses();

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      AlignBindingPoses *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

      void WriteContacts( const chemistry::FragmentComplete &FRAG, const assemble::ProteinModel &MODEL) const;

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

    }; // AlignBindingPoses

  } // namespace app
} // namespace bcl

#endif // BCL_APP_ALIGN_BINDING_POSES_H_
