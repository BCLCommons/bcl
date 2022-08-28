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
#include <util/bcl_util_undefined.h>
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "bcl_app_molecule_minimize.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_rdkit_mol_utils.h"

// external includes - sorted alphabetically
#include "GraphMol/RWMol.h"
#include "GraphMol/ForceFieldHelpers/MMFF/MMFF.h"
#include "ForceField/ForceField.h"

namespace bcl
{
  namespace app
  {

    const ApplicationType MoleculeMinimize::MoleculeMinimize_Instance
    (
      GetAppGroups().AddAppToGroup( new MoleculeMinimize(), GetAppGroups().e_Molecule)
    );

  ///////////////////////////////////
  // construction and destruction //
  ///////////////////////////////////

    //! @brief standard constructor
    MoleculeMinimize::MoleculeMinimize() :
      m_OutputFilenameBase
      (
        new command::FlagStatic
        (
          "output",
          "base name for output",
          command::Parameter( "output", "base name for output")
        )
      )
    {
    }

    //! @brief Clone function
    //! @return pointer to new MoleculeMinimize
    MoleculeMinimize *MoleculeMinimize::Clone() const
    {
      return new MoleculeMinimize( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MoleculeMinimize::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> MoleculeMinimize::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

    ////////////////////
    // common options //
    ////////////////////

      // add flags for input
      chemistry::FragmentFeed::AddFlags( *sp_cmd);

      // Output filename base
      sp_cmd->AddFlag( m_OutputFilenameBase);

      // molecule reading prefs
      sdf::AddMoleculeIOPrefFlags( *sp_cmd);

    ///////////////////
    // default flags //
    ///////////////////

      // default flags are unnecessary for this application, but message level and read-me are useful
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags
      (
        *sp_cmd,
        storage::Set< command::FlagTypeEnum>( command::e_Pthread)
      );

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string MoleculeMinimize::GetDescription() const
    {
      return "Performs molecular mechanics energy-based geometry optimization on molecules";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &MoleculeMinimize::GetReadMe() const
    {
      static std::string s_read_me =
        "MoleculeMinimize performs energy minimization of molecular geometries.";
      return s_read_me;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int MoleculeMinimize::Main() const
    {
      io::OFStream output;
      io::File::MustOpenOFStream
      (
        output,
        m_OutputFilenameBase->GetFirstParameter()->GetValue() + util::Format()( ".sdf"),
        std::ios::app
      );
      size_t feed_index( 0);
      std::shared_ptr< ::RDKit::RWMol> rdkit_mol = nullptr;
      std::shared_ptr< chemistry::FragmentComplete> min_mol = nullptr;
      for( chemistry::FragmentFeed feed; feed.NotAtEnd(); ++feed, ++feed_index)
      {
        BCL_MessageStd("Feed index: " + util::Format()( feed_index));
        rdkit_mol = chemistry::RdkitMolUtils::FragmentCompleteToRDKitRWMol( *feed);
        BCL_MessageStd("Converted from BCL to RDKit");
        ::RDKit::MMFF::MMFFOptimizeMolecule( *rdkit_mol , 1000 , "MMFF94s" );
        min_mol = chemistry::RdkitMolUtils::RDKitRWMolToFragmentComplete(*rdkit_mol);
        BCL_MessageStd("Converted from RDKit to BCL");
        BCL_MessageStd("Output final molecule");
        BCL_Debug( min_mol->GetSize());
        BCL_Debug( min_mol != nullptr);
        min_mol->WriteMDL( output);
        BCL_Debug( feed.NotAtEnd());
      }
      io::File::CloseClearFStream( output);
      BCL_MessageStd("Done!");
      return 0;
    }

  } // namespace app
} // namespace bcl

