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
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "bcl_app_make_chimeric_molecule.h"

// includes from bcl - sorted alphabetically
#include "app/bcl_app_apps.h"
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "chemistry/bcl_chemistry_fragment_split_interface.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "io/bcl_io_file.h"
#include "storage/bcl_storage_template_instantiations.h"

namespace bcl
{
  namespace app
  {
    const ApplicationType MakeChimericMolecule::MakeChimericMolecule_Instance
    (
      GetAppGroups().AddAppToGroup( new MakeChimericMolecule(), GetAppGroups().e_Molecule)
    );

    //! @brief standard constructor
    MakeChimericMolecule::MakeChimericMolecule() :
      m_RecenterFlag
      (
        new command::FlagStatic
        (
          "recenter",
          "translate molecules such that the middle of the molecule is at the origin"
        )
      ),
      m_OutputFilenameFlag
      (
        new command::FlagStatic
        (
          "output",
          "file to write split molecules into",
          command::Parameter
          (
            "output",
            "file to write split molecules into"
          )
        )
      )
    {
    }

    //! copy constructor; skips i/o streams
    MakeChimericMolecule::MakeChimericMolecule( const MakeChimericMolecule &PARENT) :
          m_RecenterFlag( PARENT.m_RecenterFlag),
          m_OutputFilenameFlag( PARENT.m_OutputFilenameFlag)
    {
    }

    //! @brief Clone function
    //! @return pointer to new FoldProtein
    MakeChimericMolecule *MakeChimericMolecule::Clone() const
    {
      return new MakeChimericMolecule( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MakeChimericMolecule::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> MakeChimericMolecule::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

    ////////////////////
    // common options //
    ////////////////////

      // ensembles containing the molecules to be edited
      chemistry::FragmentFeed::AddFlags( *sp_cmd);

      // add command line options to add/remove hydrogens
      sdf::AddMoleculeIOPrefFlags( *sp_cmd);

      //! whether to recenter the molecules
      sp_cmd->AddFlag( m_RecenterFlag);

      //! output filename
      sp_cmd->AddFlag( m_OutputFilenameFlag);

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

    //! @brief the Main function
    //! @return error code - 0 for success
    int MakeChimericMolecule::Main() const
    {
      // open the output file
      io::File::MustOpenOFStream( m_OutputFile, m_OutputFilenameFlag->GetFirstParameter()->GetValue());

      // Make pseudo molecules from each ensemble
      chemistry::AtomVector< chemistry::AtomComplete> atom_vector;
      storage::Vector< sdf::BondInfo> empty_bonds( 0);

      // read in molecules
      chemistry::FragmentFeed itr_fragments;

      // combine molecules
      for( ; itr_fragments.NotAtEnd(); ++itr_fragments)
      {
        atom_vector.AddAtomsWithConnectivity( itr_fragments->GetAtomVector(), empty_bonds);
      }
      chemistry::FragmentComplete pseudo_mol_a( atom_vector, "");

      // write supermolecule
      Write( pseudo_mol_a);
      io::File::CloseClearFStream( m_OutputFile);

      // end
      return 0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MakeChimericMolecule::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &MakeChimericMolecule::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief writes a single molecule out
    //! @param MOLECULE the molecule or fragment to write
    void MakeChimericMolecule::Write( chemistry::FragmentComplete &MOLECULE) const
    {
      if( m_RecenterFlag->GetFlag())
      {
        MOLECULE.Translate( -MOLECULE.GetCenter());
      }
      MOLECULE.WriteMDL( m_OutputFile);
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &MakeChimericMolecule::GetReadMe() const
    {
      static std::string s_read_me =
        "MakeChimericMolecule combines multiple input molecules into a disconnected superstructure.\n"
        "This is a temporary hack to align to multiple molecules simultaneously.\n";
      return s_read_me;
    }

  } // namespace app
} // namespace bcl
