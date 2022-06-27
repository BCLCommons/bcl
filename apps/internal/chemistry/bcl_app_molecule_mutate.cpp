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
#include "command/bcl_command_parameter_check_ranged.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "app/bcl_app_apps.h"
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "chemistry/bcl_chemistry_fragment_mutate_interface.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "internal/chemistry/bcl_app_molecule_mutate.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_mutate_interface.h"
#include "storage/bcl_storage_template_instantiations.h"
namespace bcl
{
  namespace app
  {
    const ApplicationType MoleculeMutate::MoleculeMutate_Instance
    (
      GetAppGroups().AddAppToGroup( new MoleculeMutate(), GetAppGroups().e_Molecule)
    );

    //! @brief standard constructor
    MoleculeMutate::MoleculeMutate() :
      m_ImplementationFlag
      (
        new command::FlagDynamic
        (
          "implementation",
          "method to mutate molecules",
          command::Parameter
          (
            "mutate",
            "",
            command::ParameterCheckSerializable
            (
              util::Implementation< chemistry::FragmentMutateInterface>()
            )
          )
        )
      ),
      m_MutableAtomsFlag
      (
        new command::FlagDynamic
        (
          "mutable_atoms",
          "sequentially mutate each atom; "
          "this flag overrides the serializable atom selection of the mutate "
          "specified with 'implementation'; this flag will cause each atom specified "
          "to be mutated and saved as an independent output; if 'accumulate' is specified "
          "then these changes are aggregated; if you only desire a single mutation pulled "
          "from a pool of potential atoms, just use the 'mutable_atoms' option within "
          "the implementation",
          command::Parameter
          (
            "atom_index",
            "",
            command::ParameterCheckRanged< size_t>( 0, std::numeric_limits< size_t>::max()),
            ""
          )
        )
      ),
      m_AccumulateFlag
      (
        new command::FlagStatic
        (
          "accumulate",
          "accumulate mutations into one molecule; be careful"
        )
      ),
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
          "file to write mutated molecules into",
          command::Parameter
          (
            "output",
            "file to write mutated molecules into",
            "bcl_moleculeMutate.sdf"
          )
        )
      )
    {
    }

    //! copy constructor; skips i/o streams
    MoleculeMutate::MoleculeMutate( const MoleculeMutate &PARENT) :
          m_ImplementationFlag( PARENT.m_ImplementationFlag),
          m_AccumulateFlag( PARENT.m_AccumulateFlag),
          m_MutableAtomsFlag( PARENT.m_MutableAtomsFlag),
          m_RecenterFlag( PARENT.m_RecenterFlag),
          m_OutputFilenameFlag( PARENT.m_OutputFilenameFlag),
          m_MoleculeIndex( 0)
    {
    }

    //! @brief Clone function
    //! @return pointer to new FoldProtein
    MoleculeMutate *MoleculeMutate::Clone() const
    {
      return new MoleculeMutate( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MoleculeMutate::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &MoleculeMutate::GetReadMe() const
    {
      static std::string s_read_me =
        "MoleculeMutate applies individual mutates to molecules to transform them into new molecules.";
      return s_read_me;
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> MoleculeMutate::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

    ////////////////////
    // common options //
    ////////////////////

      // ensembles containing the molecules to be edited
      chemistry::FragmentFeed::AddFlags( *sp_cmd);

      // add command line options to add/remove hydrogens
      sdf::AddMoleculeIOPrefFlags( *sp_cmd);

      //! mutate molecular complexes into separate small molecules
      sp_cmd->AddFlag( m_ImplementationFlag);

      //! scan each specified atom for mutation
      sp_cmd->AddFlag( m_MutableAtomsFlag);

      //! make all mutations to one molecule
      sp_cmd->AddFlag( m_AccumulateFlag);

      //! whether to recenter the molecules
      sp_cmd->AddFlag( m_RecenterFlag);

      //! output filename
      sp_cmd->AddFlag( m_OutputFilenameFlag);

    ///////////////////
    // default flags //
    ///////////////////

      // add default flags
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int MoleculeMutate::Main() const
    {
      // exit if no implementation specified
      if( !m_ImplementationFlag->GetFlag())
      {
        BCL_MessageStd( "No implementation specified.");
        return 0;
      }
      auto mutates( m_ImplementationFlag->GetParameterList());
      size_t n_mutates( mutates.GetSize());

      // open the output file
      io::File::MustOpenOFStream( m_OutputFile, m_OutputFilenameFlag->GetFirstParameter()->GetValue());

      // TODO: distribute itr_fragments
      // load molecules from application interface
      chemistry::FragmentFeed itr_fragments;
      m_MoleculeIndex = 0;

      // TODO: this flag could be handled using similar logic as the implementation,
      // i.e. allow atom selection in multiple ways
      // overwrite serialized FragmentMutateInterface mutable atoms if specified here
      if( m_MutableAtomsFlag->GetFlag())
      {
        // setup the mutable atom indices
        storage::Vector< size_t> mutable_atoms;

        // convert the functionalization points to numeric values
        storage::Vector< std::string> fxnl_pts_str( m_MutableAtomsFlag->GetStringList());

        // output mutable indices to terminal
        std::string mutable_indices_message;
        for
        (
            auto itr( fxnl_pts_str.Begin()), itr_end( fxnl_pts_str.End());
            itr != itr_end;
            ++itr
        )
        {
          mutable_indices_message.append( util::Format()( *itr));
          mutable_indices_message.append(" ");
        }
        BCL_MessageStd( "Mutable atom indices: " + util::Format()( mutable_indices_message));

        // unique atom indices
        storage::Set< size_t> fxnl_pts_set;
        for( size_t i( 0), l( fxnl_pts_str.GetSize()); i < l; ++i)
        {
          size_t point;
          if( !util::TryConvertFromString( point, fxnl_pts_str( i), util::GetLogger()))
          {
            BCL_MessageStd( "Could not parse \"" + fxnl_pts_str( i) + "\" as a number");
            continue;
          }
          fxnl_pts_set.Insert( point);
        }
        mutable_atoms = storage::Vector< size_t>( fxnl_pts_set.Begin(), fxnl_pts_set.End());

        // perform a single mutate on each mutable atom for each input molecule
        // note - if you specified more than one input molecule AND mutable atom indices,
        // you better be sure that your atom indices correspond in each structure to what
        // you think they do
        // TODO: write a class that organizes atom indices based on a reference molecule

        // perform mutate
        // iterate over all input molecules
        for( ; itr_fragments.NotAtEnd(); ++itr_fragments, ++m_MoleculeIndex)
        {
          // apply multiple mutates to the same molecule
          chemistry::FragmentComplete current_frag( *itr_fragments);

          // iterate over each mutable atom
          size_t i( 0);
          for
          (
              storage::Vector< size_t>::iterator atom_itr( mutable_atoms.Begin()),
              atom_itr_end( mutable_atoms.End());
              atom_itr != atom_itr_end;
              ++atom_itr, ++i
          )
          {
            // defeats the purpose of using this flag
            BCL_Assert( mutable_atoms.GetSize(), "No mutable atoms specified!");

            // loop over mutates
            for( size_t i( 0); i < n_mutates; ++i)
            {
              // get the current mutate
              m_Mutate = mutates( i)->GetValue();

              // set the current mutable atom
              storage::Vector< size_t> mutable_atom( size_t( 1), *atom_itr);
              m_Mutate->SetMutableAtomIndices( mutable_atom);

              // perform mutation
              math::MutateResult< chemistry::FragmentComplete> mutated_object( ( *m_Mutate)( current_frag));
              if( mutated_object.GetArgument().IsDefined())
              {
                if( m_AccumulateFlag->GetFlag() && i < mutable_atoms.GetSize() - 1)
                {
                  current_frag = *( mutated_object.GetArgument());
                  continue;
                }
                chemistry::FragmentComplete fragment( *( mutated_object.GetArgument()));
                Write( fragment);
              }
            }
          }
        }
      }
      // perform a single mutate on a random mutable atom from each input molecule
      else
      {
        // perform mutate
        // iterate over all input molecules
        for( ; itr_fragments.NotAtEnd(); ++itr_fragments, ++m_MoleculeIndex)
        {
          // apply multiple mutates to the same molecule
          chemistry::FragmentComplete current_frag( *itr_fragments);

          // loop over mutates
          for( size_t i( 0); i < n_mutates; ++i)
          {
            // perform mutation
            m_Mutate.Reset();
            m_Mutate = mutates( i)->GetValue();
            math::MutateResult< chemistry::FragmentComplete> mutated_object( ( *m_Mutate)( current_frag));
            if( mutated_object.GetArgument().IsDefined())
            {
              current_frag = *( mutated_object.GetArgument());
              chemistry::FragmentComplete fragment( *( mutated_object.GetArgument()));
              Write( fragment);
            }
          }
        }
      }

      // end
      io::File::CloseClearFStream( m_OutputFile);
      return 0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MoleculeMutate::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &MoleculeMutate::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief writes a single molecule out
    //! @param MOLECULE the molecule or fragment to write
    void MoleculeMutate::Write( chemistry::FragmentComplete &MOLECULE) const
    {
      if( m_RecenterFlag->GetFlag())
      {
        MOLECULE.Translate( -MOLECULE.GetCenter());
      }
      MOLECULE.WriteMDL( m_OutputFile);
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string MoleculeMutate::GetDescription() const
    {
      return "Generate molecules via multiple different mutating schemes";
    }

  } // namespace app
} // namespace bcl
