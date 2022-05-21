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

// include headers from the bcl - sorted alphabetically
#include "app/bcl_app_apps.h"
#include "app/bcl_app_groups.h"
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "chemistry/bcl_chemistry_collector_valence.h"
#include "chemistry/bcl_chemistry_configurational_bond_types.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_merge_fragment_complete.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "graph/bcl_graph_connectivity.h"
#include "io/bcl_io_file.h"

namespace bcl
{
  namespace app
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AddFragments
    //! @brief Application for analyzing fragments of molecules_in which are derived from a common scaffold
    //!
    //! @author geanesar
    //! @date 05/12/2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // TODO: this really should just be another FragmentSplitInterface derived class so that we can just call it from MoleculeSplit

    class AddFragments :
      public Interface
    {

    public:

    //////////
    // data //
    //////////

      util::ShPtr< command::FlagInterface> m_FragmentsToAddFlag;

      util::ShPtr< command::FlagInterface> m_OutputFilenameFlag;

      util::ShPtr< command::FlagInterface> m_InputFilenamesFlag;

    //////////////////////////////////
    // Construction and destruction //
    //////////////////////////////////

      AddFragments();

      AddFragments *Clone() const
      {
        return new AddFragments( *this);
      }

      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const
      {
        // initialize command to be returned
        util::ShPtr< command::Command> sp_cmd( new command::Command());

        // insert all the flags and params
        sp_cmd->AddFlag( m_InputFilenamesFlag);
        sp_cmd->AddFlag( m_OutputFilenameFlag);
        sp_cmd->AddFlag( m_FragmentsToAddFlag);

        // default flags
        command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

        // return assembled Command object
        return sp_cmd;
      }

      int Main() const
      {

        // Read in the molecules
        const storage::Vector< std::string> filenames( m_InputFilenamesFlag->GetStringList());
        chemistry::FragmentEnsemble molecules;

        io::IFStream molecules_in;

        for
        (
          storage::Vector< std::string>::const_iterator itr( filenames.Begin()), itr_end( filenames.End());
          itr != itr_end;
          ++itr
        )
        {
          io::File::MustOpenIFStream( molecules_in, *itr);
          molecules.ReadMoreFromMdl( molecules_in, sdf::e_Maintain);
          io::File::CloseClearFStream( molecules_in);
        }

        io::File::MustOpenIFStream( molecules_in, m_FragmentsToAddFlag->GetFirstParameter()->GetValue());
        chemistry::FragmentEnsemble fragments_add;
        fragments_add.ReadMoreFromMdl( molecules_in, sdf::e_Maintain);

        io::File::CloseClearFStream( molecules_in);

        if( fragments_add.GetSize() == 0 || molecules.GetSize() == 0)
        {
          BCL_Exit( "Not enough molecules", -1);
        }

        io::OFStream output;
        io::File::MustOpenOFStream( output, m_OutputFilenameFlag->GetFirstParameter()->GetValue());

        chemistry::CollectorValence valence_collector;

        for
        (
          storage::List< chemistry::FragmentComplete>::iterator itr_frag( fragments_add.Begin()), itr_frag_end( fragments_add.End());
          itr_frag != itr_frag_end;
          ++itr_frag
        )
        {
          // Find the open valence
          util::SiPtrList< const chemistry::AtomConformationalInterface> open_valences( valence_collector.Collect( *itr_frag));
          if( open_valences.GetSize() != 1)
          {
            BCL_Exit( "Something is going wrong, not the right number of open valences", -1);
          }

          size_t frag_index( itr_frag->GetAtomIndex( **open_valences.Begin()));

          for
          (
              storage::List< chemistry::FragmentComplete>::iterator itr_mol( molecules.Begin()), itr_mol_end( molecules.End());
              itr_mol != itr_mol_end;
              itr_mol++
          )
          {
            util::SiPtrList< const chemistry::AtomConformationalInterface> open_mol_valences( valence_collector.Collect( *itr_mol));

            if( open_valences.GetSize() != 1)
            {
              BCL_Exit( "Something is going wrong, not the right number of open valences on mol", -1);
            }

            if( open_mol_valences.Begin() == open_mol_valences.End())
            {
              BCL_MessageStd( " Continuing");
              continue;
            }
            size_t mol_index( itr_mol->GetAtomIndex( **open_mol_valences.Begin()));

            storage::Pair< bool, chemistry::FragmentComplete> new_pair
            (
              chemistry::MergeFragmentComplete::MergeFragments
              (
                *itr_frag,
                *itr_mol,
                chemistry::GetConfigurationalBondTypes().e_ConjugatedSingleBond,
                storage::Pair< size_t, size_t>( frag_index, mol_index)
              )
            );

            new_pair.Second().WriteMDL( output);

          }
        }
        io::File::CloseClearFStream( output);

        return 0;
      } // Main()

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

    private:

      static const ApplicationType AddFragments_Instance;

    }; // AddFragments

      //! @brief standard constructor
    AddFragments::AddFragments() :
      m_FragmentsToAddFlag
      (
        new command::FlagStatic
        (
          "fragments_to_add", "filename for fragments to replace removed fragments with",
          command::Parameter
          (
            "filename",
            "filename for fragments to replace removed fragments with",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_OutputFilenameFlag
      (
        new command::FlagStatic
        (
          "output_filename", "flag selecting the output file name",
          command::Parameter
          (
            "filename", "filename for output sdf"
          )
        )
      ),
      m_InputFilenamesFlag
      (
        new command::FlagDynamic
        (
          "input_filenames",
          "filenames for input sdf",
          command::Parameter
          (
            "filenames of derived structure sdfs",
            "name of files containing derived structures",
            command::ParameterCheckFileExistence()
          ),
          1,
          21
        )
      )
    {
    }

    const ApplicationType AddFragments::AddFragments_Instance
    (
      GetAppGroups().AddAppToGroup( new AddFragments(), GetAppGroups().e_Molecule)
    );

  } // namespace app
} // namespace bcl
