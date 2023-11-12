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
#include "math/bcl_math_limits.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include headers from the bcl - sorted alphabetically
#include "app/bcl_app_apps.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_reaction_complete.h"
#include "chemistry/bcl_chemistry_reaction_ensemble.h"
#include "chemistry/bcl_chemistry_reaction_search.h"
#include "chemistry/bcl_chemistry_reaction_worker.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ofstream.h"
#include "io/bcl_io_serialization.h"
#include "sdf/bcl_sdf_rxn_factory.h"
#include "sdf/bcl_sdf_rxn_handler.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_si_ptr_vector.h"

namespace bcl
{
  namespace app
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MatchReactions
    //! @brief Application for substituting pieces of a molecule with other fragments
    //! @details This can be used to e.g. substitute the I of a C-I bond with an N-R group
    //!
    //! @author geanesar
    //! @date 12/14/2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class MatchReactions :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      //! Reaction flag
      util::ShPtr< command::FlagInterface> m_ReactionFlag;

      //! Reactants flag
      util::ShPtr< command::FlagInterface> m_ReactantsFlag;

      //! the output filename flag
      util::ShPtr< command::FlagInterface> m_OutputFilenameFlag;

    //////////////////////
    // Helper functions //
    //////////////////////

    public:

    //////////////////////////////////
    // Construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MatchReactions();

      //! @brief a clone function
      //! @return a pointer to a copy of this class
      MatchReactions *Clone() const
      {
        return new MatchReactions( *this);
      }

      //! @brief gets the name of the class
      //! @return the name of the class
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief initializes the command object for that executable
      //! @return a ShPtr to a command object that contains this applications command line flags
      util::ShPtr< command::Command> InitializeCommand() const
      {
        // initialize command to be returned
        util::ShPtr< command::Command> sp_cmd( new command::Command());

        // insert all the flags and params
        sp_cmd->AddFlag( m_ReactionFlag);
        sp_cmd->AddFlag( m_ReactantsFlag);
        sp_cmd->AddFlag( m_OutputFilenameFlag);

        // default flags
        command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

        // return assembled Command object
        return sp_cmd;
      }

      //! @brief the main part of the application
      //! @return 0 on success
      int Main() const
      {
        std::string rxn_dir_name( m_ReactionFlag->GetFirstParameter()->GetValue());

        io::Directory rxn_dir( rxn_dir_name);
        BCL_Assert( rxn_dir.DoesExist(), "Reaction directory \"" + rxn_dir_name + "\" does not exist");
        storage::List< io::DirectoryEntry> files( rxn_dir.ListEntries());
        chemistry::ReactionEnsemble ens;
        for
        (
          storage::List< io::DirectoryEntry>::const_iterator itr_entry( files.Begin()), itr_entry_end( files.End());
          itr_entry != itr_entry_end;
          ++itr_entry
        )
        {
          if( !itr_entry->DoesExist() || itr_entry->GetType() != io::Directory::e_File)
          {
            BCL_MessageStd( "\"" + itr_entry->GetName() + "\" is not a file, skipping");
            continue;
          }
          else
          {
            BCL_MessageStd( "  Reading \"" + itr_entry->GetName() + "\"");
          }
          io::IFStream input;
          io::File::MustOpenIFStream( input, itr_entry->GetFullName());

          // TODO change this bit of code so that we don't have to call stuff from sdf namespace directly
          sdf::RXNHandler rxn_handler( input);
          if( rxn_handler.IsValid())
          {
            ens.PushBack( sdf::RXNFactory::MakeReactionComplete( rxn_handler));
          }
          else
          {
            BCL_MessageStd( "File \"" + itr_entry->GetName() + "\" did not contain RXN data");
          }
          io::File::CloseClearFStream( input);
        }

        chemistry::FragmentEnsemble dummy;
        chemistry::ReactionSearch rxn_search
        (
          dummy,
          ens,
          false
        );

        std::string reactant_file( m_ReactantsFlag->GetFirstParameter()->GetValue());
        BCL_MessageStd( "Reading molecules from " + reactant_file);
        io::IFStream in;
        io::File::MustOpenIFStream( in, reactant_file);
        chemistry::FragmentEnsemble mols( in);
        io::File::CloseClearFStream( in);
        BCL_MessageStd( "Read " + util::Format()( mols.GetSize()) + " molecules");

        io::OFStream out;
        io::File::MustOpenOFStream( out, m_OutputFilenameFlag->GetFirstParameter()->GetValue());
        for
        (
          chemistry::FragmentEnsemble::const_iterator itr_mol( mols.Begin()), itr_mol_end( mols.End());
          itr_mol != itr_mol_end;
          ++itr_mol
        )
        {
          util::SiPtrVector< const chemistry::ReactionComplete> rxns( rxn_search.FindReactions( *itr_mol));
          std::string rxn_descs;
          for
          (
            util::SiPtrVector< const chemistry::ReactionComplete>::const_iterator itr_rxn( rxns.Begin()), itr_rxn_end( rxns.End());
            itr_rxn != itr_rxn_end;
            ++itr_rxn
          )
          {
            std::string cur_desc( ( *itr_rxn)->GetDescription());
            // erase newlines from description because they'll mess up
            size_t pos( 0);
            while( ( pos = cur_desc.find( "\n", pos)) != std::string::npos)
            {
              cur_desc.erase( pos, 1);
            }
            rxn_descs += cur_desc + '\n';
          }
          if( !rxns.IsEmpty())
          {
            chemistry::FragmentComplete m( *itr_mol);
            m.StoreProperty( "MatchedReactions", rxn_descs);
            m.WriteMDL( out);
          }
        }
        io::File::CloseClearFStream( out);
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

      static const ApplicationType MatchReactions_Instance;

    }; // MatchReactions

      //! @brief default constructor
    MatchReactions::MatchReactions() :
      m_ReactionFlag
      (
        new command::FlagStatic
        (
          "reactions",
          "the reactions to use (rxn format)",
          command::Parameter
          (
            "filename",
            "file containing a reactions that will be screened",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_ReactantsFlag
      (
        new command::FlagDynamic
        (
          "reactants", "sdf file containing the reactants for reaction",
          command::Parameter
          (
            "filename",
            "file containing a substructure to match with reactants in the reaction",
            command::ParameterCheckFileExistence()
          ),
          0,
          math::GetHighestBoundedValue< size_t>()
        )
      ),
      m_OutputFilenameFlag
      (
        new command::FlagStatic
        (
          "output", "flag for the output file name",
          command::Parameter
          (
            "filename",
            "file containing modified molecules"
          )
        )
      )
    {
    }

    const ApplicationType MatchReactions::MatchReactions_Instance
    (
      GetAppGroups().AddAppToGroup( new MatchReactions(), GetAppGroups().e_Molecule)
    );

  } // namespace app
} // namespace bcl
