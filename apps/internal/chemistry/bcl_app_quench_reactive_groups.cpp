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
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_reaction_complete.h"
#include "chemistry/bcl_chemistry_reaction_ensemble.h"
#include "chemistry/bcl_chemistry_reaction_worker.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ofstream.h"
#include "io/bcl_io_serialization.h"
#include "sdf/bcl_sdf_rxn_factory.h"
#include "sdf/bcl_sdf_rxn_handler.h"

namespace bcl
{
  namespace app
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class QuenchReactiveGroups
    //! @brief Application for substituting pieces of a molecule with other fragments
    //! @details This can be used to e.g. substitute the I of a C-I bond with an N-R group
    //!
    //! @author geanesar
    //! @date 12/14/2014
    //!
    //! TODO: this really should just be another FragmentSplitInterface derived class so that we can just call it from MoleculeSplit
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class QuenchReactiveGroups :
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
      QuenchReactiveGroups();

      //! @brief a clone function
      //! @return a pointer to a copy of this class
      QuenchReactiveGroups *Clone() const
      {
        return new QuenchReactiveGroups( *this);
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

        io::IFStream input;
        io::File::MustOpenIFStream( input, m_ReactionFlag->GetFirstParameter()->GetValue());
        chemistry::ReactionEnsemble rxn_ensemble( input);
        io::File::CloseClearFStream( input);

        storage::Vector< std::string> reactant_filenames( m_ReactantsFlag->GetStringList());
        chemistry::FragmentEnsemble mols;
        for( size_t fn( 0), end_fn( reactant_filenames.GetSize()); fn < end_fn; ++fn)
        {
          io::IFStream input;
          io::File::MustOpenIFStream( input, reactant_filenames( fn));
          mols.ReadMoreFromMdl( input);
          io::File::CloseClearFStream( input);
        }

        chemistry::ReactionWorker worker;

        size_t n_transforms( 0);
        storage::Vector< int> transform_count;
        for
        (
          chemistry::ReactionEnsemble::const_iterator itr_rxn( rxn_ensemble.Begin()), itr_rxn_end( rxn_ensemble.End());
          itr_rxn != itr_rxn_end;
          ++itr_rxn
        )
        {
          transform_count.PushBack( 0);
          if( itr_rxn->GetNumberReactants() != 1)
          {
            transform_count.LastElement() = -1;
            continue;
          }

          for
          (
            chemistry::FragmentEnsemble::iterator itr_mol( mols.Begin()), itr_mol_end( mols.End());
            itr_mol != itr_mol_end;
            ++itr_mol
          )
          {
            if( worker.MatchesReactants( *itr_mol, *itr_rxn).GetSize())
            {
              ++n_transforms;
              ++transform_count.LastElement();
              chemistry::FragmentEnsemble reactants;
              reactants.PushBack( *itr_mol);
              chemistry::FragmentEnsemble modified( worker.ExecuteReaction( *itr_rxn, reactants));
              if( modified.GetSize() > 1)
              {
                chemistry::FragmentEnsemble::const_iterator itr_mod( modified.Begin()), itr_mod_end( modified.End());
                *itr_mol = *itr_mod;
                ++itr_mod;
                for( ; itr_mod != itr_mod_end; ++itr_mod)
                {
                  mols.PushBack( *itr_mod);
                }
                itr_mol_end = mols.End();
              }
            }
          }
        }
        BCL_MessageStd( "Performed " + util::Format()( n_transforms) + " functional group transformations");
        size_t p( 0);
        for
        (
          chemistry::ReactionEnsemble::const_iterator itr_rxn( rxn_ensemble.Begin()), itr_rxn_end( rxn_ensemble.End());
          itr_rxn != itr_rxn_end;
          ++itr_rxn, ++p
        )
        {
          std::string rxn_desc( "  Reaction " + itr_rxn->GetDescription());
          for( size_t i( 11), l( rxn_desc.length()); i < l; ++i)
          {
            if( rxn_desc[ i] == '\n')
            {
              rxn_desc[ i] = ' ';
            }
          }
          if( transform_count( p) >= 0)
          {
            rxn_desc += ": " + util::Format()( transform_count( p)) + " instances";
          }
          else
          {
            rxn_desc += ": invalid reaction";
          }
          BCL_MessageStd( rxn_desc);
        }
        io::OFStream out;
        io::File::MustOpenOFStream( out, m_OutputFilenameFlag->GetFirstParameter()->GetValue());
        mols.WriteMDL( out);
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

      static const ApplicationType React_Instance;

    }; // QuenchReactiveGroups

      //! @brief default constructor
    QuenchReactiveGroups::QuenchReactiveGroups() :
      m_ReactionFlag
      (
        new command::FlagStatic
        (
          "reaction",
          "the reaction to use (rxn format)",
          command::Parameter
          (
            "filename",
            "file containing a reaction that will be used to combine reactants together",
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
          10
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

    const ApplicationType QuenchReactiveGroups::React_Instance
    (
      GetAppGroups().AddAppToGroup( new QuenchReactiveGroups(), GetAppGroups().e_Molecule)
    );

  } // namespace app
} // namespace bcl
