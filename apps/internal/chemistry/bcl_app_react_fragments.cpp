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
    //! @class ReactLegacy
    //! @brief Application for substituting pieces of a molecule with other fragments
    //! @details This can be used to e.g. substitute the I of a C-I bond with an N-R group
    //!
    //! @author geanesar
    //! @date 12/14/2014
    //!
    //! TODO: this really should just be another FragmentSplitInterface derived class so that we can just call it from MoleculeSplit
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class ReactLegacy :
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
      ReactLegacy();

      //! @brief a clone function
      //! @return a pointer to a copy of this class
      ReactLegacy *Clone() const
      {
        return new ReactLegacy( *this);
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
        storage::Vector< chemistry::FragmentEnsemble> reactant_vector;
        reactant_vector.AllocateMemory( reactant_filenames.GetSize());
        for( size_t fn( 0), end_fn( reactant_filenames.GetSize()); fn < end_fn; ++fn)
        {
          io::IFStream input;
          io::File::MustOpenIFStream( input, reactant_filenames( fn));
          reactant_vector.PushBack( chemistry::FragmentEnsemble( input));
          io::File::CloseClearFStream( input);
        }

        chemistry::ReactionWorker worker;
        io::OFStream out;
        io::File::MustOpenOFStream( out, m_OutputFilenameFlag->GetFirstParameter()->GetValue());
        for
        (
          chemistry::ReactionEnsemble::const_iterator itr_rxn( rxn_ensemble.Begin()), itr_rxn_end( rxn_ensemble.End());
          itr_rxn != itr_rxn_end;
          ++itr_rxn
        )
        {

          size_t product_no( 0);

          if( reactant_filenames.GetSize() != itr_rxn->GetNumberReactants())
          {
            BCL_MessageStd( "Cannot execute reaction " + itr_rxn->GetDescription());
            continue;
          }

          if( itr_rxn->GetNumberReactants() == 1)
          {
            size_t mol_no( 0);
            for
            (
              chemistry::FragmentEnsemble::const_iterator itr_r1( reactant_vector( 0).Begin()), itr_r1_end( reactant_vector( 0).End());
              itr_r1 != itr_r1_end;
              ++itr_r1, ++mol_no
            )
            {
              if( !( mol_no % 1000))
              {
                util::GetLogger().LogStatus( "Reacted " + util::Format()( mol_no) + " molecules");
              }
              chemistry::FragmentEnsemble reactants;
              reactants.PushBack( *itr_r1);
              chemistry::FragmentEnsemble result( worker.React( *itr_rxn, reactants));
              if( !result.IsEmpty())
              {
                result.WriteMDL( out);
              }
            }
          }
          else if( itr_rxn->GetNumberReactants() == 2)
          {
            size_t n1( reactant_vector( 0).GetSize());
            size_t n2( reactant_vector( 1).GetSize());
            size_t total( n1 * n2);
            for
            (
              chemistry::FragmentEnsemble::const_iterator itr_r1( reactant_vector( 0).Begin()), itr_r1_end( reactant_vector( 0).End());
              itr_r1 != itr_r1_end;
              ++itr_r1
            )
            {
              for
              (
                chemistry::FragmentEnsemble::const_iterator itr_r2( reactant_vector( 1).Begin()), itr_r2_end( reactant_vector( 1).End());
                itr_r2 != itr_r2_end;
                ++itr_r2, ++product_no
              )
              {
                if( ( product_no % 1000) == 0)
                {
                  util::GetLogger().LogStatus( "Reacted " + util::Format()( product_no) + " pairs out of " + util::Format()( total));
                }
                chemistry::FragmentEnsemble reactants;
                reactants.PushBack( *itr_r1);
                reactants.PushBack( *itr_r2);
                chemistry::FragmentEnsemble res( worker.React( *itr_rxn, reactants));
                if( !res.IsEmpty())
                {
                  res.WriteMDL( out);
                }
              }
            }
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

      static const ApplicationType React_Instance;

    }; // React

      //! @brief default constructor
    ReactLegacy::ReactLegacy() :
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

    const ApplicationType ReactLegacy::React_Instance
    (
      GetAppGroups().AddAppToGroup( new ReactLegacy(), GetAppGroups().e_Molecule)
    );

  } // namespace app
} // namespace bcl
