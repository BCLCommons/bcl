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
#include "fold/bcl_fold_protocols.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "fold/bcl_fold_protocol_assembly.h"
#include "fold/bcl_fold_protocol_create.h"
#include "fold/bcl_fold_protocol_default.h"
#include "fold/bcl_fold_protocol_dock.h"
#include "fold/bcl_fold_protocol_em.h"
#include "fold/bcl_fold_protocol_ensemble.h"
#include "fold/bcl_fold_protocol_ensemble_filter.h"
#include "fold/bcl_fold_protocol_ensemble_replicate_conformation.h"
#include "fold/bcl_fold_protocol_ensemble_switch_conformation.h"
#include "fold/bcl_fold_protocol_loop_close.h"
#include "fold/bcl_fold_protocol_loop_coordinate_add.h"
#include "fold/bcl_fold_protocol_membrane.h"
#include "fold/bcl_fold_protocol_multimer.h"
#include "fold/bcl_fold_protocol_refinement.h"
#include "fold/bcl_fold_protocol_restraint.h"
#include "fold/bcl_fold_protocol_template.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Protocols::Protocols() :
      e_Default(           AddEnum( "Default",           util::SiPtr< ProtocolInterface>( ProtocolDefault::GetInstance()))),
      e_Create(            AddEnum( "Create",            util::SiPtr< ProtocolInterface>( ProtocolCreate::GetInstance()))),
      e_Assembly(          AddEnum( "Assembly",          util::SiPtr< ProtocolInterface>( ProtocolAssembly::GetInstance()))),
      e_Refinement(        AddEnum( "Refinement",        util::SiPtr< ProtocolInterface>( ProtocolRefinement::GetInstance()))),
      e_Membrane(          AddEnum( "Membrane",          util::SiPtr< ProtocolInterface>( ProtocolMembrane::GetInstance()))),
      e_Restraint(         AddEnum( "Restraint",         util::SiPtr< ProtocolInterface>( ProtocolRestraint::GetInstance()))),
      e_EM(                AddEnum( "EM",                util::SiPtr< ProtocolInterface>( ProtocolEM::GetInstance()))),
      e_Multimer(          AddEnum( "Multimer",          util::SiPtr< ProtocolInterface>( ProtocolMultimer::GetInstance()))),
      e_Template(          AddEnum( "Template",          util::SiPtr< ProtocolInterface>( ProtocolTemplate::GetInstance()))),
      e_LoopCoordinateAdd( AddEnum( "LoopCoordinateAdd", util::SiPtr< ProtocolInterface>( ProtocolLoopCoordinateAdd::GetInstance()))),
      e_LoopClose(         AddEnum( "LoopClose",         util::SiPtr< ProtocolInterface>( ProtocolLoopClose::GetInstance()))),
      e_Dock(              AddEnum( "Dock",              util::SiPtr< ProtocolInterface>( ProtocolDock::GetInstance()))),
      e_Ensemble(          AddEnum( "Ensemble",          util::SiPtr< ProtocolInterface>( ProtocolEnsemble::GetInstance()))),
      e_EnsembleSwitchConformation   ( AddEnum( "EnsembleSwitchConformation",    util::SiPtr< ProtocolInterface>( ProtocolEnsembleSwitchConformation::GetInstance()))),
      e_EnsembleReplicateConformation( AddEnum( "EnsembleReplicateConformation", util::SiPtr< ProtocolInterface>( ProtocolEnsembleReplicateConformation::GetInstance()))),
      e_EnsembleFilter               ( AddEnum( "EnsembleFilter",                util::SiPtr< ProtocolInterface>( ProtocolEnsembleFilter::GetInstance())))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief return command line flag for choosing the protocols to be used to run folding
    //! @return command line flag for choosing the protocols to be used to run folding
    util::ShPtr< command::FlagInterface> &Protocols::GetFlagProtocols()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagDynamic
        (
          "protocols",
          "\tFlag for choosing the protocols to be used to run folding",
          command::Parameter
          (
            "protocol",
            "\tprotocol to use",
            command::ParameterCheckEnumerate< Protocols>(),
            GetProtocols().e_Default.GetName()
          ), 0, GetProtocols().GetEnumCount()
        )
      );

      // end
      return s_flag;
    }

    //! @brief return command line flag for choosing the protocols to be used for mutates
    //! @return command line flag for choosing the protocols to be used for mutates
    util::ShPtr< command::FlagInterface> &Protocols::GetFlagMutateProtocols()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagDynamic
        (
          "mutate_protocols",
          "\tFlag for choosing the protocols to be used to for mutates",
          command::Parameter
          (
            "mutate_protocol",
            "\tmutate protocol to use",
            command::ParameterCheckEnumerate< Protocols>(),
            GetProtocols().e_Default.GetName()
          ), 0, GetProtocols().GetEnumCount()
        )
      );

      // end
      return s_flag;
    }

    //! @brief return command line flag for choosing the protocols to be used for scoring
    //! @return command line flag for choosing the protocols to be used for scoring
    util::ShPtr< command::FlagInterface> &Protocols::GetFlagScoreProtocols()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagDynamic
        (
          "score_protocols",
          "\tFlag for choosing the protocols to be used to for scores",
          command::Parameter
          (
            "score_protocol",
            "\tscore_protocol to use",
            command::ParameterCheckEnumerate< Protocols>(),
            GetProtocols().e_Default.GetName()
          ), 0, GetProtocols().GetEnumCount()
        )
      );

      // end
      return s_flag;
    }

    //! @brief return list of protocols specified to be used specified in command line
    //! @return list of protocols to be used specified in command line
    storage::List< Protocol> Protocols::GetCommandLineProtocolList()
    {
      // initialize list to hold protocols
      storage::List< Protocol> protocols_list( 1, GetProtocols().e_Default);

      // get the vector of protocols from command line flag
      storage::Vector< Protocol> protocols( GetFlagProtocols()->GetObjectList< Protocol>());

      // iterate over these
      for
      (
        storage::Vector< Protocol>::const_iterator
          protocol_itr( protocols.Begin()), protocol_itr_end( protocols.End());
        protocol_itr != protocol_itr_end; ++protocol_itr
      )
      {
        // if it's not already in the list
        if( std::find( protocols_list.Begin(), protocols_list.End(), *protocol_itr) == protocols_list.End())
        {
          // then pushback into the list
          protocols_list.PushBack( *protocol_itr);
        }
      }

      // end
      return protocols_list;
    }

    //! @brief return list of protocols to be used for defining scores specified in command line
    //! @return list of protocols to be used for defining scores specified in command line
    storage::List< Protocol> Protocols::GetCommandLineScoreProtocolList()
    {
      // initialize list to hold protocols
      storage::List< Protocol> protocols_list( GetCommandLineProtocolList());

      // get the vector of protocols from command line flag
      storage::Vector< Protocol> protocols( GetFlagScoreProtocols()->GetObjectList< Protocol>());

      // iterate over these
      for
      (
        storage::Vector< Protocol>::const_iterator
          protocol_itr( protocols.Begin()), protocol_itr_end( protocols.End());
        protocol_itr != protocol_itr_end; ++protocol_itr
      )
      {
        // if it's not already in the list
        if( std::find( protocols_list.Begin(), protocols_list.End(), *protocol_itr) == protocols_list.End())
        {
          // then pushback into the list
          protocols_list.PushBack( *protocol_itr);
        }
      }

      // end
      return protocols_list;
    }

    //! @brief return list of protocols to be used for defining mutates specified in command line
    //! @return list of protocols to be used for defining mutates specified in command line
    storage::List< Protocol> Protocols::GetCommandLineMutateProtocolList()
    {
      // initialize static list to hold protocols
      storage::List< Protocol> protocols_list( GetCommandLineProtocolList());

      // get the vector of protocols from command line flag
      storage::Vector< Protocol> protocols( GetFlagMutateProtocols()->GetObjectList< Protocol>());

      // iterate over these
      for
      (
        storage::Vector< Protocol>::const_iterator
          protocol_itr( protocols.Begin()), protocol_itr_end( protocols.End());
        protocol_itr != protocol_itr_end; ++protocol_itr
      )
      {
        // if it's not already in the list
        if( std::find( protocols_list.Begin(), protocols_list.End(), *protocol_itr) == protocols_list.End())
        {
          // then pushback into the list
          protocols_list.PushBack( *protocol_itr);
        }
      }

      // end
      return protocols_list;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Protocols::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief construct on access function for all Protocols
    //! @return reference to only instances of Protocols
    Protocols &GetProtocols()
    {
      return Protocols::GetEnums();
    }

  } // namespace fold

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< SiPtr< fold::ProtocolInterface>, fold::Protocols>;

  } // namespace util
} // namespace bcl
