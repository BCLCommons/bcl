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
#include "fold/bcl_fold_protocol_dock.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_locator_domain_specified.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "coord/bcl_coord_move_transform_random.h"
#include "fold/bcl_fold_default_flags.h"
#include "fold/bcl_fold_mutate_protein_model_domain.h"
#include "fold/bcl_fold_mutate_tree.h"
#include "fold/bcl_fold_mutates.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_mutate_move_wrapper.h"
#include "math/bcl_math_mutate_result.h"
#include "pdb/bcl_pdb_factory.h"
#include "util/bcl_util_string_replacement.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProtocolDock::ProtocolDock()
    {
    }

    //! @brief Clone function
    //! @return pointer to new ProtocolDock
    ProtocolDock *ProtocolDock::Clone() const
    {
      return new ProtocolDock( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief access to only instance
    //! @return reference to only instance
    ProtocolDock &ProtocolDock::GetInstance()
    {
      static ProtocolDock s_protocol_instance;
      return s_protocol_instance;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProtocolDock::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ProtocolDock::GetAlias() const
    {
      static const std::string s_name( "ProtocolDock");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProtocolDock::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "protocol to be applied to dock domains");

      return serializer;
    }

  ///////////
  // flags //
  ///////////

    //! @brief returns all flags that are specialized for this protocol
    //! @return all flags that are specialized for this protocol
    const util::ShPtrVector< command::FlagInterface> &ProtocolDock::GetAllFlags() const
    {
      // initialize static ShPtrVector of FlagInterfaces to form the comment line
      static util::ShPtrVector< command::FlagInterface> s_all_flags_vector;

      // if the flag vector is initialize for the first time
      if( s_all_flags_vector.IsEmpty())
      {
        // insert all the flags in the vector
        s_all_flags_vector.PushBack( GetFlagDomainSpecification());
        s_all_flags_vector.PushBack( GetFlagPrintDomainToPymolScript());
        s_all_flags_vector.PushBack( pdb::Factory::GetFlagWritePDBResID());
      }

      // end
      return s_all_flags_vector;
    }

    //! @brief return command line flag for specifying the domains
    //! @return command line flag for specifying domains
    util::ShPtr< command::FlagInterface> &ProtocolDock::GetFlagDomainSpecification()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "domain_specify",
          "\tspecify the domains which should be moved"
        )
      );

      // initialize parameters
      static util::ShPtr< command::ParameterInterface> s_filename
      (
        new command::Parameter( "domain_specification_filename", "\tfull path and name of file specifying domain", "domain.txt")
      );

      // if the flag is initialized for the first time
      if( s_flag->GetParameterList().IsEmpty())
      {
        util::ShPtr< command::FlagStatic> flag( s_flag);
        // insert parameters
        flag->PushBack( s_filename);
      }

      // end
      return s_flag;
    }

    //! @brief return command line flag for specifying the domains
    //! @return command line flag for specifying domains
    util::ShPtr< command::FlagInterface> &ProtocolDock::GetFlagPrintDomainToPymolScript()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "print_domain",
          "\tprint domain to pymol script so that it can be selected on a structure"
        )
      );

      // initialize parameters
      static util::ShPtr< command::ParameterInterface> s_filename
      (
        new command::Parameter( "pymol_script_filename", "\tfull path and name of file to show domain", "domain.pml")
      );

      // if the flag is initialized for the first time
      if( s_flag->GetParameterList().IsEmpty())
      {
        util::ShPtr< command::FlagStatic> flag( s_flag);
        // insert parameters
        flag->PushBack( s_filename);
      }

      // end
      return s_flag;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief modifies the start model
    //! @param START_MODEL Start model to be modified
    void ProtocolDock::ModifyStartModel( assemble::ProteinModel &START_MODEL) const
    {
      const util::ShPtrList< assemble::LocatorSSE> &remove_sses( GetDomainMutates().Second());

      // iterate through the locators of sses to remove the located sses
      for
      (
        util::ShPtrList< assemble::LocatorSSE>::const_iterator
          locator_itr( remove_sses.Begin()), locator_itr_end( remove_sses.End());
        locator_itr != locator_itr_end;
        ++locator_itr
      )
      {
        util::SiPtr< const assemble::SSE> sse( ( *locator_itr)->Locate( START_MODEL));

        // if sse could not be located
        if( !sse.IsDefined())
        {
          BCL_MessageDbg( "could not locate sse " + ( *locator_itr)->GetSSEIDString());
          continue;
        }

        // remove sse
        BCL_Assert( START_MODEL.Remove( *sse), "could not remove sse " + sse->GetIdentification());
      }
    }

    //! @brief initialize the scores and add them to Scores enumerator
    void ProtocolDock::InitializeScores()
    {
    }

    //! @brief modify the score weight set
    //! @param SCORE_WEIGHT_SET Score weight set
    void ProtocolDock::ModifyScoreWeightSet( ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
    }

    //! @brief modify the terminate object
    //! @param CRITERION which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolDock::ModifyCriterion
    (
      opti::CriterionCombine< assemble::ProteinModel, double> &CRITERION,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief modify the printer object
    //! @param PRINTER which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolDock::ModifyPrinter
    (
      mc::PrinterCombined< assemble::ProteinModel, double> &PRINTER,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief modify the pdb factory object
    //! @param FACTORY pdb factory to be modified
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolDock::ModifyFactory
    (
      util::ShPtr< pdb::Factory> &FACTORY,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief initialize the mutates and add them to Mutates enumerator
    void ProtocolDock::InitializeMutates()
    {
      const util::ShPtrList< math::MutateInterface< assemble::ProteinModel> > &mutates( GetDomainMutates().First());

      // add the mutates
      for
      (
        util::ShPtrList< math::MutateInterface< assemble::ProteinModel> >::const_iterator
          mutate_itr( mutates.Begin()), mutate_itr_end( mutates.End());
          mutate_itr != mutate_itr_end;
        ++mutate_itr
      )
      {
        BCL_MessageDbg( "adding mutate for " + ( *mutate_itr)->GetScheme());
        GetMutates().AddMutate( *mutate_itr);
      }
    }

    //! @brief modify the mutate tree used
    //! @param MUTATE_TREE MutateTree to be modified
    void ProtocolDock::ModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
      // reset all the probabilities
      MUTATE_TREE.Reset();
      MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_Domain,    1.0);

      const util::ShPtrList< math::MutateInterface< assemble::ProteinModel> > &mutates( GetDomainMutates().First());

      // add the mutates
      for
      (
        util::ShPtrList< math::MutateInterface< assemble::ProteinModel> >::const_iterator
          mutate_itr( mutates.Begin()), mutate_itr_end( mutates.End());
          mutate_itr != mutate_itr_end;
        ++mutate_itr
      )
      {
        BCL_MessageDbg( "setting mutate probability for " + ( *mutate_itr)->GetScheme());
        MUTATE_TREE.SetMutateProbability( ( *mutate_itr)->GetScheme(),   1.0);
      }
    }

    //! @brief get the mutate tree associated with this protocol
    //! @return the mutate tree associated with this protocol
    util::ShPtr< MutateTree> ProtocolDock::GetMutateTree() const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( new MutateTree());
      ModifyMutateTree( *sp_mutate_tree);
      return sp_mutate_tree;
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void ProtocolDock::MergeAndModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }

  ////////////
  // readme //
  ////////////

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ProtocolDock::GetDescription() const
    {
      // initialize string to store the description
      static const std::string s_description
      (
        "protocol to be applied to dock domains"
      );

      // end
      return s_description;
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ProtocolDock::GetReadMe() const
    {
      // TODO: add readme for this protocol
      // initialize string to store the readme information
      static const std::string s_readme
      (
        "readme for dock protocol"
      );

      // end
      return s_readme;
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProtocolDock::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProtocolDock::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief second list is for sses that should be removed from starting model to minimize clashes during docking
    const storage::Pair
    <
      util::ShPtrList< math::MutateInterface< assemble::ProteinModel> >,
      util::ShPtrList< assemble::LocatorSSE>
    > &ProtocolDock::GetDomainMutates()
    {
      static storage::Pair
      <
        util::ShPtrList< math::MutateInterface< assemble::ProteinModel> >,
        util::ShPtrList< assemble::LocatorSSE>
      > mutates_removers;

      // initialize mutates if it is empty
      if( mutates_removers.First().IsEmpty())
      {
        io::OFStream pml_write;
        if( GetFlagPrintDomainToPymolScript()->GetFlag())
        {
          const std::string &pymol_domain_filename( GetFlagPrintDomainToPymolScript()->GetFirstParameter()->GetValue());

          io::File::MustOpenOFStream( pml_write, pymol_domain_filename);
        }
        static const std::string begin_domain_specifier( "DomainSpecifier");
        static const std::string end_domain_specifier( "DomainSpecifierEnd");
        const std::string &domain_filename( GetFlagDomainSpecification()->GetFirstParameter()->GetValue());
        io::IFStream read;
        io::File::MustOpenIFStream( read, domain_filename);

        std::string current_line;
        std::getline( read, current_line);
        current_line = util::TrimString( current_line);

        size_t domain_identifier( 0);

        while( !read.eof() && !current_line.empty())
        {
          BCL_Assert
          (
            current_line == begin_domain_specifier, "first line should be " + begin_domain_specifier +
            " but is " + current_line
          );
          // max translation line formatted "translate min = 1.0 max = 3.0" (space around "=" is optional)
          std::getline( read, current_line);
          BCL_MessageDbg( "current line is |" + current_line + "|");
          current_line = util::TrimString( current_line);
          util::StringReplacement replacer( util::StringReplacement::e_Any, "=", " ");
          replacer.ReplaceEachIn( current_line);
          // now that the ='s are spaces, can split based on spaces and it doesn't matter if spaces surrounded the ='s or not
          storage::Vector< std::string> split_line( util::SplitString( current_line));
          // make sure line has correct size
          BCL_Assert
          (
            split_line.GetSize() == 5, "split line should have 5 entries but is " + util::Format()( split_line) +
            "\nmake sure your line has format \"translate min = 1.0 max = 3.0\""
          );
          // get the min and max translation allowed for this domain
          const double min_translation( util::ConvertStringToNumericalValue< double>( split_line( 2)));
          const double max_translation( util::ConvertStringToNumericalValue< double>( split_line( 4)));
          BCL_MessageDbg
          (
            "translation is min " + util::Format()( min_translation) + " max " + util::Format()( max_translation)
          );

          // max rotation line formatted "rotate min = 1.0 max = 3.0" (space around "=" is optional)
          std::getline( read, current_line);
          current_line = util::TrimString( current_line);
          replacer.ReplaceEachIn( current_line);
          // now that the ='s are spaces, can split based on spaces and it doesn't matter if spaces surrounded the ='s or not
          split_line = util::SplitString( current_line);
          // make sure line has correct size
          BCL_Assert
          (
            split_line.GetSize() == 5, "split line should have 5 entries but is " + util::Format()( split_line) +
            "\nmake sure your line has format \"rotate min = 3.0 max = 3.14\""
          );
          const double min_rotation( util::ConvertStringToNumericalValue< double>( split_line( 2)));
          const double max_rotation( util::ConvertStringToNumericalValue< double>( split_line( 4)));
          BCL_MessageDbg
          (
            "rotation is min " + util::Format()( min_rotation) + " max " + util::Format()( max_rotation)
          );

          util::ShPtrList< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> > sse_locators;
          std::string current;
          read >> current;
          BCL_MessageDbg( "current is |" + current + "|");
          bool remove( false);
          if( current == "REMOVE")
          {
            remove = true;
          }
          while( current != end_domain_specifier)
          {
            char chain_id;
            io::Serialize::Read( chain_id, read);
            int seq_id_start, seq_id_end;
            io::Serialize::Read( seq_id_start, read);
            io::Serialize::Read( seq_id_end, read);
            util::ShPtr< assemble::LocatorSSE> sse_locator
            (
              new assemble::LocatorSSE( chain_id, seq_id_start, seq_id_end, DefaultFlags::GetFlagPDBIDNumbering()->GetFlag())
            );
            if( remove)
            {
              mutates_removers.Second().PushBack( sse_locator);
            }
            else
            {
              sse_locators.PushBack( sse_locator);
            }
            BCL_MessageDbg
            (
              "remove is " + util::Format()( remove) + " domain " +
              util::Format()( domain_identifier) + " inserted locator " + util::Format()( chain_id) + " " +
              util::Format()( seq_id_start) + " " + util::Format()( seq_id_end)
            );

            read >> current;
            BCL_MessageDbg( "current is |" + current + "|");
            if( current == "REMOVE")
            {
              remove = true;
            }
            else
            {
              remove = false;
            }
          } //< end read all sses involved in this domain
          BCL_MessageDbg( "number of sses in this domain " + util::Format()( sse_locators.GetSize()));
          // locator for the domain
          const assemble::LocatorDomainSpecified domain_locator( sse_locators);

          // method to transform the domain
          const coord::MoveTransformRandom domain_move( min_translation, max_translation, min_rotation, max_rotation);

          const math::MutateMoveWrapper< assemble::Domain> domain_mutate( domain_move, false);

          const std::string scheme( "domain_transform_random_" + util::Format()( domain_identifier));

          if( GetFlagPrintDomainToPymolScript()->GetFlag())
          {
            domain_locator.WritePymolDomainFile( pml_write, scheme);
          }

          // mutate for the protein model
          util::ShPtr< math::MutateInterface< assemble::ProteinModel> > model_mutate
          (
            new MutateProteinModelDomain( domain_locator, domain_mutate, scheme)
          );

          mutates_removers.First().PushBack( model_mutate);
          ++domain_identifier;
          std::getline( read, current_line); //< get line twice since >> was used previously
          std::getline( read, current_line); //< so first getline call just returns rest of the line >> was called on
          current_line = util::TrimString( current_line);
          BCL_MessageDbg( "end current_line is |" + current_line + "|");
        } //< end read all domains

        BCL_MessageDbg( "number domains specified " + util::Format()( mutates_removers.First().GetSize()));
        BCL_MessageDbg( "number sses to be removed " + util::Format()( mutates_removers.Second().GetSize()));
      }

      // return the mutates
      return mutates_removers;
    }

  } // namespace fold
} // namespace bcl
