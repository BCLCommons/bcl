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
#include "fold/bcl_fold_protocol_multimer.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_printer_protein_model_multimer.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "coord/bcl_coord_move_rotate_random.h"
#include "coord/bcl_coord_move_translate_external_axis.h"
#include "coord/bcl_coord_move_translate_random.h"
#include "fold/bcl_fold_default_flags.h"
#include "fold/bcl_fold_mutate_protein_model.h"
#include "fold/bcl_fold_mutate_protein_model_sse_swap_multimer.h"
#include "fold/bcl_fold_mutate_tree.h"
#include "fold/bcl_fold_setup.h"
#include "mc/bcl_mc_printer_combined.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_printer_quality_multimer.h"
#include "sspred/bcl_sspred_method_handler.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProtocolMultimer::ProtocolMultimer() :
      m_CyclicSubunits( util::GetUndefined< size_t>()),
      m_IsDihedral( false),
      m_Multiplier()
    {
    }

    //! @brief Clone function
    //! @return pointer to new ProtocolMultimer
    ProtocolMultimer *ProtocolMultimer::Clone() const
    {
      return new ProtocolMultimer( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief access to only instance
    //! @return reference to only instance
    ProtocolMultimer &ProtocolMultimer::GetInstance()
    {
      static ProtocolMultimer s_protocol_instance;
      return s_protocol_instance;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProtocolMultimer::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the native multimer
    //! @return the native multimer
    const util::ShPtr< assemble::ProteinModel> &ProtocolMultimer::GetNativeMultimer() const
    {
      static util::ShPtr< assemble::ProteinModel> s_native;

      // if the native is undefined
      if( !s_native.IsDefined())
      {
        // get the native
        s_native = util::ShPtr< assemble::ProteinModel>
        (
          (
            pdb::Factory().ProteinModelFromPDBFilename( GetFlagNativeMultimer()->GetFirstParameter()->GetValue()).Clone()
          )
        );

        // if the list is not empty
        storage::Set< sspred::Method> ss_pred_methods( sspred::Methods::GetCommandLineMethods());
        if( !ss_pred_methods.IsEmpty())
        {
          //try to add the ss predictions, if it fails, copy them from the parent chains
          if
          (
            !sspred::MethodHandler::ReadPredictionsForProteinModel
            (
              ss_pred_methods,
              *s_native,
              DefaultFlags::GetFlagReadSequenceDataPath()->GetParameterList()( 1)->GetValue(),
              DefaultFlags::GetFlagReadSequenceDataPath()->GetFirstParameter()->GetValue()
            )
          )
          {
            // get the chain mappings from the multiplier
            const storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> > chain_mappings
            (
              m_Multiplier.GetTransformationMatrices()
            );

            // iterate over the mappings
            for
            (
              storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> >::const_iterator
                itr( chain_mappings.Begin()), itr_end( chain_mappings.End());
              itr != itr_end; ++itr
            )
            {
              // set the predictions
              util::ShPtr< assemble::Chain> sp_target_chain( s_native->GetChain( itr->Second()));
              SetSSPrediction
              (
                *( sp_target_chain->GetSequence()),
                *( s_native->GetChain( itr->First())->GetSequence())
              );
            }
          }
        }
      }

      // end
      return s_native;
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ProtocolMultimer::GetAlias() const
    {
      static const std::string s_name( "ProtocolMultimer");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProtocolMultimer::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "protocol for folding multimeric proteins");
      serializer.AddInitializer
        (
         "cyclic subunits",
         "number of cyclic subunits",
         io::Serialization::GetAgent( &m_CyclicSubunits)
         );
      serializer.AddInitializer
        (
         "is dihedral",
         "true if dihedral symmetry, false if cyclic symmetry",
         io::Serialization::GetAgent( &m_IsDihedral)
         );
      serializer.AddInitializer
        (
         "multiplier",
         "protein model multiplier",
         io::Serialization::GetAgent( &m_Multiplier)
         );
      serializer.AddInitializer
        (
         "mutate model global xy translate",
         "translate the model along its xy axis",
         io::Serialization::GetAgent( &e_MutateModelGlobalXYTranslate)
         );
      serializer.AddInitializer
       (
        "mutate model global xyz translate",
        "translate the model in 3D space",
        io::Serialization::GetAgent( &e_MutateModelGlobalXYZTranslate)
        );
      serializer.AddInitializer
      (
        "mutate model global rotate mult",
        "rotate the entire multimer",
        io::Serialization::GetAgent( &e_MutateModelGlobalRotateMult)
      );
      serializer.AddInitializer
      (
        "mutate swap SSE multimer",
        "swap SSEs",
        io::Serialization::GetAgent( &e_MutateSwapSSEMultimer)
      );
      serializer.AddInitializer
      (
        "mutate model global z translate",
        "translate the model along its z axis",
        io::Serialization::GetAgent( &e_MutateModelGlobalZTranslate)
      );
      serializer.AddInitializer
      (
        "mutate model global z rotation",
        "rotate the model along its z axis",
        io::Serialization::GetAgent( &e_MutateModelGlobalZRotation)
      );

      return serializer;
    }

  ///////////
  // flags //
  ///////////

    //! @brief returns all flags that are specialized for this protocol
    //! @return all flags that are specialized for this protocol
    const util::ShPtrVector< command::FlagInterface> &ProtocolMultimer::GetAllFlags() const
    {
      // initialize static ShPtrVector of FlagInterfaces to form the comment line
      static util::ShPtrVector< command::FlagInterface> s_all_flags_vector;

      // if the flag vector is initialize for the first time
      if( s_all_flags_vector.IsEmpty())
      {
        // insert all the flags in the vector
        s_all_flags_vector.PushBack( GetFlagMultimer());
        s_all_flags_vector.PushBack( GetFlagNativeMultimer());
        s_all_flags_vector.PushBack( pdb::Factory::GetFlagBiomolecule());
      }

      // end
      return s_all_flags_vector;
    }

    //! @brief return command line flag for generating a multimeric protein model
    //! @return command line flag for generating a multimeric protein model
    util::ShPtr< command::FlagInterface> &ProtocolMultimer::GetFlagMultimer()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "symmetry", "\t\tsymmetry type to use to generate a multimeric protein model",
          command::Parameter
          (
            "symmetry_type",
            "\ttype of symmetry (i.e. C3 for 3-fold cyclic symmetry or D2 for 2-fold dihedral symmetry), "
              "currently only Cn and Dn symmetry is supported",
            "C2"
          )
        )
      );

      // end
      return s_flag;
    }

    //! @brief return command line flag for getting the native multimer model
    //! @return command line flag for getting the native multimer model
    util::ShPtr< command::FlagInterface> &ProtocolMultimer::GetFlagNativeMultimer()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "native_multimer", "\t\tnative multimer model for quality calculations",
          command::Parameter
          (
            "native_multimer_filename",
            "\tpath and filename for native multimer pdb",
            ""
          )
        )
      );

      // end
      return s_flag;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief modifies the start model
    //! @param START_MODEL Start model to be modified
    void ProtocolMultimer::ModifyStartModel( assemble::ProteinModel &START_MODEL) const
    {
      // if the biomolecule flag was set
      if( pdb::Factory::GetFlagBiomolecule()->GetFlag() && DefaultFlags::GetFlagStartModel()->GetFlag())
      {
        // recalculate the multiplier
        util::ShPtr< assemble::ProteinModelMultiplier> sp_multiplier
        (
          START_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Multiplier)
        );

        // make sure the multiplier is defined
        BCL_Assert( sp_multiplier.IsDefined(), "Multiplier is not defined");

        // make a new multiplier that uses the START_MODEL chains
        util::ShPtr< assemble::ProteinModelMultiplier> new_multiplier
        (
          new assemble::ProteinModelMultiplier
          (
            sp_multiplier->GetTransformationMatrices(),
            START_MODEL,
            true
          )
        );

        // get the protein model data
        util::ShPtr< assemble::ProteinModelData> pmd( START_MODEL.GetProteinModelData());

        // set the multiplier in the protein model data
        pmd->Replace( assemble::ProteinModelData::e_Multiplier, new_multiplier);

        // end
        return;
      }

      // get the protein model data
      util::ShPtr< assemble::ProteinModelData> model_data( START_MODEL.GetProteinModelData());

      // if the start model already has a multiplier
      if( model_data->GetData( assemble::ProteinModelData::e_Multiplier).IsDefined())
      {
        // don't modify further
        return;
      }

      // get the symmetry type
      Initialize();

      // move the model along the X-axis so that other subunits can be placed around the external Z-axis
      START_MODEL.Transform( math::Inverse( START_MODEL.GetOrientation()));
      const double offset( 2 * START_MODEL.GetSSEs().GetSize());
      START_MODEL.Translate( linal::Vector3D( offset, 0.0, 0.0));

      // if using dihedral symmetry
      if( m_IsDihedral)
      {
        // translate along the Z-axis as well
        START_MODEL.Translate( linal::Vector3D( 0.0, 0.0, offset));

        // update the protein model data with the multimer
        m_Multiplier =
          assemble::ProteinModelMultiplier
            (
              coord::GetAxes().e_Z,
              coord::GetAxes().e_X,
              m_CyclicSubunits,
              START_MODEL,
              true
          );
        util::ShPtr< assemble::ProteinModelMultiplier> s_multiplier( new assemble::ProteinModelMultiplier( m_Multiplier));
        model_data->Insert( assemble::ProteinModelData::e_Multiplier, s_multiplier);
      }
      // using cyclic symmetry
      else
      {
        // update the protein model data with the multimer
        m_Multiplier =
          assemble::ProteinModelMultiplier
            (
              coord::GetAxes().e_Z,
              m_CyclicSubunits,
              START_MODEL,
              true
          );
        util::ShPtr< assemble::ProteinModelMultiplier> s_multiplier( new assemble::ProteinModelMultiplier( m_Multiplier));
        model_data->Insert( assemble::ProteinModelData::e_Multiplier, s_multiplier);
      }

      // set the protein model data
      START_MODEL.SetProteinModelData( model_data);
    }

    //! @brief initialize the scores and add them to Scores enumerator
    void ProtocolMultimer::InitializeScores()
    {
    }

    //! @brief modify the score weight set
    //! @param SCORE_WEIGHT_SET Score weight set
    void ProtocolMultimer::ModifyScoreWeightSet( ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
    }

    //! @brief modify the terminate object
    //! @param CRITERION which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolMultimer::ModifyCriterion
    (
      opti::CriterionCombine< assemble::ProteinModel, double> &CRITERION,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief modify the printer object
    //! @param PRINTER which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolMultimer::ModifyPrinter
    (
      mc::PrinterCombined< assemble::ProteinModel, double> &PRINTER,
      const mc::Stage &STAGE
    ) const
    {
      // if a native multimer pdb file was given
      if( GetFlagNativeMultimer()->GetFlag())
      {
        // add the multimer printer
        PRINTER.Insert
        (
          util::ShPtr< mc::PrintInterface< assemble::ProteinModel, double> >
          (
            new assemble::PrinterProteinModelMultimer
            (
              GetSetup().GetPrefix(),
              GetNativeMultimer(),
              GetSetup().GetStorage(),
              GetSetup().GetSuperimposeMeasure()
            )
          )
        );
      }
    }

    //! @brief modify the pdb factory object
    //! @param FACTORY pdb factory to be modified
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolMultimer::ModifyFactory
    (
      util::ShPtr< pdb::Factory> &FACTORY,
      const mc::Stage &STAGE
    ) const
    {
      // if a native multimer pdb file was given
      if( GetFlagNativeMultimer()->GetFlag())
      {
        FACTORY->AppendPrinter
        (
          util::ShPtr< util::FunctionInterface< assemble::ProteinModel, util::ShPtrList< pdb::Line> > >
          (
            new pdb::PrinterQualityMultimer( GetSetup().GetQualityMeasures(), GetNativeMultimer())
          )
        );
      }
    }

    //! @brief initialize the mutates and add them to Mutates enumerator
    void ProtocolMultimer::InitializeMutates()
    {
      // global xy translation
      e_MutateModelGlobalXYTranslate = GetMutates().AddMutate
      (
        util::ShPtr< math::MutateInterface< assemble::ProteinModel> >
        (
          new MutateProteinModel
          (
            coord::MoveTranslateExternalAxis( 0.0, 5.0, coord::GetAxes().e_Z), "model_global_xy_translate"
          )
        )
      );

      // global xyz translation
      e_MutateModelGlobalXYZTranslate = GetMutates().AddMutate
      (
        util::ShPtr< math::MutateInterface< assemble::ProteinModel> >
        (
          new MutateProteinModel
          (
            coord::MoveTranslateRandom( 0.0, 5.0), "model_global_xyz_translate"
          )
        )
      );

      // global rotation
      e_MutateModelGlobalRotateMult = GetMutates().AddMutate
      (
        util::ShPtr< math::MutateInterface< assemble::ProteinModel> >
        (
          new MutateProteinModel( coord::MoveRotateRandom( math::g_Pi, true), "model_global_rotate_mult")
        )
      );

      // swap SSEs between subunits
      e_MutateSwapSSEMultimer = GetMutates().AddMutate
      (
        util::ShPtr< math::MutateInterface< assemble::ProteinModel> >
        (
          new MutateProteinModelSSESwapMultimer( DefaultFlags::GetFlagFitSwappedSSEs()->GetFlag(), "swap_sse_multimer")
        )
      );

      // if dihedral
      if( m_IsDihedral)
      {
        // z translation
        e_MutateModelGlobalZTranslate = GetMutates().AddMutate
        (
          util::ShPtr< math::MutateInterface< assemble::ProteinModel> >
          (
            new MutateProteinModel
            (
              coord::MoveTranslateRandom( 5.0, false),
              "model_global_z_translate"
            )
          )
        );

        // rotate around z_axis
        e_MutateModelGlobalZRotation = GetMutates().AddMutate
        (
          util::ShPtr< math::MutateInterface< assemble::ProteinModel> >
          (
            new MutateProteinModel
            (
              coord::MoveRotateRandom( linal::Vector3D(), linal::Vector3D( 0, 0, math::g_Pi / 3.0), false),
              "model_global_z_rotation"
            )
          )
        );
      }
    }

    //! @brief modify the mutate tree used
    //! @param MUTATE_TREE MutateTree to be modified
    void ProtocolMultimer::ModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
      MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_Model, 0.05);
      MUTATE_TREE.SetMutateProbability( e_MutateModelGlobalRotateMult, 1.0);

      // cyclic or dihedral symmetry
      if( util::IsDefined( m_CyclicSubunits))
      {
        MUTATE_TREE.SetMutateProbability( e_MutateSwapSSEMultimer, 1.0);
        MUTATE_TREE.SetMutateProbability( e_MutateModelGlobalXYTranslate, 1.0);
      }
      // symmetry read in from start model
      else
      {
        MUTATE_TREE.SetMutateProbability( e_MutateModelGlobalXYZTranslate, 1.0);
      }

      if( m_IsDihedral)
      {
        MUTATE_TREE.SetMutateProbability( e_MutateModelGlobalZTranslate, 1.0);
        MUTATE_TREE.SetMutateProbability( e_MutateModelGlobalZRotation, 1.0);
      }
    }

    //! @brief get the mutate tree associated with this protocol
    //! @return the mutate tree associated with this protocol
    util::ShPtr< MutateTree> ProtocolMultimer::GetMutateTree() const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( new MutateTree());
      ModifyMutateTree( *sp_mutate_tree);
      return sp_mutate_tree;
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void ProtocolMultimer::MergeAndModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }

  ////////////
  // readme //
  ////////////

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ProtocolMultimer::GetDescription() const
    {
      // initialize string to store the description
      static const std::string s_description
      (
        "protocol for folding multimeric proteins"
      );

      // end
      return s_description;
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ProtocolMultimer::GetReadMe() const
    {
      // initialize string to store the readme information
      static const std::string s_readme
      (
        "Adapts the BCL::Fold method to function with homo-multimeric proteins. The protomer is replicated around the "
        "axis of symmetry prior to scoring. The final PDB contains BIOMATRIX lines required to build the multimer."
        "Speficific flags:\n"
        "-symmetry #type of symmetry to use, i.e. C4 for 4-fold cyclic symmetry"
      );

      // end
      return s_readme;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProtocolMultimer::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_CyclicSubunits, ISTREAM);
      io::Serialize::Read( m_IsDihedral, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProtocolMultimer::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_CyclicSubunits, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_IsDihedral, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief initialize symmetry type
    void ProtocolMultimer::Initialize() const
    {
      // get the symmetry type from the command line
      const std::string symmetry_type( GetFlagMultimer()->GetFirstParameter()->GetValue());

      // assert that the string is long enough
      BCL_Assert
      (
        symmetry_type.length() > 1,
        "Incorrect symmetry type given (" + symmetry_type + "), "
          "must be 'C' or 'D' followed by the number of cyclic subunits"
      );

      // set m_IsDihedral based on the first character
      switch( symmetry_type[ 0])
      {
        case 'C':
          m_IsDihedral = false;
          break;
        case 'D':
          m_IsDihedral = true;
          break;
        default:
          BCL_Exit( "Incorrect symmetry type given, must begin with 'C' or 'D'", -1);
          break;
      }

      // get the rest of the string (i.e. the number of subunits)
      m_CyclicSubunits = util::ConvertStringToNumericalValue< size_t>
      (
        symmetry_type.substr( 1, symmetry_type.length() - 1)
      );
    }

    //! @brief copies the ss predictions from TEMPLATE_SEQUENCE onto TARGET_SEQUENCE
    //! @param TARGET_SEQUENCE sequence that will get ss predictions set
    //! @param TEMPLATE_SEQUENCE sequence containing existing ss predictions
    void ProtocolMultimer::SetSSPrediction
    (
      biol::AASequence &TARGET_SEQUENCE,
      const biol::AASequence &TEMPLATE_SEQUENCE
    )
    {
      // assert the sequences are same size
      BCL_Assert
      (
        TARGET_SEQUENCE.GetSize() == TEMPLATE_SEQUENCE.GetSize(),
        "Size of chain " + util::Format()( TARGET_SEQUENCE.GetChainID()) + " does not match size of chain " +
        util::Format()( TEMPLATE_SEQUENCE.GetChainID()) + ". The ss prediction files must be created for that chain."
      );

      // iterate over the sequences
      biol::AASequence::const_iterator template_itr( TEMPLATE_SEQUENCE.Begin());
      const biol::AASequence::const_iterator template_itr_end( TEMPLATE_SEQUENCE.End());
      for
      (
        biol::AASequence::iterator target_itr( TARGET_SEQUENCE.Begin()), target_itr_end( TARGET_SEQUENCE.End());
        target_itr != target_itr_end && template_itr != template_itr_end; ++target_itr, ++template_itr
      )
      {
        // set the predictions
        ( *target_itr)->SetSSPredictions( ( *template_itr)->GetSSPredictions());
      }
    }

  } // namespace fold
} // namespace bcl
