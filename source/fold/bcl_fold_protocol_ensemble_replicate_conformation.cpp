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
#include "fold/bcl_fold_protocol_ensemble_replicate_conformation.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_collector_protein_model_conformations.h"
#include "assemble/bcl_assemble_printer_protein_model_ensemble.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "fold/bcl_fold_mutate_protein_model_replicate_conformation.h"
#include "fold/bcl_fold_protocol_ensemble.h"
#include "fold/bcl_fold_setup.h"
#include "math/bcl_math_mutate_result.h"
#include "mc/bcl_mc_printer_combined.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProtocolEnsembleReplicateConformation::s_Instance
    (
      GetObjectInstances().AddInstance( new ProtocolEnsembleReplicateConformation())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProtocolEnsembleReplicateConformation::ProtocolEnsembleReplicateConformation()
    {
    }

    //! @brief Clone function
    //! @return pointer to new ProtocolEnsembleReplicateConformation
    ProtocolEnsembleReplicateConformation *ProtocolEnsembleReplicateConformation::Clone() const
    {
      return new ProtocolEnsembleReplicateConformation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief access to only instance
    //! @return reference to only instance
    ProtocolEnsembleReplicateConformation &ProtocolEnsembleReplicateConformation::GetInstance()
    {
      static ProtocolEnsembleReplicateConformation s_protocol_instance;
      return s_protocol_instance;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProtocolEnsembleReplicateConformation::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ProtocolEnsembleReplicateConformation::GetAlias() const
    {
      static const std::string s_name( "ProtocolEnsembleReplicateConformation");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProtocolEnsembleReplicateConformation::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "protocol to be applied for replicating the conformation of a model that is going to be folded.");

      return serializer;
    }

  ///////////
  // flags //
  ///////////

    //! @brief returns all flags that are specialized for this protocol
    //! @return all flags that are specialized for this protocol
    const util::ShPtrVector< command::FlagInterface> &ProtocolEnsembleReplicateConformation::GetAllFlags() const
    {
      // initialize static ShPtrVector of FlagInterfaces to form the comment line
      static util::ShPtrVector< command::FlagInterface> s_all_flags_vector;

      // if the flag vector is initialize for the first time
      if( s_all_flags_vector.IsEmpty())
      {
      }

      // end
      return s_all_flags_vector;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief modifies the start model
    //! @param START_MODEL Start model to be modified
    void ProtocolEnsembleReplicateConformation::ModifyStartModel( assemble::ProteinModel &START_MODEL) const
    {
      // add one since the current conformation of focus is not included
      const int num_conformations( START_MODEL.GetConformationalEnsemble().GetSize() + 1);
      const int wanted_size( ProtocolEnsemble::GetFlagEnsembleSize()->GetFirstParameter()->GetNumericalValue< int>());

      // subtract one since if the division gives 1, no replications should take place
      const int num_replications( ( wanted_size - num_conformations) / num_conformations);

      BCL_MessageStd( "number of replications is " + util::Format()( num_replications));

      // create object for replicating the models in the ensemble
      MutateProteinModelReplicateConformation mutate
      (
        util::ShPtr< assemble::CollectorProteinModelConformations>
        (
          new assemble::CollectorProteinModelConformations( true)
        ),
        num_replications
      );

      // copy the models
      math::MutateResult< assemble::ProteinModel> mutate_result( mutate( START_MODEL));

      BCL_MessageStd
      (
        "ensemble size is " +
        util::Format()( mutate_result.GetArgument()->GetConformationalEnsemble().GetSize())
      );

      // set the start model to the result of the mutate - ensemble should have desired size and number of copies
      START_MODEL = ( *mutate_result.GetArgument());
    }

    //! @brief initialize the scores and add them to Scores enumerator
    void ProtocolEnsembleReplicateConformation::InitializeScores()
    {
    }

    //! @brief modify the score weight set
    //! @param SCORE_WEIGHT_SET Score weight set
    void ProtocolEnsembleReplicateConformation::ModifyScoreWeightSet( ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
    }

    //! @brief modify the terminate object
    //! @param CRITERION which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolEnsembleReplicateConformation::ModifyCriterion
    (
      opti::CriterionCombine< assemble::ProteinModel, double> &CRITERION,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief modify the printer object
    //! @param PRINTER which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolEnsembleReplicateConformation::ModifyPrinter
    (
      mc::PrinterCombined< assemble::ProteinModel, double> &PRINTER,
      const mc::Stage &STAGE
    ) const
    {
      util::ShPtr< mc::PrintInterface< assemble::ProteinModel, double> > ensemble_printer
      (
        new assemble::PrinterProteinModelEnsemble
        (
          GetSetup().GetPrefix(),
          GetSetup().GetStorage(),
          GetSetup().GetSuperimposeMeasure()
        )
      );
      PRINTER.Insert( ensemble_printer);
    }

    //! @brief modify the pdb factory object
    //! @param FACTORY pdb factory to be modified
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolEnsembleReplicateConformation::ModifyFactory
    (
      util::ShPtr< pdb::Factory> &FACTORY,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief initialize the mutates and add them to Mutates enumerator
    void ProtocolEnsembleReplicateConformation::InitializeMutates()
    {
    }

    //! @brief modify the mutate tree used
    //! @param MUTATE_TREE MutateTree to be modified
    void ProtocolEnsembleReplicateConformation::ModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
    }

    //! @brief get the mutate tree associated with this protocol
    //! @return the mutate tree associated with this protocol
    util::ShPtr< MutateTree> ProtocolEnsembleReplicateConformation::GetMutateTree() const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( new MutateTree());
      ModifyMutateTree( *sp_mutate_tree);
      return sp_mutate_tree;
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void ProtocolEnsembleReplicateConformation::MergeAndModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }

  ////////////
  // readme //
  ////////////

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ProtocolEnsembleReplicateConformation::GetDescription() const
    {
      // initialize string to store the description
      static const std::string s_description
      (
        "protocol to be applied for replicating the conformation of a model that is going to be folded."
      );

      // end
      return s_description;
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ProtocolEnsembleReplicateConformation::GetReadMe() const
    {
      // TODO: add readme for this protocol
      // initialize string to store the readme information
      static const std::string s_readme
      (
        "readme for replicating the conformation of a model that is being folded when the ensemble protocol is used."
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
    std::istream &ProtocolEnsembleReplicateConformation::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProtocolEnsembleReplicateConformation::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace fold
} // namespace bcl
