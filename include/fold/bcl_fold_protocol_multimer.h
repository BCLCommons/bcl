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

#ifndef BCL_FOLD_PROTOCOL_MULTIMER_H_
#define BCL_FOLD_PROTOCOL_MULTIMER_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_fold_mutates.h"
#include "bcl_fold_protocol_interface.h"
#include "assemble/bcl_assemble_protein_model_multiplier.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProtocolMultimer
    //! @brief protocol for folding multimeric proteins
    //! @details this protocol wraps the scores used in default protocol so that coordinates for all subunits are
    //! calculated on the go and are considered when scoring. It also introduces global moves as well as special
    //! modifications to the starting model
    //!
    //! @remarks example unnecessary
    //! @author karakam
    //! @date Nov 17, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProtocolMultimer :
      public ProtocolInterface
    {

    private:

    //////////
    // data //
    //////////

      //! number of cyclic subunits
      mutable size_t m_CyclicSubunits;

      //! true if dihedral symmetry, false if cyclic symmetry
      mutable bool m_IsDihedral;

      //! protein model multiplier
      mutable assemble::ProteinModelMultiplier m_Multiplier;

    public:

      Mutate e_MutateModelGlobalXYTranslate;
      Mutate e_MutateModelGlobalXYZTranslate;
      Mutate e_MutateModelGlobalRotateMult;
      Mutate e_MutateSwapSSEMultimer;
      Mutate e_MutateModelGlobalZTranslate;
      Mutate e_MutateModelGlobalZRotation;

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ProtocolMultimer();

    public:

      //! @brief Clone function
      //! @return pointer to new ProtocolMultimer
      ProtocolMultimer *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief access to only instance
      //! @return reference to only instance
      static ProtocolMultimer &GetInstance();

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the native multimer
      //! @return the native multimer
      const util::ShPtr< assemble::ProteinModel> &GetNativeMultimer() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ///////////
    // flags //
    ///////////

      //! @brief returns all flags that are specialized for this protocol
      //! @return all flags that are specialized for this protocol
      const util::ShPtrVector< command::FlagInterface> &GetAllFlags() const;

      //! @brief return command line flag for generating a multimeric protein model
      //! @return command line flag for generating a multimeric protein model
      static util::ShPtr< command::FlagInterface> &GetFlagMultimer();

      //! @brief return command line flag for getting the native multimer model
      //! @return command line flag for getting the native multimer model
      static util::ShPtr< command::FlagInterface> &GetFlagNativeMultimer();

    ////////////////
    // operations //
    ////////////////

      //! @brief modifies the start model
      //! @param START_MODEL Start model to be modified
      void ModifyStartModel( assemble::ProteinModel &START_MODEL) const;

      //! @brief initialize the scores and add them to Scores enumerator
      void InitializeScores();

      //! @brief modify the score weight set
      //! @param SCORE_WEIGHT_SET Score weight set
      void ModifyScoreWeightSet( ScoreWeightSet &SCORE_WEIGHT_SET) const;

      //! @brief modify the terminate object
      //! @param CRITERION which will be modified by protocols
      //! @param STAGE Stage in which this terminate will be used
      void ModifyCriterion
      (
        opti::CriterionCombine< assemble::ProteinModel, double> &CRITERION,
        const mc::Stage &STAGE
      ) const;

      //! @brief modify the printer object
      //! @param PRINTER which will be modified by protocols
      //! @param STAGE Stage in which this terminate will be used
      void ModifyPrinter
      (
        mc::PrinterCombined< assemble::ProteinModel, double> &PRINTER,
        const mc::Stage &STAGE
      ) const;

      //! @brief modify the pdb factory object
      //! @param FACTORY pdb factory to be modified
      //! @param STAGE Stage in which this terminate will be used
      void ModifyFactory
      (
        util::ShPtr< pdb::Factory> &FACTORY,
        const mc::Stage &STAGE
      ) const;

      //! @brief initialize the mutates and add them to Mutates enumerator
      void InitializeMutates();

      //! @brief modify the mutate tree used
      //! @param MUTATE_TREE MutateTree to be modified
      void ModifyMutateTree( MutateTree &MUTATE_TREE) const;

      //! @brief get the mutate tree associated with this protocol
      //! @return the mutate tree associated with this protocol
      util::ShPtr< MutateTree> GetMutateTree() const;

      //! @brief merges this protocol's mutate tree into given mutate tree
      //! @param MUTATE_TREE tree into which to merge this protocol's tree
      void MergeAndModifyMutateTree( MutateTree &MUTATE_TREE) const;

    ////////////
    // readme //
    ////////////

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetDescription() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief initialize symmetry type
      void Initialize() const;

    private:

      //! @brief copies the ss predictions from TEMPLATE_SEQUENCE onto TARGET_SEQUENCE
      //! @param TARGET_SEQUENCE sequence that will get ss predictions set
      //! @param TEMPLATE_SEQUENCE sequence containing existing ss predictions
      static void SetSSPrediction( biol::AASequence &TARGET_SEQUENCE, const biol::AASequence &TEMPLATE_SEQUENCE);

    }; // class ProtocolMultimer

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_PROTOCOL_MULTIMER_H_ 
