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

#ifndef BCL_FOLD_PROTOCOL_TEMPLATE_H_
#define BCL_FOLD_PROTOCOL_TEMPLATE_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_fold_mutates.h"
#include "bcl_fold_protocol_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProtocolTemplate
    //! @brief protocol for de novo folding of proteins using fold templates
    //! @details this protocol uses structural templates to efficiently sample complex topologies. therefore it uses a
    //! database of structures obtained from the protein data bank the target sequence can be fitted to and subsequently
    //! added to the protein model. This class is implemented as singleton, therefore only one instance of this class
    //! can exist at a given time.
    //!
    //! @example unnecessary
    //! @author weinerbe, fischea
    //! @date April 17, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProtocolTemplate :
      public ProtocolInterface
    {

    //////////
    // data //
    //////////

    private:

      //! selects a sheet from the protein model and fits it to a sheet template from the database
      Mutate e_MutateSheetFitToTemplate;

      //! selects strands from the SSE pool and fits them into sheet template from the database
      Mutate e_MutateAddSheetFromTemplate;

      //! selects a SSE from the protein model and bends it according to template from the database
//      Mutate e_MutateSSEBendTemplate;

    ///////////
    // flags //
    ///////////

    public:

      //! @brief returns all flags that are specialized for this protocol
      //! @return all flags that are specialized for this protocol
      const util::ShPtrVector< command::FlagInterface> &GetAllFlags() const;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    private:

      //! @brief default constructor
      ProtocolTemplate();

    public:

      //! @brief returns a pointer to a new ProtocolTemplate
      //! @return pointer to a new ProtocolTemplate
      ProtocolTemplate *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the only instance of this class
      //! @return the only instance of this class
      static ProtocolTemplate &GetInstance();

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

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
      void ModifyPrinter( mc::PrinterCombined< assemble::ProteinModel, double> &PRINTER, const mc::Stage &STAGE) const;

      //! @brief modify the pdb factory object
      //! @param FACTORY pdb factory to be modified
      //! @param STAGE Stage in which this terminate will be used
      void ModifyFactory( util::ShPtr< pdb::Factory> &FACTORY, const mc::Stage &STAGE) const;

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

      //! @brief returns a short description about this protocol
      //! @return string containing a short description about this protocol
      const std::string &GetDescription() const;

      //! @brief returns readme information about this protocol
      //! @return string containing readme information about this protocol
      const std::string &GetReadMe() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read object from input stream
      //! @param ISTREAM input stream to read object from
      //! @return input stream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write object into  output stream
      //! @param OSTREAM output stream to write object into
      //! @param INDENT number of indentations to separate members
      //! @return output stream object was written into
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class ProtocolTemplate

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_PROTOCOL_TEMPLATE_H_
