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

#ifndef BCL_APP_TEMPLATE_MODEL_H_
#define BCL_APP_TEMPLATE_MODEL_H_

// include the interface for all apps
#include "app/bcl_app_apps.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_mutate_protein_model_thread_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class TemplateModel
    //! @brief application for making a template model
    //! @details Uses an alignment to assign coordinates from a known structure onto a sequence of unknown structure.
    //! A sequence alignment between a sequence with known coordinates and a sequence with unknown coordinates is
    //! used to assign the coordinates from the given protein model to the second sequence in the alignment. The
    //! protein model must have the same sequence as the first sequence in the alignment. The provided start
    //! model is used as the template whose coordinates are used. The first sequence in the alignment must
    //! correspond to the sequence of the template start model. The output file generated is
    //! "<prefix>threaded_<template_pdb_name_without_path>"
    //!
    //! @see @link example_app_TemplateModel.cpp @endlink
    //! @author alexanns
    //! @date Sep 20, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API TemplateModel :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      //! if set, structurally problematic residues will be printed to logger
      util::ShPtr< command::FlagInterface> m_FlagPrintProblemAAs;

    public:

      // instantiate enumerator for GenerateDataset class
      static const ApplicationType TemplateModel_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      TemplateModel();

      //! @brief Clone function
      //! @return pointer to new TemplateModel
      TemplateModel *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

    ///////////////
    // operators //
    ///////////////

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

    private:

      //! @brief gives the mutate for threading sequence
      //! @return mutate for threading sequence onto template structure
      fold::MutateProteinModelThreadSequence GetMutate() const;

    }; // class TemplateModel

  } // namespace app
} // namespace bcl

#endif // BCL_APP_TEMPLATE_MODEL_H_
