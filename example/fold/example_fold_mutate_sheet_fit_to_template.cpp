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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "fold/bcl_fold_mutate_sheet_fit_to_template.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_collector_sheet.h"
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_sheet_fit_to_template.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateSheetFitToTemplate :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateSheetFitToTemplate *Clone() const
    {
      return new ExampleFoldMutateSheetFitToTemplate( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

      // construct the map of min sse lengths
      storage::Map< biol::SSType, size_t> min_sse_lengths;
      min_sse_lengths[ biol::GetSSTypes().HELIX] = 7;
      min_sse_lengths[ biol::GetSSTypes().STRAND] = 4;
      min_sse_lengths[ biol::GetSSTypes().COIL] = 999;

      // read the idealized 1UBI model
      assemble::ProteinModel model
      (
        Proteins::GetModel
        (
          AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"),
          biol::GetAAClasses().e_AABackBone,
          min_sse_lengths
        )
      );

      // make a copy and idealize the model
      util::ShPtr< assemble::ProteinModel> ideal_model_ptr( model.HardCopy());
      assemble::ProteinModel &ideal_model( *ideal_model_ptr);
      ideal_model.SetToIdealConformation();

      // get the sheet from the ideal model
      util::ShPtrVector< assemble::Domain> sheet_vector( assemble::CollectorSheet().Collect( ideal_model));
      const assemble::Domain &ideal_sheet( *sheet_vector.FirstElement());

      // print out the sheet
      BCL_MessageStd( "sheet selected \n:" + ideal_sheet.GetTopology()->GetOrderedIdentification());

      // number of strands in the sheet
      const size_t nr_strands( 4);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      fold::MutateSheetFitToTemplate mutate;

      // test clone constructor
      util::ShPtr< fold::MutateSheetFitToTemplate> sp_mutate( mutate.Clone());

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck( GetStaticClassName< fold::MutateSheetFitToTemplate>(), mutate.GetClassIdentifier());

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // test operator()
      const math::MutateResult< assemble::Domain> result( mutate( ideal_sheet));

      // make sure the returned result is defined
      BCL_ExampleAssert( result.GetArgument().IsDefined(), true);

      // create a reference on the returned domain
      const assemble::Domain &new_sheet( *result.GetArgument());

      // make sure it has the same number of strands
      BCL_ExampleCheck( new_sheet.GetNumberSSEs(), nr_strands);

      // construct a chain from the returned domain
      util::ShPtr< assemble::Chain> sp_new_chain
      (
        new assemble::Chain( ideal_model.GetChains().FirstElement()->GetSequence(), new_sheet)
      );

      // construct a model from the new chain
      assemble::ProteinModel new_model( sp_new_chain);

      // write out the fitted sheet
      const std::string pdb_filename( AddExampleOutputPathToFilename( mutate, "mutate_sheet_fit_to_template.pdb"));
      Proteins::WriteModelToPDB( new_model, pdb_filename);

      BCL_ExampleIndirectCheck
      (
        new_model.GetSSEs( biol::GetSSTypes().STRAND).GetSize(),
        4,
        "model should have 4 strands but has " +
          util::Format()( new_model.GetSSEs( biol::GetSSTypes().STRAND).GetSize())
      );

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateSheetFitToTemplate

  const ExampleClass::EnumType ExampleFoldMutateSheetFitToTemplate::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateSheetFitToTemplate())
  );

} // namespace bcl
