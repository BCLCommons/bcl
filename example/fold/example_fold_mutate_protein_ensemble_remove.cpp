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
#include "fold/bcl_fold_mutate_protein_ensemble_remove.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_protein_ensemble_remove.cpp
  //!
  //! @author alexanns
  //! @date Feb 12, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinEnsembleRemove :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinEnsembleRemove *Clone() const
    {
      return new ExampleFoldMutateProteinEnsembleRemove( *this);
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
      // create protein ensemble
      assemble::ProteinEnsemble ensemble;

      ensemble.InsertElement
      (
        util::ShPtr< assemble::ProteinModel>
        (
          Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb")).Clone()
        )
      );
      ensemble.InsertElement
      (
        util::ShPtr< assemble::ProteinModel>
        (
          Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1IE9.pdb")).Clone()
        )
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      const fold::MutateProteinEnsembleRemove default_constr;
      BCL_ExampleCheck( default_constr( ensemble).GetArgument()->GetSize(), 1);

      // clone constructor
      const util::ShPtr< fold::MutateProteinEnsembleRemove> clone_constr( default_constr.Clone());
      BCL_ExampleCheck( clone_constr->operator()( ensemble).GetArgument()->GetSize(), 1);

    /////////////////
    // data access //
    /////////////////

      // check GetStaticClassName
      const std::string correct_static_class_name( "bcl::fold::MutateProteinEnsembleRemove");
      BCL_ExampleCheck( GetStaticClassName< fold::MutateProteinEnsembleRemove>(), correct_static_class_name);

      // check GetClassIdentifier
      BCL_ExampleCheck( GetStaticClassName< fold::MutateProteinEnsembleRemove>(), clone_constr->GetClassIdentifier());

    ///////////////
    // operators //
    ///////////////

      // test operator to remove an element from the ensemble
      BCL_ExampleCheck( clone_constr->operator()( ensemble).GetArgument()->GetSize(), 1);

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities");
      // write the object
      WriteBCLObject( default_constr);

      // read the object back in
      fold::MutateProteinEnsembleRemove read;
      ReadBCLObject( read);

      // compare the objects
      BCL_MessageStd( "compare written and read objects");
      BCL_ExampleCheck( read( ensemble).GetArgument()->GetSize(), 1);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateProteinEnsembleRemove

  const ExampleClass::EnumType ExampleFoldMutateProteinEnsembleRemove::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinEnsembleRemove())
  );

} // namespace bcl
