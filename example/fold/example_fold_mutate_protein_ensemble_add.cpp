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
#include "fold/bcl_fold_mutate_protein_ensemble_add.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_protein_ensemble_add.cpp
  //! TODO: add an detailed description for this example
  //!
  //!
  //! @author alexanns
  //! @date Feb 11, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinEnsembleAdd :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinEnsembleAdd *Clone() const
    {
      return new ExampleFoldMutateProteinEnsembleAdd( *this);
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
      // create two proteins
      const util::ShPtr< assemble::ProteinModel> ubi
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb")).Clone()
      );
      const util::ShPtr< assemble::ProteinModel> ie
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1IE9.pdb")).Clone()
      );

      // create protein ensemble and insert the two proteins
      assemble::ProteinEnsemble ensemble;
      ensemble.InsertElement( ubi);
      ensemble.InsertElement( ie);

      // create protein ensemble pool which can be used to add additional proteins to the above ensemble
      util::ShPtr< assemble::ProteinEnsemble> pool_ensemble( new assemble::ProteinEnsemble());
      pool_ensemble->InsertElement( ubi);
      pool_ensemble->InsertElement( ie);
      const util::ShPtr< assemble::ProteinModel> lzm
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "2LZM.pdb")).Clone()
      );
      pool_ensemble->InsertElement( lzm);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      BCL_MessageStd( "test default constructor");
      const fold::MutateProteinEnsembleAdd default_constr;
      BCL_ExampleCheck( default_constr( ensemble).GetArgument().IsDefined(), false);

      // constructor taking an ensemble that can be added from
      BCL_MessageStd( "test param constructor");
      const fold::MutateProteinEnsembleAdd param_constr( pool_ensemble);
      BCL_ExampleCheck( param_constr( ensemble).GetArgument()->GetSize(), 3);

      // clone constructor
      BCL_MessageStd( "test clone constructor");
      const util::ShPtr< fold::MutateProteinEnsembleAdd> clone_constr( param_constr.Clone());
      BCL_ExampleCheck( clone_constr->operator()( ensemble).GetArgument()->GetSize(), 3);

    /////////////////
    // data access //
    /////////////////

      BCL_MessageStd( "test GetStaticClassName GetClassIdentifier");
      // check GetStaticClassName
      const std::string correct_static_class_name( "bcl::fold::MutateProteinEnsembleAdd");
      BCL_ExampleCheck( GetStaticClassName< fold::MutateProteinEnsembleAdd>(), correct_static_class_name);

      // check GetClassIdentifier
      BCL_ExampleCheck( GetStaticClassName< fold::MutateProteinEnsembleAdd>(), clone_constr->GetClassIdentifier());

    ///////////////
    // operators //
    ///////////////

      // test operator to add protein to the ensemble
      {
        BCL_MessageStd( "test operator")
        const util::ShPtr< assemble::ProteinEnsemble> mutated_ensemble( param_constr( ensemble).GetArgument());
        BCL_ExampleCheck( mutated_ensemble->GetSize(), 3);
        // make sure 2lzm was added since it was the only one in the pool not in the ensemble already
        const assemble::ProteinEnsemble::const_iterator lzm_itr
        (
          std::find( mutated_ensemble->Begin(), mutated_ensemble->End(), lzm)
        );
        BCL_ExampleCheck( lzm_itr != mutated_ensemble->End(), true);
      }

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities");
      BCL_ExampleCheck( TestBCLObjectIOForSymmetry( param_constr, fold::MutateProteinEnsembleAdd()), true);

      std::stringstream input_output;
      // write the object
      input_output << param_constr;

      // read the object back in
      fold::MutateProteinEnsembleAdd read;
      input_output >> read;

      // compare the objects
      BCL_MessageStd( "compare written and read objects");
      {
        const util::ShPtr< assemble::ProteinEnsemble> mutated_ensemble( read( ensemble).GetArgument());
        BCL_ExampleCheck( mutated_ensemble->GetSize(), 3);
      }

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateProteinEnsembleAdd

  const ExampleClass::EnumType ExampleFoldMutateProteinEnsembleAdd::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinEnsembleAdd())
  );

} // namespace bcl
