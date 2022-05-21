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
#include "fold/bcl_fold_mutate_multimer.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_protein_model_multiplier.h"
#include "quality/bcl_quality_rmsd.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_multimer.cpp
  //!
  //! @author weinerbe
  //! @date Dec 3, 2012
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateMultimer :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateMultimer *Clone() const
    {
      return new ExampleFoldMutateMultimer( *this);
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
      // get a pdb
      assemble::ProteinModel protein_model( Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb")));

      // create the multimer
      util::ShPtr< assemble::ProteinModelMultiplier> sp_multiplier
      (
        new assemble::ProteinModelMultiplier
        (
          coord::GetAxes().e_Z,
          4,
          protein_model
        )
      );

      // insert the multiplier
      util::ShPtr< assemble::ProteinModelData> sp_data( protein_model.GetProteinModelData());
      sp_data->Insert( assemble::ProteinModelData::e_Multiplier, sp_multiplier);
      protein_model.SetProteinModelData( sp_data);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test constructor
      const std::string scheme( "test_mutate");
      fold::MutateMultimer def_construct;
      const fold::MutateMultimer mutate( false, scheme);

      // initialize from the string
      const std::string parameters( "MutateMultimer(is dihedral=False)");
      util::Implementation< math::MutateInterface< assemble::ProteinModel> > mutate_ser;
      BCL_ExampleAssert( mutate_ser.TryRead( parameters, util::GetLogger()), true);

    /////////////////
    // data access //
    /////////////////

      // test GetScheme
      BCL_ExampleCheck( mutate.GetAlias(), "MutateMultimer");

    ////////////////
    // operations //
    ////////////////

      // test () operator
      const util::ShPtr< assemble::ProteinModel> sp_mutated_model( ( *mutate_ser)( protein_model).GetArgument());
      const util::SiPtrVector< const linal::Vector3D> orig_coords( protein_model.GetAtomCoordinates());
      const util::SiPtrVector< const linal::Vector3D> new_coords( sp_mutated_model->GetAtomCoordinates());
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinAbsoluteTolerance
        (
          quality::RMSD( false).CalculateMeasure( orig_coords, new_coords),
          0.0,
          0.00001
        ),
        false,
        "() operator"
      );

    //////////////////////
    // input and output //
    //////////////////////

      WriteBCLObject( mutate);
      ReadBCLObject( def_construct);
      BCL_ExampleCheck( def_construct.GetLabel(), mutate.GetLabel());

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateMultimer

  const ExampleClass::EnumType ExampleFoldMutateMultimer::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateMultimer())
  );

} // namespace bcl
