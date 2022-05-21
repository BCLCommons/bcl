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
#include "fold/bcl_fold_mutate_protein_model_sse_split.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "math/bcl_math_mutate_result.h"
#include "sspred/bcl_sspred_method_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_protein_model_sse_split.cpp
  //!
  //! @author karakam, woetzen
  //! @date Aug 3, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinModelSSESplit :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinModelSSESplit *Clone() const
    {
      return new ExampleFoldMutateProteinModelSSESplit( *this);
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
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      //build models from pdb sequences
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 3;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 2;
      ssetype_min_size[ biol::GetSSTypes().COIL] = 999;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

      // method
      sspred::Method method( sspred::GetMethods().e_JUFO);
      sspred::MethodHandler::ReadPredictionsForProteinModel
      (
        storage::Set< sspred::Method>( method), protein_model, "1ubi", AddExampleInputPathToFilename( e_Biology, "")
      );

      // locator for helix and strand
      util::ShPtr< assemble::LocatorSSE> sp_locator_helix( new assemble::LocatorSSE( 'A', 23, 34));
      util::ShPtr< assemble::LocatorSSE> sp_locator_strand( new assemble::LocatorSSE( 'A', 64, 72));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      const math::Range< size_t> coil_length_range( 1, 1);

      fold::MutateProteinModelSSESplit splitter_helix( method, sp_locator_helix, ssetype_min_size, coil_length_range);
      fold::MutateProteinModelSSESplit splitter_strand( method, sp_locator_strand, ssetype_min_size, coil_length_range);

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

      // calling operator to mutate helix
      const math::MutateResult< assemble::ProteinModel> result_helix( splitter_helix( protein_model));

      // check if result is defined
      BCL_ExampleAssert( result_helix.GetArgument().IsDefined(), true);

      // write the mutated result
      Proteins::WriteModelToPDB
      (
        *result_helix.GetArgument(), AddExampleOutputPathToFilename( splitter_helix, "1ubi_split_helix.pdb")
      );

      // calling operator to mutate strand
      const math::MutateResult< assemble::ProteinModel> result_strand( splitter_strand( protein_model));

      // check if result is defined
      BCL_ExampleAssert( result_strand.GetArgument().IsDefined(), true);

      // write the mutated result
      Proteins::WriteModelToPDB
      (
        *result_strand.GetArgument(), AddExampleOutputPathToFilename( splitter_strand, "1ubi_split_strand.pdb")
      );

    //////////////////////
    // input and output //
    //////////////////////

      // read and write the object
      WriteBCLObject( splitter_helix);
      fold::MutateProteinModelSSESplit splitter_read
      (
        sspred::GetMethods().e_Undefined,
        util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> >(),
        storage::Map< biol::SSType, size_t>(),
        math::Range< size_t>()
      );
      ReadBCLObject( splitter_read);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateProteinModelSSESplit

  const ExampleClass::EnumType ExampleFoldMutateProteinModelSSESplit::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinModelSSESplit())
  );

} // namespace bcl
