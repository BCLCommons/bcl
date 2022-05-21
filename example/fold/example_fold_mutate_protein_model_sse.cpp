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
#include "fold/bcl_fold_mutate_protein_model_sse.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_locator_sse_random.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_quality.h"
#include "coord/bcl_coord_move_rotate_defined.h"
#include "coord/bcl_coord_move_transform_random.h"
#include "math/bcl_math_mutate_move_wrapper.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_protein_model_sse.cpp
  //!
  //! @author karakam
  //! @date Jan 25, 2010
  //! @remarks status empty
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinModelSSE :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinModelSSE *Clone() const
    { return new ExampleFoldMutateProteinModelSSE( *this);}

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
      // initialize pdb filename
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1C1D.pdb"));

      //build model
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 5;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create mutate object from max translation and rotation and sse locator
      fold::MutateProteinModelSSE mutate
      (
        util::CloneToShPtr( assemble::LocatorSSE( 'A', 1, 7)),
        util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( coord::MoveTransformRandom( 3.0, math::g_Pi), false))
      );

      // create mutate object from a defined rotation and rotation axis, simulating a flip around x axis
      fold::MutateProteinModelSSE mutate_flip
      (
        util::CloneToShPtr( assemble::LocatorSSE( 'A', 48, 70)),
        util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( coord::MoveRotateDefined( math::g_Pi, coord::GetAxes().e_X), false))
      );

      // create mutate object "mutate_random" from max translation and rotation and SSELocatorRandom
      fold::MutateProteinModelSSE mutate_random
      (
        util::CloneToShPtr( assemble::LocatorSSERandom()),
        util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( coord::MoveTransformRandom( 3.0, math::g_Pi), false))
      );

    ///////////////
    // operators //
    ///////////////

      // mutate the protein model (random sse will be moved
      assemble::ProteinModel mutated_model( *mutate( protein_model).GetArgument());

      // create ProteinModel "mutated_model_random" mutate the protein model (random sse will be moved)
      assemble::ProteinModel mutated_model_random( *mutate_random( protein_model).GetArgument());

      // create ProteinModel "mutated_model_random" mutate the protein model (random sse will be moved)
      assemble::ProteinModel mutated_model_flip( *mutate_flip( protein_model).GetArgument());

      // calculate the RMSDs for each model
      const storage::Set< biol::AtomType> backbone_atoms( biol::GetAtomTypes().GetBackBoneAtomTypes());
      const double rmsd_mutate
      (
        assemble::Quality::Calculate( quality::GetMeasures().e_RMSD, mutated_model, backbone_atoms)
      );
      const double rmsd_random
      (
        assemble::Quality::Calculate( quality::GetMeasures().e_RMSD, mutated_model_random, backbone_atoms)
      );
      const double rmsd_flip
      (
        assemble::Quality::Calculate( quality::GetMeasures().e_RMSD, mutated_model_flip, backbone_atoms)
      );

      BCL_MessageStd
      (
        "rmsd after random perturbation of sse 1-7 " + util::Format()( rmsd_mutate)
      );
      BCL_MessageStd
      (
        "rmsd after random sse random perturbation " + util::Format()( rmsd_random)
      );
      BCL_MessageStd
      (
        "rmsd after flipping sse 48-70 " + util::Format()( rmsd_flip)
      );

      // write mutated proteinmodel to an example pdb
      BCL_MessageStd( "write mutated.pdb");
      std::string out_filename( AddExampleOutputPathToFilename( mutate, "mutated.pdb"));
      Proteins::WriteModelToPDB( mutated_model, out_filename);

      // write mutated proteinmodel to an example pdb
      BCL_MessageStd( "mutated_random.pdb");
      out_filename = AddExampleOutputPathToFilename( mutate, "mutated_random.pdb");
      Proteins::WriteModelToPDB( mutated_model_random, out_filename);

      // write mutated proteinmodel to an example pdb
      BCL_MessageStd( "write mutated_flip.pdb");
      out_filename = AddExampleOutputPathToFilename( mutate, "mutated_flip.pdb");
      Proteins::WriteModelToPDB( mutated_model_flip, out_filename);

    //////////////////////
    // input and output //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateProteinModelSSE

  const ExampleClass::EnumType ExampleFoldMutateProteinModelSSE::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinModelSSE())
  );

} // namespace bcl
