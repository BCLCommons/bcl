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
#include "math/bcl_math_mutate_combine.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_quality.h"
#include "coord/bcl_coord_move_rotate_defined.h"
#include "coord/bcl_coord_move_translate_defined.h"
#include "fold/bcl_fold_mutate_protein_model_sse.h"
#include "math/bcl_math_mutate_move_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_mutate_combine.cpp
  //!
  //! @author woetzen
  //! @date
  //! @remarks status empty
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathMutateCombine :
    public ExampleInterface
  {
  public:

    ExampleMathMutateCombine *Clone() const
    { return new ExampleMathMutateCombine( *this);}

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
      // read start argument necessary for the constructor
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1C1D.pdb"));

      //build models from pdbsequences
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 5;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel protein_model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size));

      // list of mutates
      util::ShPtrList< math::MutateInterface< assemble::ProteinModel> > mutate_list;

      // create mutate object from defined translation and sse locator
      util::ShPtr< math::MutateInterface< assemble::ProteinModel> > mutate_trans_x
      (
        new fold::MutateProteinModelSSE
        (
          util::CloneToShPtr( assemble::LocatorSSE( 'A', 1, 7)),
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( coord::MoveTranslateDefined( linal::Vector3D( 0.0, 0.0, 3.0)), false))
        )
      );

      // create mutate object from defined translation and sse locator
      util::ShPtr< math::MutateInterface< assemble::ProteinModel> > mutate_trans_y
      (
        new fold::MutateProteinModelSSE
        (
          util::CloneToShPtr( assemble::LocatorSSE( 'A', 12, 18)),
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( coord::MoveTranslateDefined( linal::Vector3D( 0.0, 1.0, 0.0)), false))
        )
      );

      // create mutate object from a defined rotation and rotation axis, simulating a flip around x axis
      util::ShPtr< math::MutateInterface< assemble::ProteinModel> > mutate_flip_x
      (
        new fold::MutateProteinModelSSE
        (
          util::CloneToShPtr( assemble::LocatorSSE( 'A', 23, 31)),
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( coord::MoveRotateDefined( math::g_Pi, coord::GetAxes().e_X), false))
        )
      );

      // insert all mutates
      mutate_list.PushBack( mutate_trans_x);
      mutate_list.PushBack( mutate_trans_y);
      mutate_list.PushBack( mutate_flip_x);

    //////////////////
    // construction //
    //////////////////

      // default constructor
      math::MutateCombine< assemble::ProteinModel> empty_combine;

      // construct from list of mutates
      math::MutateCombine< assemble::ProteinModel> combine
      (
        mutate_list, false, math::MutateCombine< assemble::ProteinModel>::GetDefaultScheme()
      );

    ////////////////
    // operations //
    ////////////////

      // mutate
      assemble::ProteinModel mutated_model( *combine( protein_model).GetArgument());

      // write the current model to file
      const std::string out_filename( AddExampleOutputPathToFilename( empty_combine, "mutate_combined.pdb"));
      Proteins::WriteModelToPDB( mutated_model, out_filename);

      // print rmsd100
      // calculate RMSD of proteinmodel sses to coordinates given in chain
      const double rmsd( assemble::Quality::Calculate( quality::GetMeasures().e_RMSD, mutated_model, biol::GetAtomTypes().GetBackBoneAtomTypes()));
      BCL_MessageStd
      (
        "rmsd:    " + util::Format()( rmsd) +
        " rmsd100: " + util::Format()( assemble::Quality::RMSD100( rmsd, mutated_model.GetNumberAAs()))
      );

    //////////////////////
    // input and output //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathMutateCombine

  const ExampleClass::EnumType ExampleMathMutateCombine::s_Instance
  (
    GetExamples().AddEnum( ExampleMathMutateCombine())
  );

} // namespace bcl
