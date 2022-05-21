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
#include "math/bcl_math_mutate_perturbation.h"

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
  //! @example example_math_mutate_perturbation.cpp
  //!
  //! @author woetzen
  //! @date
  //! @remarks status empty
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathMutatePerturbation :
    public ExampleInterface
  {

  public:

  //////////
  // data //
  //////////

    //! single instance of this class
    static const ExampleClass::EnumType s_Instance;

    ExampleMathMutatePerturbation *Clone() const
    {
      return new ExampleMathMutatePerturbation( *this);
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
      // read start argument necessary for the constructor
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1C1D.pdb"));
      //build models from pdbsequences
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 5;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel protein_model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size));

    //////////////////
    // construction //
    //////////////////

      // default constructor
      math::MutatePerturbation< assemble::ProteinModel> empty_perturbation;

      // construct from start argument
      math::MutatePerturbation< assemble::ProteinModel> perturbation( protein_model);

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

    ////////////////
    // operations //
    ////////////////

      // insert all mutates with the considered number of repetitions
      empty_perturbation.InsertPerturbation( mutate_trans_x, 2);
      empty_perturbation.InsertPerturbation( mutate_trans_y, 2);
      empty_perturbation.InsertPerturbation( mutate_flip_x, 2);

      // set the starting value
      empty_perturbation.SetStartArgument( protein_model);

      size_t iteration_nr( 0);
      while( !empty_perturbation.IsFinished())
      {
        // perturb
        empty_perturbation.Perturb();

        // get current
        const assemble::ProteinModel &current( empty_perturbation.GetCurrent());

        // write all models if in debug mode
        if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
        {
          // write the current model to file
          const std::string out_filename
          (
            AddExampleOutputPathToFilename( empty_perturbation, "perturbed_" + util::Format()( iteration_nr) + ".pdb")
          );
          Proteins::WriteModelToPDB( current, out_filename);
        }

        // print rmsd100
        // calculate RMSD of protein model sses to coordinates given in chain
        const double rmsd( assemble::Quality::Calculate( quality::GetMeasures().e_RMSD, current, biol::GetAtomTypes().GetBackBoneAtomTypes()));
        BCL_MessageStd
        (
          "rmsd:    " + util::Format()( rmsd) +
          " rmsd100: " + util::Format()( assemble::Quality::RMSD100( rmsd, current.GetNumberAAs()))
        );
        ++iteration_nr;
      }

    //////////////////////
    // input and output //
    //////////////////////

      return 0;
    }

  }; //end ExampleMathMutatePerturbation

  const ExampleClass::EnumType ExampleMathMutatePerturbation::s_Instance
  (
    GetExamples().AddEnum( ExampleMathMutatePerturbation())
  );

} // namespace bcl
