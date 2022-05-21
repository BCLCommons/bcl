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
#include "fold/bcl_fold_mutate_protein_model_domain.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_collector_sheet.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "coord/bcl_coord_move_rotate_defined.h"
#include "coord/bcl_coord_move_transform_random.h"
#include "coord/bcl_coord_move_translate_defined.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_protein_model_domain.cpp
  //!
  //! @author karakam
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinModelDomain :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinModelDomain *Clone() const
    {
      return new ExampleFoldMutateProteinModelDomain( *this);
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
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      BCL_MessageStd( "test default constructor");
      fold::MutateProteinModelDomain mutate_def;

      // initialize pdb filename
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "2qv3_sheet_ideal.pdb"));

      // get the protein model
      assemble::ProteinModel model( Proteins::GetModel( pdb_filename));

      // create sheet collector
      util::ShPtr< find::CollectorInterface< util::ShPtrVector< assemble::Domain>, assemble::ProteinModel> >
        collector_sheet( new assemble::CollectorSheet());

      coord::MoveTransformRandom transform( linal::Vector3D( 5, 0, 0), linal::Vector3D( 0, 0, math::g_Pi / 2), true);
      coord::MoveRotateDefined rotate( math::g_Pi, coord::GetAxes().e_X, true);
      coord::MoveTranslateDefined translate_x( linal::Vector3D( 15, 0, 0), true);
      coord::MoveTranslateDefined translate_y( linal::Vector3D( 0, 15, 0), true);
      coord::MoveTranslateDefined translate_z( linal::Vector3D( 0, 0, 15), true);

      fold::MutateProteinModelDomain
        mutate_transform( *collector_sheet, transform);

      fold::MutateProteinModelDomain
        mutate_translate_x( *collector_sheet, translate_x);

      fold::MutateProteinModelDomain
        mutate_translate_y( *collector_sheet, translate_y);

      fold::MutateProteinModelDomain
        mutate_translate_z( *collector_sheet, translate_z);

      fold::MutateProteinModelDomain
        mutate_rotate( *collector_sheet, rotate);

      math::MutateResult< assemble::ProteinModel> result_a( mutate_transform( model));
      Proteins::WriteModelToPDB
      (
        *result_a.GetArgument(), AddExampleOutputPathToFilename( mutate_def, "mutate_domain_wrap_transform.pdb")
      );

      math::MutateResult< assemble::ProteinModel> result_bx( mutate_translate_x( model));
      Proteins::WriteModelToPDB
      (
        *result_bx.GetArgument(), AddExampleOutputPathToFilename( mutate_def, "mutate_domain_wrap_translate_x.pdb")
      );

      math::MutateResult< assemble::ProteinModel> result_by( mutate_translate_y( model));
      Proteins::WriteModelToPDB
      (
        *result_by.GetArgument(), AddExampleOutputPathToFilename( mutate_def, "mutate_domain_wrap_translate_y.pdb")
      );

      math::MutateResult< assemble::ProteinModel> result_bz( mutate_translate_z( model));
      Proteins::WriteModelToPDB
      (
        *result_bz.GetArgument(), AddExampleOutputPathToFilename( mutate_def, "mutate_domain_wrap_translate_z.pdb")
      );

      math::MutateResult< assemble::ProteinModel> result_c( mutate_rotate( model));
      Proteins::WriteModelToPDB
      (
        *result_c.GetArgument(), AddExampleOutputPathToFilename( mutate_def, "mutate_domain_wrap_rotate.pdb")
      );

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateProteinModelDomain

  const ExampleClass::EnumType ExampleFoldMutateProteinModelDomain::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinModelDomain())
  );
  
} // namespace bcl
