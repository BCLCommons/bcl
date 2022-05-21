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
#include "restraint/bcl_restraint_mutate_transformation_matrix_3d_null.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_mutate_transformation_matrix_3d_null.cpp
  //!
  //! @author alexanns
  //! @date
  //! @remarks status empty
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintMutateTransformationMatrix3DNull :
    public ExampleInterface
  {
  public:

    ExampleRestraintMutateTransformationMatrix3DNull *Clone() const
    { return new ExampleRestraintMutateTransformationMatrix3DNull( *this);}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    virtual const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      // create string "pdb_filename" which has path for example pdb file
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi_ideal_model.pdb"));

      // create ProteinModel "protein_model" from "pdb"
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 0;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 0;
      assemble::ProteinModel protein_model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size));

      // get the helix secondary structure elements of "protein_model"
      util::SiPtrVector< const assemble::SSE> helix_sses
      (
        protein_model.GetChains()( 0)->GetSSEs( biol::GetSSTypes().HELIX)
      );
      BCL_MessageStd( "helix_sses size is " + util::Format()( helix_sses.GetSize()));

      // create ShPtr< assemble::SSE> sse_copy and initialize with the first helix of "helix_sses"
      util::ShPtr< assemble::SSE> sse_copy( ( *helix_sses.Begin())->Clone());

      BCL_MessageStd
      (
        "\n\n the transformation matrix of sse_copy before mutation is: \n"
        + util::Format()( sse_copy->GetOrientation())
      );

      // create restraint::MutateTransformationMatrix3DNull "mutate"
      restraint::MutateTransformationMatrix3DNull mutate;

      // mutate the transformation matrix of "sse_copy" with "mutate"
      sse_copy->Transform( *mutate( sse_copy->GetOrientation()).GetArgument());

      BCL_MessageStd
      (
        "\n\nthe transformation matrix of sse_copy after mutation is: \n"
        + util::Format()( sse_copy->GetOrientation())
      );

      // put the mutated "sse_copy" into "protein_model"
      protein_model.Replace( sse_copy);

      // write mutated proteinmodel to an example pdb
      BCL_MessageStd( "write mutated_transformation_matrix_3d_null.pdb");
      Proteins::WriteModelToPDB
      (
        protein_model, AddExampleOutputPathToFilename( mutate, "mutated_transformation_matrix_3d_null.pdb")
      );

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintMutateTransformationMatrix3DNull

  const ExampleClass::EnumType ExampleRestraintMutateTransformationMatrix3DNull::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintMutateTransformationMatrix3DNull())
  );

} // namespace bcl
