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
#include "assemble/bcl_assemble_sse_transformer.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_sse_transformer.cpp
  //!
  //! @author weinerbe
  //! @date Nov 16, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleSSETransformer :
    public ExampleInterface
  {
  public:

    ExampleAssembleSSETransformer *Clone() const
    {
      return new ExampleAssembleSSETransformer( *this);
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
      // create string "pdb_filename" which has path for example pdb file
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      // create protein model from pdb
      BCL_MessageStd( "building model");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 9;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

      // create sequence and transformation for new sses
      const util::ShPtr< math::TransformationMatrix3D> sp_transform
      (
        new math::TransformationMatrix3D
        (
          math::RotationMatrix3D( linal::Vector3D( 0.0, 0.0, 1.0), 2.0 * math::g_Pi / double( 5))
        )
      );
      util::ShPtr< biol::AASequence> sp_sequence( protein_model.GetSequences().FirstElement()->HardCopy());
      sp_sequence->SetChainID( 'B');
      sp_sequence->Transform( *sp_transform);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      assemble::SSETransformer def_construct;

      // test constructor from sequence and transformation
      assemble::SSETransformer seq_construct( sp_sequence, sp_transform);

      // test clone
      util::ShPtr< assemble::SSETransformer> clone_construct( seq_construct.Clone());

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck( clone_construct->GetClassIdentifier(), "bcl::assemble::SSETransformer");

    ///////////////
    // operators //
    ///////////////

      // check () operator
      util::ShPtr< assemble::SSE> transformed_sse( seq_construct( *( protein_model.GetSSEs().FirstElement())));
      BCL_ExampleIndirectCheck
      (
        transformed_sse->GetChainID() == 'B' &&
          transformed_sse->GetOrientation() != protein_model.GetSSEs().FirstElement()->GetOrientation(),
        true,
        "() operator"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // test write and read
      WriteBCLObject( seq_construct);
      assemble::SSETransformer read_construct;
      ReadBCLObject( read_construct);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleSSETransformer

  const ExampleClass::EnumType ExampleAssembleSSETransformer::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleSSETransformer())
  );
  
} // namespace bcl
