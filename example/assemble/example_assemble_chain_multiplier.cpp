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
#include "assemble/bcl_assemble_chain_multiplier.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_transformer.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_chain_multiplier.cpp
  //!
  //! @author weinerbe
  //! @date Nov 16, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleChainMultiplier :
    public ExampleInterface
  {
  public:

    ExampleAssembleChainMultiplier *Clone() const
    {
      return new ExampleAssembleChainMultiplier( *this);
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
      util::ShPtr< assemble::SSETransformer> sp_sse_transformer
      (
        new assemble::SSETransformer( sp_sequence, sp_transform)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      assemble::ChainMultiplier def_construct;

      // test constructor from sse transformer
      assemble::ChainMultiplier transformer_construct( sp_sse_transformer, 'A', sp_transform, sp_sequence);

      // test clone
      util::ShPtr< assemble::ChainMultiplier> clone_construct( transformer_construct.Clone());
      BCL_ExampleIndirectCheck
      (
        clone_construct->GetInitialChainID(),
        transformer_construct.GetInitialChainID(),
        "Clone"
      );

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck( clone_construct->GetClassIdentifier(), "bcl::assemble::ChainMultiplier");

      // check GetInitialChainID
      BCL_ExampleCheck( transformer_construct.GetInitialChainID(), 'A');

      // check GetNewChainID
      BCL_ExampleCheck( transformer_construct.GetNewChainID(), 'B');

      // check GetTransformationMatrix
      BCL_ExampleCheck( transformer_construct.GetTransformationMatrix(), sp_transform);

    ///////////////
    // operators //
    ///////////////

      // check () operator
      util::ShPtr< assemble::Chain> multiplied_chain
      (
        transformer_construct( *( protein_model.GetChains().FirstElement()))
      );
      BCL_ExampleIndirectCheck
      (
        multiplied_chain->GetChainID() == 'B' &&
        multiplied_chain->GetOrientation() != protein_model.GetChains().FirstElement()->GetOrientation(),
        true,
        "() operator"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // test write and read
      WriteBCLObject( transformer_construct);
      assemble::ChainMultiplier read_construct;
      ReadBCLObject( read_construct);
      BCL_ExampleIndirectCheck
      (
        transformer_construct.GetInitialChainID() == read_construct.GetInitialChainID() &&
          transformer_construct.GetNewChainID() == read_construct.GetNewChainID(),
        true,
        "Read and Write"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleChainMultiplier

  const ExampleClass::EnumType ExampleAssembleChainMultiplier::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleChainMultiplier())
  );
  
} // namespace bcl
