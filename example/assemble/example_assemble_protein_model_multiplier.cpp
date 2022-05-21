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
#include "assemble/bcl_assemble_protein_model_multiplier.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_protein_model_multiplier.cpp
  //!
  //! @author weinerbe
  //! @date Nov 12, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleProteinModelMultiplier :
    public ExampleInterface
  {
  public:

    ExampleAssembleProteinModelMultiplier *Clone() const
    {
      return new ExampleAssembleProteinModelMultiplier( *this);
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
      const std::string small_pdb_filename( AddExampleInputPathToFilename( e_Biology, "2HAC.pdb"));

      // create protein model from pdb
      BCL_MessageStd( "building model");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 9;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );
      assemble::ProteinModel small_model
      (
        Proteins::GetModel( small_pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );
      protein_model.Transform( math::Inverse( protein_model.GetOrientation()));
      protein_model.Translate( linal::Vector3D( 16.0, 0.0, 0.0));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      assemble::ProteinModelMultiplier def_construct;

      // test cyclic constructor
      assemble::ProteinModelMultiplier cyclic_construct
      (
        linal::Vector3D( 0.0, 0.0, 1.0),
        5,
        protein_model
      );

      // test cached construct
      assemble::ProteinModelMultiplier cached_construct
      (
        linal::Vector3D( 0.0, 0.0, 1.0),
        5,
        protein_model,
        true
      );

      // test clone
      util::ShPtr< assemble::ProteinModelMultiplier> clone_construct( cyclic_construct.Clone());
      BCL_ExampleIndirectCheck
      (
        clone_construct->GetChainMultipliers().GetSize(),
        cyclic_construct.GetChainMultipliers().GetSize(),
        "Clone"
      );

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck( clone_construct->GetClassIdentifier(), "bcl::assemble::ProteinModelMultiplier");

    ////////////////
    // operations //
    ////////////////

      // mapping of chain ids for from original to multiplied model
      const storage::Table< char> mapping( cyclic_construct.ChainIDMapping());
      mapping.WriteFormatted( util::GetLogger());
      BCL_ExampleCheck( cyclic_construct.GetChainMultipliers().GetSize(), mapping.GetNumberRows());

    ///////////////
    // operators //
    ///////////////

      // test () operator
      assemble::ProteinModel multimer( cyclic_construct( protein_model));
      BCL_ExampleIndirectCheck( multimer.GetNumberOfChains(), 5, "() operator");
      assemble::ProteinModel multimer_cached( cached_construct( protein_model));
      BCL_ExampleIndirectCheck( multimer_cached.GetNumberOfChains(), 5, "() operator cached");

    //////////////////////
    // input and output //
    //////////////////////

      // test read and write using a small protein multimer
      assemble::ProteinModelMultiplier small_construct
      (
        linal::Vector3D( 0.0, 0.0, 1.0),
        2,
        small_model
      );

      WriteBCLObject( small_construct);
      assemble::ProteinModelMultiplier read_construct;
      ReadBCLObject( read_construct);
      BCL_ExampleIndirectCheck
      (
        small_construct.GetChainMultipliers().GetSize(),
        read_construct.GetChainMultipliers().GetSize(),
        "Read and Write"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleProteinModelMultiplier

  const ExampleClass::EnumType ExampleAssembleProteinModelMultiplier::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleProteinModelMultiplier())
  );

} // namespace bcl
