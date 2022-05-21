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
#include "assemble/bcl_assemble_biomolecule.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "math/bcl_math_mutate_result.h"
#include "quality/bcl_quality_superimpose_measures.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_biomolecule.cpp
  //! @brief TODO: add an detailed description for this example
  //!
  //! @author woetzen
  //! @date Jul 27, 2012
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleBiomolecule :
    public ExampleInterface
  {
  public:

    ExampleAssembleBiomolecule *Clone() const
    {
      return new ExampleAssembleBiomolecule( *this);
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
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1j4n_multimer.pdb"));
      // create ProteinModel "protein_model" from "pdb"
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 0;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 0;
      ssetype_min_size[ biol::GetSSTypes().COIL] = 0;

      const assemble::ProteinModel native_protein
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AAComplete, ssetype_min_size)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      assemble::Biomolecule biomolecule( biol::GetAtomTypes().GetBackBoneAtomTypes(), *quality::GetSuperimposeMeasures().e_RMSD, 1.0);

    /////////////////
    // data access //
    /////////////////

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      math::MutateResult< assemble::ProteinModel> result( biomolecule( native_protein));
      BCL_ExampleCheck( result.GetArgument()->GetNumberOfChains(), 1);
      Proteins::WriteModelToPDB( *result.GetArgument(), AddExampleOutputPathToFilename( biomolecule, "biomolecule1j4n.pdb"));

    //////////////////////
    // input and output //
    //////////////////////

      // test write and read
      WriteBCLObject( biomolecule);
      assemble::Biomolecule
      read_biomolecule;
      ReadBCLObject( biomolecule);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleBiomolecule

  const ExampleClass::EnumType ExampleAssembleBiomolecule::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleBiomolecule())
  );

} // namespace bcl
