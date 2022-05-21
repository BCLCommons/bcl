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
#include "sspred/bcl_sspred_pdb.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_aa.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_membrane.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_sspred_pdb.cpp
  //!
  //! @author weinerbe
  //! @date Jan 14, 2011
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleSspredPDB :
    public ExampleInterface
  {
  public:

    ExampleSspredPDB *Clone() const
    {
      return new ExampleSspredPDB( *this);
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
      // read in protein model
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "2K73A.pdb"))
      );

      // create a membrane object
      const util::ShPtr< biol::Membrane> sp_membrane( new biol::Membrane());

      // add the membrane to protein data
      util::ShPtr< assemble::ProteinModelData> sp_data( protein_model.GetProteinModelData());
      sp_data->Insert( assemble::ProteinModelData::e_Membrane, sp_membrane);
      protein_model.SetProteinModelData( sp_data);

      // create a few locators
      const assemble::LocatorAA locator_24( 'A', 24);
      const assemble::LocatorAA locator_109( 'A', 109);
      const assemble::LocatorAA locator_136( 'A', 136);

    ////////////////
    // operations //
    ////////////////

      // check SetEnvironmentTypes
      sspred::PDB::SetEnvironmentTypes( protein_model);
      util::SiPtr< const sspred::MethodInterface> sp_pred_24
      (
        locator_24.Locate( protein_model)->GetSSPrediction( sspred::GetMethods().e_PDB)
      );
      BCL_ExampleCheck( sp_pred_24->GetOneStateSSPrediction(), biol::GetSSTypes().HELIX);
      BCL_ExampleCheck( sp_pred_24->GetOneStateTMPrediction(), biol::GetEnvironmentTypes().e_MembraneCore);
      util::SiPtr< const sspred::MethodInterface> sp_pred_109
      (
        locator_109.Locate( protein_model)->GetSSPrediction( sspred::GetMethods().e_PDB)
      );
      BCL_ExampleCheck( sp_pred_109->GetOneStateSSPrediction(), biol::GetSSTypes().COIL);
      BCL_ExampleCheck( sp_pred_109->GetOneStateTMPrediction(), biol::GetEnvironmentTypes().e_GapCoreTransition);
      util::SiPtr< const sspred::MethodInterface> sp_pred_136
      (
        locator_136.Locate( protein_model)->GetSSPrediction( sspred::GetMethods().e_PDB)
      );
      BCL_ExampleCheck( sp_pred_136->GetOneStateSSPrediction(), biol::GetSSTypes().STRAND);
      BCL_ExampleCheck( sp_pred_136->GetOneStateTMPrediction(), biol::GetEnvironmentTypes().e_Transition);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSspredPDB

  const ExampleClass::EnumType ExampleSspredPDB::s_Instance
  (
    GetExamples().AddEnum( ExampleSspredPDB())
  );

} // namespace bcl
