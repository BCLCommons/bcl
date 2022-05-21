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
#include "example_proteins.h"
// include the header of the class which this example is for
#include "density/bcl_density_protein_agreement_ccc.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_atom.h"
#include "density/bcl_density_map.h"
#include "density/bcl_density_simulators.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_density_protein_agreement_ccc.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDensityProteinAgreementCCC :
    public ExampleInterface
  {
  public:

    ExampleDensityProteinAgreementCCC *Clone() const
    { return new ExampleDensityProteinAgreementCCC( *this);}

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
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      BCL_MessageStd( "reading pdb: " + pdb_filename);

      //set parameters
      const double resolution( 6.6), voxelsize( 2.2);
      const util::ShPtr< density::SimulateInterface> simulator
      (
        density::GetSimulators().CreateSimulator
        (
          density::GetSimulators().e_Gaussian, linal::Vector3D( voxelsize), resolution
        )
      );

      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 0;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 0;
      assemble::ProteinModel model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size));
      ssetype_min_size[ biol::GetSSTypes().COIL] = 0;
      assemble::ProteinModel complete_model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AAComplete, ssetype_min_size));

      // density map of original model
      util::ShPtr< density::Map> sp_density_map( new density::Map());
      *sp_density_map = simulator->operator ()( complete_model.GetAtoms());

      density::ProteinAgreementCCC agreement_ccc( false, false);
      agreement_ccc.SetDensityMap( sp_density_map);
      agreement_ccc.SetSimulator( simulator);

      //calculate correlation between given densitymap and model densitymap and simulated densitymap
      const double expected_ccc_complete( -1);
      const double calculated_ccc_complete( agreement_ccc( complete_model));
      BCL_ExampleIndirectCheckWithinTolerance( calculated_ccc_complete, expected_ccc_complete, 0.001, "correlation original model AAComplete");

      //calculate correlation between given densitymap and model densitymap and simulated densitymap
      const double expected_ccc_backbone( -0.751326);
      const double calculated_ccc_backbone( agreement_ccc( model));
      BCL_ExampleIndirectCheckWithinTolerance( calculated_ccc_backbone, expected_ccc_backbone, 0.001, "correlation original model AABackBone");

      //idealize model
      model.SetToIdealConformation();

      //calculate correlation between given density map and model density map and simulated density map
      const double expected_ccc_ideal( -0.736202);
      const double calculated_ccc_ideal( agreement_ccc( model));
      BCL_ExampleIndirectCheckWithinTolerance( calculated_ccc_ideal, expected_ccc_ideal, 0.001, "correlation ideal model AABackBone");

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDensityProteinAgreementCCC

  const ExampleClass::EnumType ExampleDensityProteinAgreementCCC::s_Instance
  (
    GetExamples().AddEnum( ExampleDensityProteinAgreementCCC())
  );

} // namespace bcl
