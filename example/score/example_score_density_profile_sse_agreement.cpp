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
#include "score/bcl_score_density_profile_sse_agreement.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_atom.h"
#include "density/bcl_density_map.h"
#include "density/bcl_density_simulators.h"
#include "io/bcl_io_file.h"
#include "restraint/bcl_restraint_body.h"
#include "restraint/bcl_restraint_contains_body_origin.h"
#include "restraint/bcl_restraint_handler_body.h"
#include "score/bcl_score_body_assignment.h"
#include "score/bcl_score_restraint_body_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_density_profile_sse_agreement.cpp
  //!
  //! @author linders
  //! @date
  //! @remarks status empty
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreDensityProfileSSEAgreement :
    public ExampleInterface
  {
  public:

    ExampleScoreDensityProfileSSEAgreement *Clone() const
    {
      return new ExampleScoreDensityProfileSSEAgreement( *this);
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
      // create string "restraint_filename" which has path for example restraint file
      const std::string restraint_filename( AddExampleInputPathToFilename( e_Biology, "1ubi_ideal_model.pdb"));

      // create stream to restraint file
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, restraint_filename);
      BCL_MessageStd( "read restraint file: " + restraint_filename);

      // create body restraints
      // create HandlerAtomDistanceAssigned "handler" as method for determining if a restraint body is occupied
      restraint::HandlerBody handler
      (
        util::ShPtr< util::BinaryFunctionInterface< assemble::SSEGeometryInterface, assemble::SSE, bool> >
        (
          new restraint::ContainsBodyOrigin()
        )
      );

      // create restraints
      const util::ShPtr< restraint::Body> restraints( handler.CreateRestraintsBody( read).FirstElement());
      io::File::CloseClearFStream( read);

      // create ideal protein model
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 9;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      const assemble::ProteinModel ideal_model( Proteins::GetModel( restraint_filename, biol::GetAAClasses().e_AAComplete, ssetype_min_size));

      // get pdb
      // create string "pdb_filename" which has path for example pdb file
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      // create ProteinModel "protein_model" from "pdb"
      BCL_MessageStd( "building models from pdb chains and sse information");
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 0;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 0;
      ssetype_min_size[ biol::GetSSTypes().COIL] = 0;
      assemble::ProteinModel protein_model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AAComplete, ssetype_min_size));

      // make density map and sse's from pdb
      const double resolution( 6.6), voxelsize( 2.2);
      const util::ShPtr< density::SimulateInterface> simulator
      (
        density::GetSimulators().CreateSimulator
        (
          density::GetSimulators().e_Gaussian, linal::Vector3D( voxelsize), resolution
        )
      );
      // density map of original model
      util::ShPtr< density::Map> sp_density_map( new density::Map());
      *sp_density_map = simulator->operator ()( protein_model.GetAtoms());

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct from simulator, density and restraint for each of the three measures
      util::ShPtr< score::DensityProfileSSEAgreement> angle_1D_score_profile_agreement
      (
        new score::DensityProfileSSEAgreement( simulator, *sp_density_map, *restraints, score::DensityProfileSSEAgreement::e_Angle)
      );

      util::ShPtr< score::DensityProfileSSEAgreement> height_1D_score_profile_agreement
      (
        new score::DensityProfileSSEAgreement( simulator, *sp_density_map, *restraints, score::DensityProfileSSEAgreement::e_Height)
      );

      util::ShPtr< score::DensityProfileSSEAgreement> radius_1D_score_profile_agreement
      (
        new score::DensityProfileSSEAgreement( simulator, *sp_density_map, *restraints, score::DensityProfileSSEAgreement::e_Radius)
      );

    /////////////////
    // data access //
    /////////////////

      // print static classname
      BCL_MessageStd
      (
        "static classname is " + GetStaticClassName< score::DensityProfileSSEAgreement>()
      );

    ///////////////
    // operators //
    ///////////////

      // create an assignment score for angle_1D
      const score::BodyAssignment angle_1D_assignment_score( angle_1D_score_profile_agreement);

      // create an assignment score for height_1D
      const score::BodyAssignment height_1D_assignment_score( height_1D_score_profile_agreement);

      // create an assignment score for radius_1D
      const score::BodyAssignment radius_1D_assignment_score( radius_1D_score_profile_agreement);

      // initialize ShPtr with a new ShPtrVector of restraint interfaces which are body restraints
      util::ShPtr
      <
        util::ShPtrVector< restraint::Body>
      > restraints_sp
      (
        new util::ShPtrVector< restraint::Body>( size_t( 1), restraints)
      );

      // create angle_1D protein model score
      const score::RestraintBodyProteinModel angle_1D_score_model( restraints_sp, angle_1D_assignment_score);

      // create height_1D protein model score
      const score::RestraintBodyProteinModel height_1D_score_model( restraints_sp, height_1D_assignment_score);

      // create radius_1D protein model score
      const score::RestraintBodyProteinModel radius_1D_score_model( restraints_sp, radius_1D_assignment_score);

      const double angle_1D_ideal_model_score( angle_1D_score_model( ideal_model));

      const double height_1D_ideal_model_score( height_1D_score_model( ideal_model));

      const double radius_1D_ideal_model_score( radius_1D_score_model( ideal_model));

      BCL_MessageStd( "angle_1D profile score of ideal model: " + util::Format()( angle_1D_ideal_model_score));

      BCL_MessageStd( "height_1D profile score of ideal model: " + util::Format()( height_1D_ideal_model_score));

      BCL_MessageStd( "radius_1D profile score of ideal model: " + util::Format()( radius_1D_ideal_model_score));

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

  }; //end ExampleScoreDensityProfileSSEAgreement

  const ExampleClass::EnumType ExampleScoreDensityProfileSSEAgreement::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreDensityProfileSSEAgreement())
  );

} // namespace bcl
