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
#include "density/bcl_density_protein_agreement_likelihood.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_locator_sse_random.h"
#include "biol/bcl_biol_atom.h"
#include "coord/bcl_coord_move_transform_random.h"
#include "density/bcl_density_simulate_gaussian_sphere.h"
#include "fold/bcl_fold_mutate_protein_model_sse.h"
#include "math/bcl_math_quadratic_function.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_density_protein_agreement_likelihood.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDensityProteinAgreementLikelihood :
    public ExampleInterface
  {
  public:

    static const math::FunctionInterfaceSerializable< double, double> &GetMeanLowFit()
    {
      static const math::LinearFunction s_mean_low_fit( 0.00548655, 0.0613017);

      return s_mean_low_fit;
    }

    static const math::FunctionInterfaceSerializable< double, double> &GetSdLowFit()
    {
      static const math::QuadraticFunction s_sd_low_fit( 0.00027, 0.0041, 0.039);

      return s_sd_low_fit;
    }

    ExampleDensityProteinAgreementLikelihood *Clone() const
    { return new ExampleDensityProteinAgreementLikelihood( *this);}

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

    //////////////////
    // preparations //
    //////////////////

      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      BCL_MessageStd( "reading pdb: " + pdb_filename);

      // set parameters
      const double resolution( 6.6);
      const double gridspacing( resolution / 3.0);

      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 0;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 0;
      assemble::ProteinModel model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AAComplete, ssetype_min_size)
      );
      ssetype_min_size[ biol::GetSSTypes().COIL] = 0;
      assemble::ProteinModel density_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AAComplete, ssetype_min_size)
      );

      // density map of original model
      BCL_MessageStd( "Simulating the density map");
      util::ShPtr< density::SimulateInterface> simulator
      (
        new density::SimulateGaussianSphere( linal::Vector3D( gridspacing), resolution)
      );
      util::ShPtr< density::Map> sp_density_map( simulator->operator()( density_model.GetAtoms()).Clone());

      density::ProteinAgreementLikelihood score_density_likelihood
      (
        false,
        storage::Set< biol::AtomType>( biol::GetAtomTypes().CA),
        GetMeanLowFit(),
        GetSdLowFit()
      );
      score_density_likelihood.SetDensityMap( sp_density_map);
      score_density_likelihood.SetSimulator( simulator);

      const double score_model( score_density_likelihood( model));
      const double score_model_expected( -505.354);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( score_model, score_model_expected),
        "Score for the unchanged model " + util::Format()( score_model) + " doesn't match the expected one "
        + util::Format()( score_model_expected)
      );

      // generate mutated Protein models to get some more random values
      // create mutate object from max translation and rotation and sse locator
      fold::MutateProteinModelSSE mutate
      (
        util::CloneToShPtr( assemble::LocatorSSE( 'A', 1, 7)),
        util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( coord::MoveTransformRandom( 3.0, math::g_Pi), false))
      );
      // create mutate object "mutate_random" from max translation and rotation and SSELocatorRandom
      fold::MutateProteinModelSSE mutate_random
      (
        util::CloneToShPtr( assemble::LocatorSSERandom()),
        util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( coord::MoveTransformRandom( 3.0, math::g_Pi), false))
      );

      // mutate the protein model (random sse will be moved)
      assemble::ProteinModel mutated_model( *mutate( model).GetArgument());

      // create ProteinModel "mutated_model_random" mutate the protein model (random sse will be moved)
      assemble::ProteinModel mutated_model_random( *mutate_random( model).GetArgument());

      const double score_mutated_model( score_density_likelihood( mutated_model));
      const double score_mutated_model_expected( -385.479);
      const double score_mutated_model_random( score_density_likelihood( mutated_model_random));
      const double score_mutated_model_random_expected( -461.644);

      BCL_Example_Check
      (
        math::EqualWithinTolerance( score_mutated_model, score_mutated_model_expected, 0.02),
        "Score for the unchanged model " + util::Format()( score_mutated_model) + " doesn't match the expected one "
        + util::Format()( score_mutated_model_expected)
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( score_mutated_model_random, score_mutated_model_random_expected, 0.05),
        "Score for the unchanged model " + util::Format()( score_mutated_model_random)
        + " doesn't match the expected one " + util::Format()( score_mutated_model_random_expected)
      );
      assemble::ProteinModel mutated_model_second( *mutate( mutated_model_random).GetArgument());
      for( double res( 5.0); res < 21.0; res += 5.0)
      {
        util::ShPtr< density::SimulateInterface> simulator_res
        (
          new density::SimulateGaussianSphere( linal::Vector3D( res / 3.0), res)
        );
        util::ShPtr< density::Map> density_res_dependent( simulator_res->operator()( density_model.GetAtoms()).Clone());
        density::ProteinAgreementLikelihood score_density_likelihood_calc
        (
          false,
          storage::Set< biol::AtomType>( biol::GetAtomTypes().CA),
          GetMeanLowFit(),
          GetSdLowFit()
        );
        score_density_likelihood_calc.SetDensityMap( density_res_dependent);
        score_density_likelihood_calc.SetSimulator( simulator_res);
        BCL_MessageStd
        (
          "Value for likelihood score at resolution of " + util::Format()( res) + " A equals "
          + util::Format()( score_density_likelihood_calc( mutated_model_second))
        );
      }
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDensityProteinAgreementLikelihood

  const ExampleClass::EnumType ExampleDensityProteinAgreementLikelihood::s_Instance
  (
    GetExamples().AddEnum( ExampleDensityProteinAgreementLikelihood())
  );

} // namespace bcl
