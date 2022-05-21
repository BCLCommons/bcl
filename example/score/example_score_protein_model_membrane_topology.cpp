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
#include "score/bcl_score_protein_model_membrane_topology.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_domain_sse_pool_overlapping.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_protein_model_membrane_topology.cpp
  //! @brief TODO: add an detailed description for this example
  //!
  //! @author alexanns
  //! @date Jun 15, 2012
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreProteinModelMembraneTopology :
    public ExampleInterface
  {
  public:

    ExampleScoreProteinModelMembraneTopology *Clone() const
    {
      return new ExampleScoreProteinModelMembraneTopology( *this);
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
      // create protein model
      BCL_MessageStd( "building model");
      std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "4A2N.pdb"));
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 5;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

      // create sse pool
      BCL_MessageStd( "read .pool");
      assemble::SSEPool sse_pool;
      io::IFStream read;
      std::string pool_filename( AddExampleInputPathToFilename( e_Biology, "4A2N.SSPredHighest_OCTOPUS.pool"));
      BCL_ExampleMustOpenInputFile( read, pool_filename);
      sse_pool.ReadSSEPool( read, protein_model, 3, 3);
      io::File::CloseClearFStream( read);

      // make a domain locator
      util::ShPtr< assemble::LocatorDomainSSEPoolOverlapping> domain_locator
      (
        new assemble::LocatorDomainSSEPoolOverlapping( sse_pool)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      score::ProteinModelMembraneTopology def_constr;
      BCL_ExampleCheck( def_constr.GetScheme(), GetStaticClassName< score::ProteinModelMembraneTopology>());

      // constructor taking parameters
      const std::string scheme( "topology_score");
      score::ProteinModelMembraneTopology param_constr( domain_locator, sse_pool, scheme);
      BCL_ExampleCheck( param_constr.GetScheme(), scheme);

      // clone constructor
      util::ShPtr< score::ProteinModelMembraneTopology> clone_constr( param_constr.Clone());
      BCL_ExampleCheck( clone_constr->GetScheme(), scheme);

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck
      (
        GetStaticClassName< score::ProteinModelMembraneTopology>(), clone_constr->GetClassIdentifier()
      );

      // GetScheme
      BCL_ExampleCheck( param_constr.GetScheme(), scheme);

    ////////////////
    // operations //
    ////////////////

      // operator()
      {
        const double score( param_constr( protein_model));
        const double expected_score( 0);
        BCL_ExampleCheck( score, expected_score);
      }

      // operator()
      {
        // create protein model
        BCL_MessageStd( "building model");
        std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "4A2N_74_100_flipped.pdb"));
        storage::Map< biol::SSType, size_t> ssetype_min_size;
        ssetype_min_size[ biol::GetSSTypes().HELIX] = 5;
        ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
        assemble::ProteinModel protein_model
        (
          Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
        );

        const double score( param_constr( protein_model));
        const double expected_score( 8);
        BCL_ExampleCheck( score, expected_score);
      }

      // check score with multi chain protein / tm regions
      {
        // create protein model
        BCL_MessageStd( "building model");
        std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1L0V.pdb"));
        storage::Map< biol::SSType, size_t> ssetype_min_size;
        ssetype_min_size[ biol::GetSSTypes().HELIX] = 5;
        ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
        assemble::ProteinModel protein_model
        (
          Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
        );

        // create sse pool
        BCL_MessageStd( "read .pool");
        assemble::SSEPool sse_pool;
        io::IFStream read;
        std::string pool_filename( AddExampleInputPathToFilename( e_Biology, "1L0V.SSPredMC_OCTOPUS.pool"));
        BCL_ExampleMustOpenInputFile( read, pool_filename);
        sse_pool.ReadSSEPool( read, protein_model, 3, 3);
        io::File::CloseClearFStream( read);

        // make a domain locator
        util::ShPtr< assemble::LocatorDomainSSEPoolOverlapping> domain_locator
        (
          new assemble::LocatorDomainSSEPoolOverlapping( sse_pool)
        );

        score::ProteinModelMembraneTopology param_constr( domain_locator, sse_pool);
        const double score( param_constr( protein_model));
        const double expected_score( 0);
        BCL_ExampleCheck( score, expected_score);
      }

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

      // test read and write
      WriteBCLObject( *clone_constr);
      score::ProteinModelMembraneTopology read_score;
      ReadBCLObject( read_score);
      {
        const double score( read_score( protein_model));
        const double expected_score( 0);
        BCL_ExampleCheck( score, expected_score);
      }
      BCL_ExampleCheck( read_score.GetScheme(), scheme);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreProteinModelMembraneTopology

  const ExampleClass::EnumType ExampleScoreProteinModelMembraneTopology::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreProteinModelMembraneTopology())
  );
  
} // namespace bcl
