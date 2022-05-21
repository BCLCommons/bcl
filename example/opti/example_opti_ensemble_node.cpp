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
#include "opti/bcl_opti_ensemble_node.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "mc/bcl_mc_optimization_mcm.h"
#include "pdb/bcl_pdb_factory.h"
#include "score/bcl_score_protein_model_completeness.h"

// external includes - sorted alphabetically

namespace bcl
{

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opti_ensemble_node.cpp
  //! @brief tests the implementation of EnsembleNode
  //!
  //! @author fischea
  //! @date Oct 30, 2016
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleOptiEnsembleNode :
    public ExampleInterface
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    static const ExampleClass::EnumType s_Instance;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleEnsembleNode
    ExampleOptiEnsembleNode *Clone() const
    {
      return new ExampleOptiEnsembleNode( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief run routine
    //! this is performing the execution of the example
    int Run() const
    {

    //////////////////////
    // data preparation //
    //////////////////////

      // read in the list of protein models to optimize
      const std::string pdb_list_filename( AddExampleInputPathToFilename( e_Mc, "loop_hash_test_set.ls"));
      io::IFStream pdb_list;
      BCL_ExampleMustOpenInputFile( pdb_list, pdb_list_filename);
      storage::Vector< std::string> pdb_paths( util::StringLineListFromIStream( pdb_list));
      io::File::CloseClearFStream( pdb_list);

      // create the initial ensemble
      pdb::Factory pdb_factory;
      assemble::Ensemble< assemble::ProteinModel> ensemble;
      const size_t ensemble_size( 3);
      for( auto pdb_it( pdb_paths.Begin()), pdb_it_end( pdb_paths.End()); pdb_it != pdb_it_end; ++pdb_it)
      {
        // create a protein model from the current PDB file and add it three times to the ensemble
        // TODO create a more diverse set of structures instead of adding the same one multiple times
        const std::string &pdb_path( *pdb_it);
        assemble::ProteinModel model( pdb_factory.ProteinModelFromPDBFilename( pdb_path));
        ensemble.AddElement( model);
      }

      // read in the loop hash optimizer
      mc::OptimizationMCM optimizer_hash;
      const std::string optimizer_hash_file_name( AddExampleInputPathToFilename( e_Mc, "optimizer_loop_hash"));
      io::IFStream optimizer_hash_file;
      BCL_ExampleMustOpenInputFile( optimizer_hash_file, optimizer_hash_file_name);
      optimizer_hash_file >> optimizer_hash;
      io::File::CloseClearFStream( optimizer_hash_file);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      opti::EnsembleNode< assemble::ProteinModel> default_node;

      // construct from optimization implementation
      opti::EnsembleNode< assemble::ProteinModel> node( optimizer_hash, ensemble_size);

    /////////////////
    // data access //
    /////////////////

      // test getter for class name identifier
      BCL_ExampleCheck
      (
        default_node.GetClassIdentifier(), ( GetStaticClassName< opti::EnsembleNode< assemble::ProteinModel> >())
      );

    ///////////////
    // operators //
    ///////////////

      // construct loop regions for the protein models in the ensemble
      node( ensemble);

      // check if the resulting ensemble has the correct number of elements
      BCL_ExampleCheck( ensemble.GetSize(), ensemble_size);

      // check if loop regions have been constructed
      score::ProteinModelCompleteness compl_score( true);
      for( auto model_it( ensemble.Begin()), model_it_end( ensemble.End()); model_it != model_it_end; ++model_it)
      {
        // check the current model for completeness
        assemble::Ensemble< assemble::ProteinModel>::Element &element( *model_it);
        const assemble::ProteinModel model( element.GetElement());
        const double completeness( compl_score( model));
        BCL_ExampleCheckWithinTolerance( completeness, -1.0, 0.01);
      }

      return 0;
    }

  }; // class ExampleEnsembleNode

  //! single instance of this class
  const ExampleClass::EnumType ExampleOptiEnsembleNode::s_Instance
  (
    GetExamples().AddEnum( ExampleOptiEnsembleNode())
  );

} // namespace bcl
