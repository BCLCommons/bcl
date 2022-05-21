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
#include "opti/bcl_opti_optimization_identity.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_ensemble.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically

namespace bcl
{

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opti_optimization_identity.cpp
  //! @brief tests the implementation of OptimizationIdentity
  //!
  //! @author fischea
  //! @date Nov 08, 2016
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleOptimizationIdentity :
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
    //! @return pointer to a new ExampleOptimizationIdentity
    ExampleOptimizationIdentity *Clone() const
    {
      return new ExampleOptimizationIdentity( *this);
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
      for( auto pdb_it( pdb_paths.Begin()), pdb_it_end( pdb_paths.End()); pdb_it != pdb_it_end; ++pdb_it)
      {
        // create a protein model from the current PDB file and add it three times to the ensemble
        // TODO create a more diverse set of structures instead of adding the same one multiple times
        const std::string &pdb_path( *pdb_it);
        assemble::ProteinModel model( pdb_factory.ProteinModelFromPDBFilename( pdb_path));
        for( size_t model_number( 0); model_number < 5; ++model_number)
        {
          ensemble.AddElement( model);
        }
      }

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      opti::OptimizationIdentity< assemble::ProteinModel> optimization_model;
      opti::OptimizationIdentity< assemble::Ensemble< assemble::ProteinModel> > optimization_ensemble;

    /////////////////
    // data access //
    /////////////////

      // test getter for class name identifier
      BCL_ExampleCheck
      (
        optimization_model.GetClassIdentifier(), ( GetStaticClassName< opti::OptimizationIdentity< assemble::ProteinModel> >())
      );

    ///////////////
    // operators //
    ///////////////

      // get the ensemble size before and after applying the algorithm
      const size_t ensemble_size_before( ensemble.GetSize());
      optimization_ensemble( ensemble);
      const size_t ensemble_size_after( ensemble.GetSize());

      // check if the resulting ensemble has the correct number of elements
      BCL_ExampleCheck( ensemble_size_before, ensemble_size_after);

      return 0;
    }

  }; // class ExampleOptimizationIdentity

  //! single instance of this class
  const ExampleClass::EnumType ExampleOptimizationIdentity::s_Instance
  (
    GetExamples().AddEnum( ExampleOptimizationIdentity())
  );

} // namespace bcl
