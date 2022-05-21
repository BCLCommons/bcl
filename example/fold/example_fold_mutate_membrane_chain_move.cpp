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
#include "fold/bcl_fold_mutate_membrane_chain_move.h"

// includes from bcl - sorted alphabetically
#include "example_check_macros.h"
#include "example_interface.h"
#include "example_proteins.h"
#include "biol/bcl_biol_aa_classes.h"
#include "math/bcl_math_mutate_result.h"
#include "util/bcl_util_enumerate.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example bcl_fold_mutate_membrane_chain_move.cpp
  //!
  //! @author lib14
  //! @date Sept 26, 2017
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateMembraneChainMove :
      public ExampleInterface
  {

  public:

    ExampleFoldMutateMembraneChainMove *Clone() const
    {
      return new ExampleFoldMutateMembraneChainMove( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get the name of this class
    //! @return the name of this class
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      // read in a testing membrane complex
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "4o6y.pdb"));

      // create an object of ProteinModel
      const assemble::ProteinModel native_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AAComplete)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct from default constructor
      fold::MutateMembraneChainMove def_mem_chain_move;

      // test SetReceptorLigandChainIDs()
      BCL_ExampleCheck( def_mem_chain_move.GetLigandChainIDs(), "");
      def_mem_chain_move.SetReceptorLigandChainIDs( native_model);
      BCL_ExampleCheck( def_mem_chain_move.GetLigandChainIDs(), "A");

      // construct from given arguments
      const std::string receptor_chain_ids( "A");
      const std::string ligand_chain_ids( "B");
      double max_translation_x( 1.1);
      double max_translation_y( 1.2);
      double max_translation_z( 1.3);
      double max_rotation_phi( 0.1);
      double max_rotation_theta( 0.2);
      double max_rotation_psi( 0.3);
      fold::MutateMembraneChainMove mem_chain_move
      (
        receptor_chain_ids,
        ligand_chain_ids,
        max_translation_x,
        max_translation_y,
        max_translation_z,
        max_rotation_phi,
        max_rotation_theta,
        max_rotation_psi
      );

      // test getters
      BCL_ExampleCheck( mem_chain_move.GetReceptorChainIDs(), "A");
      BCL_ExampleCheck( mem_chain_move.GetLigandChainIDs(), "B");
      BCL_ExampleCheck( mem_chain_move.GetMaxTranslationX(), 1.1);
      BCL_ExampleCheck( mem_chain_move.GetMaxTranslationY(), 1.2);
      BCL_ExampleCheck( mem_chain_move.GetMaxTranslationZ(), 1.3);
      BCL_ExampleCheck( mem_chain_move.GetMaxPhi(), 0.1);
      BCL_ExampleCheck( mem_chain_move.GetMaxTheta(), 0.2);
      BCL_ExampleCheck( mem_chain_move.GetMaxPsi(), 0.3);

      // test mutate operator
      // assemble::ProteinModel mutated_model( *mem_chain_move( native_model).GetArgument());
      // Proteins::WriteModelToPDB( mutated_model, AddExampleOutputPathToFilename( mem_chain_move, "4o6y_mutated.pdb"));

    //////////////////////
    // input and output //
    //////////////////////

      return 0;
    }

    // single instance of this class
    static const ExampleClass::EnumType s_Instance;
  }; // end of class ExampleFoldMutateMembraneChainMove

  // single instance of this class
  const ExampleClass::EnumType ExampleFoldMutateMembraneChainMove::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateMembraneChainMove())
  );
} // namespace bcl
