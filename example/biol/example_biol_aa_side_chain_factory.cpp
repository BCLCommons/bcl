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
#include "biol/bcl_biol_aa_side_chain_factory.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "quality/bcl_quality_rmsd.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_aa_side_chain_factory.cpp
  //!
  //! @author karakam
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolAASideChainFactory :
    public ExampleInterface
  {
  public:

      ExampleBiolAASideChainFactory *Clone() const
    {
      return new ExampleBiolAASideChainFactory( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      BCL_MessageStd( "Testing default constructor");

      // create object from constructor that takes two bools and a string
      biol::AASideChainFactory this_object( false, true);

      // create copy of this_object
      biol::AASideChainFactory object_copy( this_object);

      // create object from clone constructor
      util::ShPtr< util::ObjectInterface> object_clone( this_object.Clone());

    /////////////////
    // data access //
    /////////////////

      // check the class identifier
      BCL_Example_Check
      (
        GetStaticClassName( this_object) == "bcl::biol::AASideChainFactory",
        "unexpected static class name: " + GetStaticClassName( this_object)
          + " should be: bcl::biol::AASideChainFactory"
      );

      // class identifier
      BCL_Example_Check
      (
        object_clone->GetClassIdentifier() == GetStaticClassName< biol::AASideChainFactory>(),
        "unexpected class identifier class name: " + object_clone->GetClassIdentifier() + " should be: "
        + GetStaticClassName< biol::AASideChainFactory>()
      );

    ////////////////
    // operations //
    ////////////////

      // test adding side chains to a protein from pdb with output to a new pdb
      const std::string pdb( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
      assemble::ProteinModel this_model( Proteins::GetModel( pdb, biol::GetAAClasses().e_AAComplete));
      util::ShPtr< assemble::ProteinModel> model_ptr( this_object.ProteinModelWithSideChains( this_model));

      //write models which have only the missing loops and only the BackBone coordinates
      BCL_MessageStd( "writing 1ubi_with_side_chains.pdb");
      const std::string outfile( AddExampleOutputPathToFilename( this_object, "1ubi_with_side_chains.pdb"));
      Proteins::WriteModelToPDB( *model_ptr, outfile);

      // get the SSEs in original protein model
      util::SiPtrVector< const assemble::SSE> some_sses_a( this_model.GetSSEs());

      // get the atom coordinates for the first SSE in the vector of SSEs of original model
      util::SiPtrVector< const linal::Vector3D> atom_coordinates_a( some_sses_a( 0)->GetAtomCoordinates());

      // get the SSEs in new model with side chains
      util::SiPtrVector< const assemble::SSE> some_sses_b( model_ptr->GetSSEs());

      // get the atom coordinates for the first SSE in the vector of SSEs in new model
      util::SiPtrVector< const linal::Vector3D> atom_coordinates_b( some_sses_b( 0)->GetAtomCoordinates());

      // Calculate RMSD between model with side chains and old model from pdb
      quality::RMSD rmsd( false);
      BCL_MessageStd( "Calculate RMSD");
      const double calculated_rmsd( rmsd.CalculateMeasure( atom_coordinates_a, atom_coordinates_b));
      const double correct_rmsd( 1.70398);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( correct_rmsd, calculated_rmsd),
        "The calculated RMSD is " + util::Format()( calculated_rmsd) + " but should be " + util::Format()( correct_rmsd)
      );

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBiolAASideChainFactory

  const ExampleClass::EnumType ExampleBiolAASideChainFactory::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolAASideChainFactory())
  );

} // namespace bcl
