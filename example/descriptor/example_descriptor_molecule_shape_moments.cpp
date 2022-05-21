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
#include "descriptor/bcl_descriptor_molecule_shape_moments.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_molecule_shape_moments.cpp
  //!
  //! @author mendenjl
  //! @date Mar 09, 2015
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMoleculeShapeMoments :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMoleculeShapeMoments *Clone() const
    {
      return new ExampleDescriptorMoleculeShapeMoments( *this);
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
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // Read in molecule
      io::IFStream input_sdf;
      BCL_ExampleMustOpenInputFile( input_sdf, AddExampleInputPathToFilename( e_Chemistry, "benzene_moe_bcl_out.sdf"));

      // load information into small_mol
      chemistry::FragmentComplete small_mol = sdf::FragmentFactory::MakeFragment( input_sdf);
      io::File::CloseClearFStream( input_sdf);

      // prepare atom properties vector
      descriptor::CheminfoProperty atom_property( "Abs(Atom_SigmaCharge)");

      // test constructor
      descriptor::MoleculeShapeMoments rdf_code( atom_property, atom_property);

    /////////////////
    // data access //
    /////////////////

      //! check class identifier
      BCL_ExampleCheck( rdf_code.GetClassIdentifier(), "bcl::descriptor::MoleculeShapeMoments");

      BCL_MessageStd
      (
        "Charges: "
        + util::Format()( atom_property->CollectValuesOnEachElementOfObject( small_mol))
      );

      rdf_code.SetDimension( 0);

      // Generate the  2da code and store as storage::Vector
      linal::Vector< float> vector_rdf_code( rdf_code.SumOverObject( small_mol));

      float expected_vec_array [] =
      {
        1.93787, 0.641455, 0.0,      2.20665, 1.25347, -0.274771,
        2.88246, 1.62729, -0.373003, 2.88213, 1.62735, -0.372777
      };
      linal::Vector< float> expected_vector( 12, expected_vec_array);
      BCL_ExampleCheckWithinAbsTolerance( vector_rdf_code, expected_vector, 0.001);

      BCL_MessageStd( "vector: " + util::Format()( vector_rdf_code));
      BCL_ExampleCheck( vector_rdf_code.GetSize(), 12);

      // Write object fire for 2da code
      WriteBCLObject( rdf_code);

    //////////////////////
    // input and output //
    //////////////////////

      // check i/o
      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( rdf_code, descriptor::MoleculeShapeMoments()),
        true,
        "MoleculeShapeMoments I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
      } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMoleculeShapeMoments

  const ExampleClass::EnumType ExampleDescriptorMoleculeShapeMoments::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMoleculeShapeMoments())
  );

} // namespace bcl
