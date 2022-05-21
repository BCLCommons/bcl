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
#include "descriptor/bcl_descriptor_molecule_rdf_grid_code.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_molecule_rdf_grid_code.cpp
  //!
  //! @author kothiwsk
  //! @date Feb 18, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMoleculeRDFGridCode :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMoleculeRDFGridCode *Clone() const
    {
      return new ExampleDescriptorMoleculeRDFGridCode( *this);
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
      descriptor::CheminfoProperty atom_property( "Atom_SigmaCharge");

      // prepare atom properties vector
      descriptor::CheminfoProperty weight_property( "Atom_VDWaalsRadius");

      // test constructor
      descriptor::MoleculeRDFGridCode rdf_code( atom_property, weight_property);

    /////////////////
    // data access //
    /////////////////

      //! check class identifier
      BCL_ExampleCheck( rdf_code.GetClassIdentifier(), "bcl::descriptor::MoleculeRDFGridCode");

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
        12.2343, 4.50075, 0.224079, 0.00150983, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        9.43998, 5.52999, 0.438418, 0.00470396, 0, 0, 0, 0, 0, 0, 0, 0,
        17.34, 6.37903, 0.317593, 0.00213993, 0, 0, 0, 0, 0, 0, 0, 0,
        18.8787, 11.0592, 0.876776, 0.00940728, 0, 0, 0, 0, 0, 0, 0, 0,
        24.6287, 9.06067, 0.451127, 0.0030399, 0, 0, 0, 0, 0, 0, 0, 0,
        8.51118, 3.13109, 0.155888, 0.00105037, 0, 0, 0, 0, 0, 0, 0, 0,
        18.88, 11.06, 0.876837, 0.00940793, 0, 0, 0, 0, 0, 0, 0, 0,
        9.50464, 5.55377, 0.439601, 0.00471193, 0, 0, 0, 0, 0, 0, 0, 0,
        7.0639, 2.59866, 0.12938, 0.000871755, 0, 0, 0, 0, 0, 0, 0, 0,
        3.5643, 1.31123, 0.0652824, 0.00043987, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0
      };
      linal::Vector< float> expected_vector( 288, expected_vec_array);
      BCL_ExampleCheckWithinAbsTolerance( vector_rdf_code, expected_vector, 0.001);

      BCL_MessageStd( "vector: " + util::Format()( vector_rdf_code));
      BCL_ExampleCheck( vector_rdf_code.GetSize(), 288);

    //////////////////////
    // input and output //
    //////////////////////

      // check i/o
      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( rdf_code, descriptor::MoleculeRDFGridCode()),
        true,
        "MoleculeRDFGridCode I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
      } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMoleculeRDFGridCode

  const ExampleClass::EnumType ExampleDescriptorMoleculeRDFGridCode::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMoleculeRDFGridCode())
  );

} // namespace bcl
