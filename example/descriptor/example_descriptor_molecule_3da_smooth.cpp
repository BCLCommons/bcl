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
#include "descriptor/bcl_descriptor_molecule_3da_smooth.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "descriptor/bcl_descriptor_atom_vcharge.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_molecule_3da_smooth.cpp
  //!
  //! @author kothiwsk
  //! @date Feb 12, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMolecule3DASmooth :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMolecule3DASmooth *Clone() const
    {
      return new ExampleDescriptorMolecule3DASmooth( *this);
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
    /////////////////
    // preparation //
    /////////////////

      // Read in small molecule
      io::IFStream input_sdf;
      BCL_ExampleMustOpenInputFile( input_sdf, AddExampleInputPathToFilename( e_Chemistry, "taxol.sdf"));

      // load information into small_mol
      chemistry::FragmentComplete small_mol( sdf::FragmentFactory::MakeFragment( input_sdf));

      io::File::CloseClearFStream( input_sdf);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // prepare atom properties vector
      descriptor::CheminfoProperty atom_property( "Atom_Vcharge");

      // test empty constructor
      descriptor::Molecule3DASmooth rdf_code( atom_property, 48, 0.25, 100);

      // test copy constructor
      descriptor::Molecule3DASmooth copy_rdf_code( rdf_code);

      BCL_ExampleCheck( copy_rdf_code.GetAtomProperty().GetString(), "Atom_Vcharge");

    /////////////////
    // data access //
    /////////////////

      // Check all member variables for accuracy
      BCL_ExampleCheck
      (
        rdf_code.GetAtomProperty().GetAlias(),
        descriptor::AtomVcharge().GetAlias()
      );
      BCL_ExampleCheck( rdf_code.GetNumberSteps(), 48);
      BCL_ExampleCheck( rdf_code.GetStepSize(), float( 0.25));
      BCL_ExampleCheck( rdf_code.GetTemperature(), 100.0);

      rdf_code.SetDimension( 0);

      // Generate the VDW RDF code
      const linal::Vector< float> vector_vdw_rdf_code( rdf_code.SumOverObject( small_mol));
      BCL_ExampleCheck( vector_vdw_rdf_code.GetSize(), 48);
      float expected_vec_array [] =
      {
        4.1104, 0.00793494, -2.93049e-06, -0.00369857, -1.13379, -2.20777, -2.369, 0.0289561, 0.345824, 1.30845,
        -0.0250681, -0.252743, 0.378685, 0.444698, -0.702761, 0.0809379, 0.300381, 0.948199, 0.18207, -0.363488,
        -1.34847, -1.1184, -0.0696177, 0.343426, 0.986443, 1.66519, -1.18518, -0.0992696, -1.14451, -0.250477, 0.475224,
        0.123601, 1.58167, 0.0879383, -0.186015, -0.455844, -1.08736, -0.410073, 0.269828, 0.969113, 0.0935092,
        0.580742, -0.610249, 0.36287, -0.241158, -0.933049, -0.284854, 0.698873
      };
      linal::Vector< float> expected_vector( 48, expected_vec_array);
      BCL_ExampleCheckWithinAbsTolerance( vector_vdw_rdf_code, expected_vector, 0.001);

      BCL_MessageStd( "vector: " + util::Format()( vector_vdw_rdf_code));

      WriteBCLObject( rdf_code);

    //////////////////////
    // input and output //
    //////////////////////

      // check i/o
      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( rdf_code, descriptor::Molecule3DASmooth( atom_property)),
        true,
        "Molecule3DASmooth I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
      } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMolecule3DASmooth

  const ExampleClass::EnumType ExampleDescriptorMolecule3DASmooth::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMolecule3DASmooth())
  );

} // namespace bcl
