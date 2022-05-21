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
#include "descriptor/bcl_descriptor_molecule_3da_soft_max_sign.h"

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
  //! @example example_descriptor_molecule_3da_soft_max_sign.cpp
  //!
  //! @author kothiwsk
  //! @date Feb 12, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMolecule3DASoftMaxSign :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMolecule3DASoftMaxSign *Clone() const
    {
      return new ExampleDescriptorMolecule3DASoftMaxSign( *this);
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
      descriptor::Molecule3DASoftMaxSign rdf_code( atom_property, 48, 0.25, 100);

      // test copy constructor
      descriptor::Molecule3DASoftMaxSign copy_rdf_code( rdf_code);

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
      const linal::Vector< float> vector_vdw_rdf_code( rdf_code( small_mol));
      BCL_ExampleCheck( vector_vdw_rdf_code.GetSize(), 144);
      float expected_vec_array [] =
      {
          0.435072, 0.49324, 8.24473e-26, 0.000839886, 0.000952177, 3.07524e-12, 6.04225e-12, 1.27552e-05, 0.000427465,
          1.64406e-13, 0.00662011, 0.222184, 2.28528e-05, 0.00664118, 0.390687, 0.0118487, 0.00443493, 0.441989,
          0.00553109, 0.135714, 0.101381, 0.00031393, 0.0975556, 0.104289, 0.00070484, 0.225155, 0.0287338, 0.36528,
          0.257858, 0.376031, 0.243088, 0.168205, 0.303917, 0.363469, 0.298115, 0.438851, 0.310162, 0.17659, 0.140678,
          0.311923, 0.1658, 0.143946, 0.373627, 0.169792, 0.398371, 0.0746837, 0.480737, 0.191292, 0.31461, 0.273704,
          0.385211, 0.242281, 0.445046, 0.380022, 0.304828, 0.265705, 0.216155, 0.249793, 0.15064, 0.387017, 0.242471,
          0.166736, 0.37946, 0.182727, 0.169458, 0.369288, 0.305597, 0.306555, 0.194058, 0.313804, 0.26038, 0.377173,
          0.284091, 0.467384, 0.388412, 0.314613, 0.27296, 0.376919, 0.319559, 0.0736262, 0.387345, 0.312924, 0.454385,
          0.388577, 0.321335, 0.266317, 0.451912, 0.302965, 0.457925, 0.240781, 0.244376, 0.273072, 0.379604, 0.35262,
          0.140901, 0.144818, 0.374919, 0.480718, 0.386838, 0.363724, 0.292283, 0.357707, 0.310144, 0.144896, 0.368963,
          0.3024, 0.0751094, 0.464699, 0.314344, 0.264976, 0.387908, 0.282692, 0.143057, 0.379639, 0.307351, 0.46734,
          0.117516, 0.0745493, 0.0900015, 0.376376, 0.311464, 0.265637, 0.376401, 0.371847, 0.470252, 0.130824,
          0.237336, 0.0728111, 0.389431, 0.302947, 0.153283, 0.387526, 0.318504, 0.0796017, 0.256352, 0.0670464,
          0.150561, 0.116173, 0.250872, 0.259332, 0.14919, 0.31109, 0.0799911, 0.219168
      };
      linal::Vector< float> expected_vector( 144, expected_vec_array);
      BCL_ExampleCheckWithinAbsTolerance( vector_vdw_rdf_code, expected_vector, 0.001);

      BCL_MessageStd( "vector: " + util::Format()( vector_vdw_rdf_code));

      WriteBCLObject( rdf_code);

    //////////////////////
    // input and output //
    //////////////////////

      // check i/o
      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( rdf_code, descriptor::Molecule3DASoftMaxSign( atom_property)),
        true,
        "Molecule3DASoftMaxSign I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
      } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMolecule3DASoftMaxSign

  const ExampleClass::EnumType ExampleDescriptorMolecule3DASoftMaxSign::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMolecule3DASoftMaxSign())
  );

} // namespace bcl
