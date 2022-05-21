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
#include "descriptor/bcl_descriptor_molecule_3da_smooth_sign_code.h"

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
  //! @example example_descriptor_molecule_3da_smooth_sign_code.cpp
  //!
  //! @author kothiwsk
  //! @date Feb 11, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMolecule3DASmoothSignCode :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMolecule3DASmoothSignCode *Clone() const
    {
      return new ExampleDescriptorMolecule3DASmoothSignCode( *this);
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
      descriptor::Molecule3DASmoothSignCode rdf_code( atom_property, 48, 0.25, 100);

      // test copy constructor
      descriptor::Molecule3DASmoothSignCode copy_rdf_code( rdf_code);

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
      BCL_ExampleCheck( vector_vdw_rdf_code.GetSize(), 144);
      float expected_vec_array [] =
      {
        2.20409, 1.90631, 5.82135e-28, 0.0042549, 0.00368004, 2.17133e-14, 3.06103e-11, 8.76898e-08, 3.01821e-06, 7.875e-13,
        0.000148239, 0.00384681, 0.000109464, 0.0532669, 1.18716, 0.0569446, 0.000970815, 2.26568, 0.124826, 0.41032, 2.90415,
        0.000652868, 0.088114, 0.0598109, 0.00310012, 0.870877, 0.528153, 1.60753, 0.432763, 0.731848, 1.29513, 1.65335, 2.97355,
        1.06318, 1.07789, 2.39382, 1.49988, 0.716546, 1.83774, 1.10215, 0.812232, 1.46968, 0.251004, 0.686746, 1.64051, 0.942419,
        1.61213, 2.47361, 1.00619, 1.25213, 1.95793, 1.36967, 1.6411, 2.06257, 1.23846, 1.75564, 2.81203, 1.31748, 1.36852, 3.04949,
        1.12343, 1.39878, 3.87069, 1.22876, 0.834661, 3.18182, 1.71993, 1.4801, 3.26965, 1.62962, 1.48951, 2.7757, 2.09884, 1.57746,
        2.68985, 1.46109, 2.48374, 2.27964, 2.11092, 1.42978, 4.72589, 1.65906, 1.77006, 3.5284, 1.6299, 1.21721, 3.99162, 1.1503,
        1.35671, 2.75748, 1.269, 1.98873, 2.78251, 1.57121, 1.07963, 2.52724, 1.96458, 1.13978, 1.52269, 1.46104, 1.66206, 3.03516,
        1.11531, 0.947067, 2.24839, 0.769896, 0.688474, 1.91421, 0.658805, 1.1309, 2.87706, 1.15995, 0.595443, 2.16547, 1.12039,
        0.861634, 1.7122, 1.13394, 1.74964, 1.91447, 0.98431, 1.32069, 2.21149, 0.926004, 1.09811, 1.44337, 1.01091, 0.489723,
        2.11088, 1.00033, 0.853371, 1.49083, 1.09888, 0.862122, 2.20216, 0.293676, 0.277373, 1.5041, 0.342873, 0.696823, 1.32455,
        0.698578, 0.870644, 0.87035
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
        TestBCLObjectIOForSymmetry( rdf_code, descriptor::Molecule3DASmoothSignCode( atom_property)),
        true,
        "Molecule3DASmoothSignCode I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
      } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMolecule3DASmoothSignCode

  const ExampleClass::EnumType ExampleDescriptorMolecule3DASmoothSignCode::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMolecule3DASmoothSignCode())
  );

} // namespace bcl
