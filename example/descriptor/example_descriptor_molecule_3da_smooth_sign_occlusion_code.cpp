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
#include "descriptor/bcl_descriptor_molecule_3da_smooth_sign_occlusion_code.h"

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
  //! @example example_descriptor_molecule_3da_smooth_sign_occlusion_code.cpp
  //!
  //! @author mendenjl
  //! @date Apr 22, 2016
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMolecule3DASmoothSignOcclusionCode :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMolecule3DASmoothSignOcclusionCode *Clone() const
    {
      return new ExampleDescriptorMolecule3DASmoothSignOcclusionCode( *this);
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
      descriptor::Molecule3DASmoothSignOcclusionCode rdf_code( atom_property, 12, 0.5);

      // test copy constructor
      descriptor::Molecule3DASmoothSignOcclusionCode copy_rdf_code( rdf_code);

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
      BCL_ExampleCheck( rdf_code.GetNumberSteps(), 12);
      BCL_ExampleCheck( rdf_code.GetStepSize(), float( 0.5));

      rdf_code.SetDimension( 0);

      // Generate the VDW RDF code
      const linal::Vector< float> vector_vdw_rdf_code( rdf_code.SumOverObject( small_mol));
      BCL_ExampleCheck( vector_vdw_rdf_code.GetSize(), 12 * 3 * 3);
      float expected_vec_array [] =
      {
        2.20409, 1.90631, 0, 0, 1.48155e-06, 3.79602e-05, 0.038678, 0.0535076, 2.60632, 0.142506, 0.44659, 3.76038,
        0.737572, 1.0903, 0.992962, 2.71968, 2.48177, 4.31955, 2.51116, 1.58202, 3.6196, 0.906828, 1.44751, 2.69353,
        1.40728, 2.22933, 2.63373, 1.83586, 2.01221, 3.82129, 1.65557, 1.7147, 5.26411, 2.23837, 1.81111, 4.14697,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00126188, 0, 0, 0.0648108, 0.157838, 0.419437, 0.227004,
        0.743976, 0.378204, 0.343803, 0.991617, 0.487298, 0.531608, 0.821204, 0.504609, 0.410683, 1.0056, 0.578471,
        0.434584, 1.3204, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00660693, 0.00690304, 0, 0.0701456,
        0.215175, 0.256431, 0.361911, 0.716488, 0.293093, 0.687436, 0.68397, 0.267067, 0.386552, 0.562506, 0.270552,
        0.35212, 0.865678
      };
      linal::Vector< float> expected_vector( 108, expected_vec_array);
      BCL_ExampleCheckWithinAbsTolerance( vector_vdw_rdf_code, expected_vector, 0.001);

      BCL_MessageStd( "vector: " + util::Format()( vector_vdw_rdf_code));

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
      } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMolecule3DASmoothSignOcclusionCode

  const ExampleClass::EnumType ExampleDescriptorMolecule3DASmoothSignOcclusionCode::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMolecule3DASmoothSignOcclusionCode())
  );

} // namespace bcl
