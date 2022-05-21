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
#include "descriptor/bcl_descriptor_molecule_3da_soft_max.h"

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
  //! @example example_descriptor_molecule_3da_soft_max.cpp
  //!
  //! @author kothiwsk
  //! @date Feb 12, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMolecule3DASoftMax :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMolecule3DASoftMax *Clone() const
    {
      return new ExampleDescriptorMolecule3DASoftMax( *this);
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
      descriptor::Molecule3DASoftMax rdf_code( atom_property, 48, 0.25, 100);

      // test copy constructor
      descriptor::Molecule3DASoftMax copy_rdf_code( rdf_code);

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
      BCL_ExampleCheck( vector_vdw_rdf_code.GetSize(), 48);
      float expected_vec_array [] =
      {
        0.49324, -0.463244, -0.463244, 0.00662011, 0.00665601, 0.0121125, 0.135729, 0.0975556, 0.225363, 0.365713,
        0.243088, 0.363469, 0.310162, 0.311923, 0.374412, 0.481213, 0.315787, 0.445203, 0.30522, 0.249793, 0.242471,
        0.182729, 0.306684, 0.31416, 0.467568, 0.314967, 0.319833, 0.454967, 0.321908, 0.458033, 0.273481, 0.352881,
        0.481265, 0.363929, 0.310144, 0.3024, 0.314344, 0.283001, 0.46761, 0.0900904, 0.311684, 0.470658, 0.237527,
        0.302947, 0.318666, 0.151022, 0.259778, 0.311107
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
        TestBCLObjectIOForSymmetry( rdf_code, descriptor::Molecule3DASoftMax( atom_property)),
        true,
        "Molecule3DASoftMax I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
      } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMolecule3DASoftMax

  const ExampleClass::EnumType ExampleDescriptorMolecule3DASoftMax::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMolecule3DASoftMax())
  );

} // namespace bcl
