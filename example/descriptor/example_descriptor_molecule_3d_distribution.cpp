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
#include "descriptor/bcl_descriptor_molecule_3d_distribution.h"

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
  //! @example example_descriptor_molecule_3d_distribution.cpp
  //!
  //! @author mendenjl
  //! @date Aug 09, 2019
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMolecule3DDistribution :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMolecule3DDistribution *Clone() const
    {
      return new ExampleDescriptorMolecule3DDistribution( *this);
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
      descriptor::CheminfoProperty atom_property_center( "Abs(Atom_Vcharge)");

      // test empty constructor
      descriptor::Molecule3DDistribution rdf_code( atom_property, atom_property_center, 48, 0.25);

      // test copy constructor
      descriptor::Molecule3DDistribution copy_rdf_code( rdf_code);

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

      rdf_code.SetDimension( 0);

      // Generate the VDW RDF code
      const linal::Vector< float> vector_vdw_rdf_code( rdf_code.SumOverObject( small_mol));
      BCL_ExampleCheck( vector_vdw_rdf_code.GetSize(), 48);
      float expected_vec_array [] =
      {
        0, 0, -0.0204896, -0.0455322, 0.239821, 0.0118741, -0.0369431, -0.0703164, -0.0982207, -0.263836, 0.247425, 0.230885,
        0.23881, 0.35325, -0.731441, -0.584265, -0.135815, 0.205715, 0.587552, 0.794961, 0.708338, -0.643542, -0.988858,
        -0.278111, -0.159551, 0.147028, 0.222885, 0.348159, 0.263035, -0.0309518, -0.087607, -0.186346, -0.402818, 0.0703569,
        -0.109401, 0.0183591, 0.0334466, -0.0189894, -0.0386644, 0.0365314, -0.0489164, 0.0928426, 0.00210669, -0.0903002,
        0.0674334, 0.0413337, 0.026791, 0.0819761
      };
      linal::Vector< float> expected_vector( 48, expected_vec_array);
      BCL_ExampleCheckWithinAbsTolerance( vector_vdw_rdf_code, expected_vector, 0.001);

      WriteBCLObject( rdf_code);

    //////////////////////
    // input and output //
    //////////////////////

      // check i/o
      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( rdf_code, descriptor::Molecule3DDistribution( atom_property, atom_property_center)),
        true,
        "Molecule3DDistribution I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
      } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMolecule3DDistribution

  const ExampleClass::EnumType ExampleDescriptorMolecule3DDistribution::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMolecule3DDistribution())
  );

} // namespace bcl
