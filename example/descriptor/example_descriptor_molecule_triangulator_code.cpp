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
#include "descriptor/bcl_descriptor_molecule_triangulator_code.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_molecule_triangulator_code.cpp
  //!
  //! @author kothiwsk
  //! @date Feb 18, 2014
  //!
  //! @see @link example_descriptor_molecule_triangulator_code.cpp @endlink
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMoleculeTriangulatorCode :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMoleculeTriangulatorCode *Clone() const
    {
      return new ExampleDescriptorMoleculeTriangulatorCode( *this);
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

      // test empty constructor
      descriptor::MoleculeTriangulatorCode triangulator_code( atom_property, size_t( 12), float( 0.01), float( 0.05), float( 0.0));

    /////////////////
    // data access //
    /////////////////

      BCL_MessageStd
      (
        "Charges: "
        + util::Format()( atom_property->CollectValuesOnEachElementOfObject( small_mol))
      );

      // Check member variables for accuracy
      BCL_ExampleCheck( triangulator_code.GetAtomProperty(), atom_property);
      BCL_ExampleCheck( triangulator_code.GetNumberSteps(), 12);

      linal::Vector< float> vector_3da_code( triangulator_code( small_mol));

      float expected_vec_array [] =
      {
        4.62765e-05,4.33398e-05,4.03947e-05,3.74412e-05,3.44789e-05,3.15063e-05,2.85282e-05,2.55389e-05,2.25471e-05,
        1.95489e-05,1.65431e-05,1.35247e-05
      };
      linal::Vector< float> expected_vector( 12, expected_vec_array);
      BCL_ExampleCheckWithinAbsTolerance( vector_3da_code, expected_vector, 0.001);

    //////////////////////
    // input and output //
    //////////////////////

      // check i/o
      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( triangulator_code, descriptor::MoleculeTriangulatorCode()),
        true,
        "SmallMoleculeTriangulatorCode I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
      } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMoleculeTriangulatorCode

  const ExampleClass::EnumType ExampleDescriptorMoleculeTriangulatorCode::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMoleculeTriangulatorCode())
  );

} // namespace bcl
