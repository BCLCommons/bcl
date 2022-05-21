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
#include "descriptor/bcl_descriptor_molecule_2da_code.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_molecule_2da_code.cpp
  //!
  //! @author kothiwsk
  //! @date Feb 11, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMolecule2DACode :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMolecule2DACode *Clone() const
    {
      return new ExampleDescriptorMolecule2DACode( *this);
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

      // test constructor
      descriptor::Molecule2DACode two_d_a_code( 12, atom_property);

    /////////////////
    // data access //
    /////////////////

      //! check class identifier
      BCL_ExampleCheck( two_d_a_code.GetClassIdentifier(), "bcl::descriptor::Molecule2DACode");

      BCL_MessageStd
      (
        "Charges: "
        + util::Format()( atom_property->CollectValuesOnEachElementOfObject( small_mol))
      );

      two_d_a_code.SetDimension( 0);
      BCL_ExampleIndirectCheck( two_d_a_code.GetNumberSteps(), 12, "SetNumberSteps");

      // Generate the  2da code and store as storage::Vector
      linal::Vector< float> vector_2da_code( two_d_a_code.SumOverObject( small_mol));

      float expected_vec_array [] =
      {
          0.0405826,0,-0.0202913,-0.0101457,4.65661e-10,0.0101457,0,0,0,0,0,0
      };
      linal::Vector< float> expected_vector( 12, expected_vec_array);
      BCL_ExampleCheckWithinAbsTolerance( vector_2da_code, expected_vector, 0.001);

      BCL_MessageStd( "vector: " + util::Format()( vector_2da_code));

      // Write object fire for 2da code
      WriteBCLObject( two_d_a_code);

    //////////////////////
    // input and output //
    //////////////////////

      // check i/o
      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( two_d_a_code, descriptor::Molecule2DACode()),
        true,
        "Molecule2DACode I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
      } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMolecule2DACode

  const ExampleClass::EnumType ExampleDescriptorMolecule2DACode::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMolecule2DACode())
  );

} // namespace bcl
