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
#include "descriptor/bcl_descriptor_molecule_2da_smooth_sign_code.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_molecule_2da_smooth_sign_code.cpp
  //!
  //! @author raftersa, mendenjl
  //! @date Mar 13, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMolecule2DASmoothSignCode :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMolecule2DASmoothSignCode *Clone() const
    {
      return new ExampleDescriptorMolecule2DASmoothSignCode( *this);
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
      BCL_ExampleMustOpenInputFile( input_sdf, AddExampleInputPathToFilename( e_Chemistry, "C2H7BO3.sdf"));

      // load information into small_mol
      chemistry::FragmentComplete small_mol = sdf::FragmentFactory::MakeFragment( input_sdf);
      io::File::CloseClearFStream( input_sdf);

      // prepare atom properties vector
      descriptor::CheminfoProperty atom_property( "Atom_SigmaCharge");

      // test constructor
      descriptor::Molecule2DASmoothSignCode two_d_a_code( atom_property, 12, 100, descriptor::Molecule2DASmoothSignCode::e_2DASign);

    /////////////////
    // data access //
    /////////////////

      //! check class identifier
      BCL_ExampleCheck( two_d_a_code.GetClassIdentifier(), "bcl::descriptor::Molecule2DASmoothSignCode");

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
          0.455865,0.28694,1.4013e-44,2.10195e-44,0.00297883,0.38281,0.122982,0.139022,0.0302188,
          1.12104e-44,0.0492185,0.277245,0.199635,0.109105,0.05896,7.00649e-45,0.0286025,0.19991,
          0,0.0499979,7.00649e-45,0,1.4013e-45,0,0,0,0,0,0,0,0,0,0,0,0,0
      };
      linal::Vector< float> expected_vector( 36, expected_vec_array);
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
        TestBCLObjectIOForSymmetry
        (
          two_d_a_code,
          descriptor::Molecule2DASmoothSignCode( descriptor::Molecule2DASmoothSignCode::e_2DASign)
        ),
        true,
        "Molecule2DASmoothSignCode I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
      } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMolecule2DASmoothSignCode

  const ExampleClass::EnumType ExampleDescriptorMolecule2DASmoothSignCode::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMolecule2DASmoothSignCode())
  );

} // namespace bcl
