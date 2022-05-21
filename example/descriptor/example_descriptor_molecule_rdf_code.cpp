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
#include "descriptor/bcl_descriptor_molecule_rdf_code.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_molecule_rdf_code.cpp
  //!
  //! @author kothiwsk
  //! @date Feb 12, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMoleculeRDFCode :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMoleculeRDFCode *Clone() const
    {
      return new ExampleDescriptorMoleculeRDFCode( *this);
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
      descriptor::MoleculeRDFCode rdf_code( atom_property, 128, 0.1, 100);

    /////////////////
    // data access //
    /////////////////

      //! check class identifier
      BCL_ExampleCheck( rdf_code.GetClassIdentifier(), "bcl::descriptor::MoleculeRDFCode");

      BCL_MessageStd
      (
        "Charges: "
        + util::Format()( atom_property->CollectValuesOnEachElementOfObject( small_mol))
      );

      rdf_code.SetDimension( 0);
      BCL_ExampleIndirectCheck( rdf_code.GetNumberSteps(), 128, "SetNumberSteps");
      BCL_ExampleCheck( rdf_code.GetStepSize(), float( 0.1));
      BCL_ExampleCheck( rdf_code.GetTemperature(), 100.0);

      // Generate the  2da code and store as storage::Vector
      linal::Vector< float> vector_rdf_code( rdf_code.SumOverObject( small_mol));

      float expected_vec_array [] =
      {
          0,-1.12104e-44,-1.66599e-36,-3.03592e-29,-7.48746e-23,-2.49926e-17,-1.12907e-12,-6.90327e-09,-5.71245e-06,
          -0.000639765,-0.00969722,-0.0198898,-0.00506874,0.00803632,0.0202374,0.00672454,0.000302404,1.84052e-06,
          -1.44935e-07,-6.40078e-05,-0.0037863,-0.0303125,-0.0326485,0.00123154,0.0302361,0.0295901,0.00586929,0.00470133,
          0.0100417, 0.00301305, 0.00012236, -4.3e-06, -0.000739705, -0.0148928, -0.040582, -0.0149667, -0.000757173, -0.00092977,
          -0.0114331,-0.0191315,-0.00432977,0.000286985,0.00792924,0.0202722,0.00701385,0.000328421,2.10241e-06,1.07825e-05,
          0.000743314,0.00693606,0.00875934,0.00149709,3.46283e-05,1.08402e-07,4.59265e-11,2.63333e-15,2.04346e-20,
          2.14592e-26,3.05003e-33,5.86682e-41,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
      };
      linal::Vector< float> expected_vector( 128, expected_vec_array);
      BCL_ExampleCheckWithinAbsTolerance( vector_rdf_code, expected_vector, 0.001);

      BCL_MessageStd( "vector: " + util::Format()( vector_rdf_code));
      BCL_ExampleCheck( vector_rdf_code.GetSize(), 128);

      // Write object fire for 2da code
      WriteBCLObject( rdf_code);

    //////////////////////
    // input and output //
    //////////////////////

      // check i/o
      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( rdf_code, descriptor::MoleculeRDFCode()),
        true,
        "MoleculeRDFCode I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
      } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMoleculeRDFCode

  const ExampleClass::EnumType ExampleDescriptorMoleculeRDFCode::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMoleculeRDFCode())
  );

} // namespace bcl
