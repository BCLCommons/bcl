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
#include "descriptor/bcl_descriptor_molecule_rdf_sign_code.h"

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
  //! @example example_descriptor_molecule_rdf_sign_code.cpp
  //!
  //! @author kothiwsk
  //! @date Feb 13, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMoleculeRDFSignCode :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMoleculeRDFSignCode *Clone() const
    {
      return new ExampleDescriptorMoleculeRDFSignCode( *this);
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
      descriptor::MoleculeRDFSignCode rdf_code( atom_property, 48, 0.25, 100);

      // test copy constructor
      descriptor::MoleculeRDFSignCode copy_rdf_code( rdf_code);

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
            2.20409,1.90631,0,0.0042549,0.00368004,0,3.06103e-11,0,1.63544e-11,0,0.00010303,0.00227194,0,0.0533696,
            1.18488,0.0248787,0.000431703,2.15183,0.0572205,0.367496,2.21963,0.000387334,0.0444597,0.0320323,0.000606765,
            0.66041,0.377214,1.1823,0.241678,0.393226,0.805799,1.27218,2.07385,0.678119,0.711604,1.88277,1.19085,0.478695,
            1.15817,1.02494,0.478528,1.08129,0.149945,0.458669,1.16699,0.675534,1.09177,1.84119,0.598512,0.874313,1.40868,
            1.02294,1.12279,1.42811,0.813547,1.17203,1.9206,0.937848,0.875067,2.11845,0.736935,1.03009,2.45778,0.754679,
            0.545546,2.21944,1.32468,0.963111,2.23841,1.13825,0.915389,1.87034,1.321,1.04658,2.07909,0.889537,1.58181,
            1.64861,1.78773,1.02135,3.27033,1.1221,1.27881,2.36116,1.09112,0.928843,2.84404,0.779536,1.04524,2.08738,
            0.934572,1.52888,2.09111,1.24244,0.700473,1.92948,1.45685,0.759479,1.03026,0.793692,1.12335,2.18402,0.654879,
            0.766947,1.54008,0.562814,0.471559,1.28148,0.477052,0.783498,2.1452,0.927419,0.428924,1.46028,0.96015,0.58842,
            1.30645,0.814317,1.1387,1.37236,0.686128,0.924154,1.63951,0.562927,0.744697,1.00789,0.633766,0.364225,1.5854,
            0.781471,0.608953,1.02413,0.831184,0.560266,1.54819,0.224578,0.199341,0.997006,0.209363,0.470629,0.827748,
          0.570979, 0.716419, 0.739803
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
        TestBCLObjectIOForSymmetry( rdf_code, descriptor::MoleculeRDFSignCode( atom_property)),
        true,
        "MoleculeRDFSignCode I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
      } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMoleculeRDFSignCode

  const ExampleClass::EnumType ExampleDescriptorMoleculeRDFSignCode::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMoleculeRDFSignCode())
  );

} // namespace bcl
