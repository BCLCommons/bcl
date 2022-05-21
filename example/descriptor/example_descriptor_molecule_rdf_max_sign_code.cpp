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
#include "descriptor/bcl_descriptor_molecule_rdf_max_sign_code.h"

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
  //! @example example_descriptor_molecule_rdf_max_sign_code.cpp
  //!
  //! @author kothiwsk
  //! @date Feb 13, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMoleculeRDFMaxSignCode :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMoleculeRDFMaxSignCode *Clone() const
    {
      return new ExampleDescriptorMoleculeRDFMaxSignCode( *this);
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
      descriptor::MoleculeRDFMaxSignCode rdf_code( atom_property, 48, 0.25, 100);

      // test copy constructor
      descriptor::MoleculeRDFMaxSignCode copy_rdf_code( rdf_code);

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
            0.217536,0.24662,0,0.000419943,0.000476089,0,0,0,0,0,1.27651e-05,0.000427884,3.14333e-09,0.0066204,0.221432,
          0.00166831,0.000166119,0.381941,0.00476596,0.134444,0.378399,0.000301461,0.0224303,0.00987728,0.000450618,
            0.109057,0.096406,0.320257,0.0820786,0.067193,0.180689,0.257045,0.111479,0.138833,0.072364,0.325261,0.309336,
            0.155437,0.247711,0.305216,0.0686165,0.132608,0.0318441,0.119621,0.195776,0.234582,0.292465,0.194417,0.135161,
            0.271454,0.352673,0.232581,0.202716,0.246859,0.234459,0.0798443,0.258571,0.275106,0.0970293,0.142898,0.0604441,
            0.106929,0.283127,0.115704,0.0740058,0.27122,0.290794,0.141939,0.26128,0.294745,0.191366,0.205961,0.250353,
            0.242662,0.368423,0.124583,0.175643,0.37229,0.318329,0.163743,0.367409,0.197599,0.361858,0.360316,0.257286,
            0.240752,0.351964,0.158375,0.255647,0.43439,0.226299,0.45522,0.362278,0.326355,0.0777621,0.327598,0.291577,
            0.0971583,0.0881306,0.162262,0.235276,0.295953,0.112741,0.285231,0.257146,0.28054,0.0867511,0.250906,0.120675,
            0.162176,0.432062,0.301646,0.129258,0.375354,0.297682,0.0655743,0.330184,0.281628,0.172203,0.375054,0.288392,
            0.175328,0.307113,0.224725,0.343976,0.109448,0.179459,0.0758988,0.383062,0.287009,0.11805,0.194181,0.264815,
            0.0776482,0.292432,0.0445708,0.0619357,0.196035,0.0568433,0.128195,0.113888,0.194118,0.226438,0.0728552
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
        TestBCLObjectIOForSymmetry( rdf_code, descriptor::MoleculeRDFMaxSignCode( atom_property)),
        true,
        "MoleculeRDFMaxSignCode I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
      } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMoleculeRDFMaxSignCode

  const ExampleClass::EnumType ExampleDescriptorMoleculeRDFMaxSignCode::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMoleculeRDFMaxSignCode())
  );

} // namespace bcl
