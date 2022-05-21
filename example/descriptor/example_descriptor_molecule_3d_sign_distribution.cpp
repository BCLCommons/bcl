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
#include "descriptor/bcl_descriptor_molecule_3d_sign_distribution.h"

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
  //! @example example_descriptor_molecule_3d_sign_distribution.cpp
  //!
  //! @author mendenjl
  //! @date Aug 09, 2019
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMolecule3DSignDistribution :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMolecule3DSignDistribution *Clone() const
    {
      return new ExampleDescriptorMolecule3DSignDistribution( *this);
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
      descriptor::CheminfoProperty atom_center_property( "Abs(Atom_Vcharge)");

      // test empty constructor
      descriptor::Molecule3DSignDistribution rdf_code( atom_property, atom_center_property, 48, 0.25);

      // test copy constructor
      descriptor::Molecule3DSignDistribution copy_rdf_code( rdf_code);

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
      BCL_ExampleCheck( vector_vdw_rdf_code.GetSize(), 96);
      float expected_vec_array [] =
      {
        0, 0, 0, 0, 0.239821, 0.0118741, 0.0168901, 0.00936127, 0, 0.139597, 0.301595, 0.230885,
        0.338258, 0.96017, 0.272155, 0.0840725, 0.637054, 0.406149, 0.587552, 0.824075, 0.814973,
        0.166545, 0.230296, 0.273291, 0.271665, 0.328601, 0.500552, 0.689724, 0.723108, 0.125334,
        0.0637381, 0.228564, 0.121376, 0.0742193, 0.0485876, 0.0765882, 0.231122, 0.054062,
        0.0147095, 0.0940576, 0, 0.161518, 0.0560158, 0, 0.0674334, 0.0413337, 0.026791, 0.0819761,
        0, 0, 0.0204896, 0.0455322, 0, 0, 0.0538331, 0.0796776, 0.0982207, 0.403433, 0.0541705, 0,
        0.0994478, 0.606919, 1.0036, 0.668337, 0.772869, 0.200434, 0, 0.0291144, 0.106635, 0.810087,
        1.21915, 0.551402, 0.431216, 0.181572, 0.277667, 0.341565, 0.460072, 0.156286, 0.151345,
        0.41491, 0.524193, 0.00386239, 0.157989, 0.0582291, 0.197675, 0.0730514, 0.053374, 0.0575262,
        0.0489164, 0.0686758, 0.0539091, 0.0903002, 0, 0, 0, 0
      };
      linal::Vector< float> expected_vector( 96, expected_vec_array);
      BCL_ExampleCheckWithinAbsTolerance( vector_vdw_rdf_code, expected_vector, 0.001);

    //////////////////////
    // input and output //
    //////////////////////

      // check i/o
      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( rdf_code, descriptor::Molecule3DSignDistribution( atom_property, atom_center_property)),
        true,
        "Molecule3DSignDistribution I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
      } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMolecule3DSignDistribution

  const ExampleClass::EnumType ExampleDescriptorMolecule3DSignDistribution::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMolecule3DSignDistribution())
  );

} // namespace bcl
