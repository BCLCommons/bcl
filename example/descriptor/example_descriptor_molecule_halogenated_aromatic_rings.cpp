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
#include "descriptor/bcl_descriptor_molecule_halogenated_aromatic_rings.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "descriptor/bcl_descriptor_molecule_misc_property.h"
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_linear_least_squares.h"
#include "math/bcl_math_statistics.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_molecule_halogenated_aromatic_rings.cpp
  //!
  //! @author brownbp1
  //! @date Aug 23, 2021
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMoleculeHalogenatedAromaticRings :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMoleculeHalogenatedAromaticRings *Clone() const
    {
      return new ExampleDescriptorMoleculeHalogenatedAromaticRings( *this);
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
      // set allowed halogen types
      storage::Vector< chemistry::AtomType> allowed_halogens;
      allowed_halogens.PushBack( chemistry::GetAtomTypes().F_SP2P2P2);
      allowed_halogens.PushBack( chemistry::GetAtomTypes().Cl_SP2P2P2);
      allowed_halogens.PushBack( chemistry::GetAtomTypes().Br_SP2P2P2);
      allowed_halogens.PushBack( chemistry::GetAtomTypes().I_SP2P2P2);

      // Construction
      descriptor::MoleculeHalogenatedAromaticRings total_default; //( descriptor::MoleculeHalogenatedAromaticRings::e_Total);
      descriptor::MoleculeHalogenatedAromaticRings total( descriptor::MoleculeHalogenatedAromaticRings::e_Total, allowed_halogens);
      descriptor::MoleculeHalogenatedAromaticRings max( descriptor::MoleculeHalogenatedAromaticRings::e_Max, allowed_halogens);

      // Data access
      BCL_ExampleCheck( total.GetAlias(), "NAromaticRingHalogensTotal");
      BCL_ExampleCheck( total.GetAlias(), total_default.GetAlias());
      BCL_ExampleCheck( max.GetAlias(), "NAromaticRingHalogensMaxFragment");

      // read in sample molecule
      std::string filename_1( AddExampleInputPathToFilename( e_Chemistry, "aryl_halide.sdf"));
      io::IFStream input_1( filename_1.c_str());

      // make sure the file exists and can be opened
      BCL_ExampleAssert( input_1.is_open(), true);
      BCL_MessageStd( "Reading file " + filename_1);
      chemistry::FragmentEnsemble mol_1( input_1);

      //close stream
      io::File::CloseClearFStream( input_1);

      // check total halogens on aromatic rings
      float total_result( total( mol_1.GetMolecules().FirstElement())( 0));
      BCL_ExampleCheck
      (
        total_result,
        mol_1.GetMolecules().FirstElement().GetMDLPropertyAsVector( "NAromaticRingHalogensTotal")( 0)
      );

      // check fragment with max halogens count
      float max_result( max( mol_1.GetMolecules().FirstElement())( 0));
      BCL_ExampleCheck
      (
        max_result,
        mol_1.GetMolecules().FirstElement().GetMDLPropertyAsVector( "NAromaticRingHalogensMaxFragment")( 0)
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMoleculeHalogenatedAromaticRings

  const ExampleClass::EnumType ExampleDescriptorMoleculeHalogenatedAromaticRings::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMoleculeHalogenatedAromaticRings())
  );

} // namespace bcl
