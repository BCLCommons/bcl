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
#include "descriptor/bcl_descriptor_molecule_total_bond_energy.h"

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
  //! @example example_descriptor_molecule_total_bond_energy.cpp
  //!
  //! @author brownbp1
  //! @date Sep 03, 2019
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMoleculeTotalBondEnergy :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMoleculeTotalBondEnergy *Clone() const
    {
      return new ExampleDescriptorMoleculeTotalBondEnergy( *this);
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

      // default constructor
      descriptor::MoleculeTotalBondEnergy bonde;

      // copy constructor
      descriptor::MoleculeTotalBondEnergy bonde_copy( bonde);

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( bonde.GetAlias(), "MoleculeTotalBondEnergy");

    ///////////////
    // operators //
    ///////////////

      std::string filename( "taxol.sdf");

      // create input stream for reading a small molecule ensemble
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, filename));

      // read in ensemble
      chemistry::FragmentEnsemble ensemble( input);

      // close stream
      io::File::CloseClearFStream( input);

      // make property object
      descriptor::CheminfoProperty bondes_drug( "MoleculeTotalBondEnergy");

      // compute maximum bond energy for taxol according to each statistical potential
      linal::Vector< float> bonde_drug_energies( bondes_drug->SumOverObject( *ensemble.GetMolecules().Begin()));

      BCL_ExampleIndirectCheckWithinAbsTolerance
      (
        bonde_drug_energies( 1), -2.81321, 0.10, "Average bond propensity in taxol"
      );
      BCL_ExampleIndirectCheckWithinAbsTolerance
      (
        bonde_drug_energies( 2), 0.413212, 0.10, "Standard deviation of the bond propensity in taxol"
      );
      BCL_ExampleIndirectCheckWithinAbsTolerance
      (
        bonde_drug_energies( 3), -1.42095, 0.10, "Maximum bond propensity (i.e. score of the least likely bond) in taxol"
      );

    //////////////////////
    // input and output //
    //////////////////////

      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( bonde, bonde_copy),
        true,
        "bonde I/O"
      );

      BCL_ExampleCheck( bondes_drug.GetLabel(), bonde_copy.GetLabel());

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
      } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMoleculeTotalBondEnergy

  const ExampleClass::EnumType ExampleDescriptorMoleculeTotalBondEnergy::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMoleculeTotalBondEnergy())
  );

} // namespace bcl
