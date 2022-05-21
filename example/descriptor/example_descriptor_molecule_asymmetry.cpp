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
#include "descriptor/bcl_descriptor_molecule_asymmetry.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_molecule_ensemble.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_molecule_asymmetry.cpp
  //!
  //! @author kothiwsk
  //! @date Feb 13, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMoleculeAsymmetry :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMoleculeAsymmetry *Clone() const
    {
      return new ExampleDescriptorMoleculeAsymmetry( *this);
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
      // setup input stream
      io::IFStream input_sdf;

      // read in molecule
      const std::string filename_in( AddExampleInputPathToFilename( e_Chemistry, "corina_taxol_out.sdf"));
      BCL_ExampleMustOpenInputFile( input_sdf, filename_in);

      chemistry::MoleculeEnsemble small_mol_ensemble( input_sdf);
      io::File::CloseClearFStream( input_sdf);
      BCL_MessageStd
      (
        "Size of Ensemble: " + util::Format()( small_mol_ensemble.GetSize())
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      descriptor::CheminfoProperty atom_property( "Constant(1)");

      // test default constructor
      descriptor::MoleculeAsymmetry mol_molecular_asymmetry( atom_property, 24, 0.1, 100, false);

      // copy constructor
      descriptor::MoleculeAsymmetry copy_mol_molecular_asymmetry( mol_molecular_asymmetry);

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( mol_molecular_asymmetry.GetAlias(), "MolecularAsymmetry");

    ////////////////
    // operations //
    ////////////////

      //Check setting and getting atom properties
      BCL_ExampleCheck( copy_mol_molecular_asymmetry.GetChemInfoProperty(), atom_property);
      BCL_ExampleCheck( mol_molecular_asymmetry.GetNumberSteps(), 24);
      BCL_ExampleCheckWithinTolerance( mol_molecular_asymmetry.GetStepSize(), 0.1, 0.001);
      BCL_ExampleCheck( mol_molecular_asymmetry.GetTemperature(), 100.0);

      // test the actual code generation
      for
      (
          storage::List< chemistry::MoleculeComplete>::const_iterator
            itr_mols( small_mol_ensemble.GetMolecules().Begin()), itr_mols_end( small_mol_ensemble.GetMolecules().End());
          itr_mols != itr_mols_end;
          ++itr_mols
      )
      {
        linal::Vector< float> stored_molecular_asymmetry( mol_molecular_asymmetry( *itr_mols));
        BCL_ExampleCheckWithinTolerance( stored_molecular_asymmetry.First(), 0.311581, 0.001);
        // test the stream operators
        BCL_MessageStd
        (
          "Outputting mol_molecular_asymmetry: " + util::Format()( stored_molecular_asymmetry)
        );
      }

    //////////////////////
    // input and output //
    //////////////////////

      // Write object file for molecular asymmetry
      WriteBCLObject( mol_molecular_asymmetry);

      return 0;

    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMoleculeAsymmetry

  const ExampleClass::EnumType ExampleDescriptorMoleculeAsymmetry::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMoleculeAsymmetry())
  );

} // namespace bcl
