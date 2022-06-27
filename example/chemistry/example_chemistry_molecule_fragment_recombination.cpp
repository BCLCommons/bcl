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
#include "chemistry/bcl_chemistry_molecule_fragment_recombination.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "graph/bcl_graph_const_graph.h"
#include "io/bcl_io_file.h"
#include "random/bcl_random_uniform_distribution.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_molecule_fragment_recombination.cpp
  //!
  //! @author brownbp1
  //! @date Mar 04, 2021
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryMoleculeFragmentRecombination :
    public ExampleInterface
  {
  public:

    ExampleChemistryMoleculeFragmentRecombination *Clone() const
    {
      return new ExampleChemistryMoleculeFragmentRecombination( *this);
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

      // construct from properties
      chemistry::MoleculeFragmentRecombination recombinator
      (
        AddExampleInputPathToFilename( e_Chemistry, "recombinator_b.sdf")
      );

      // clone function
      util::ShPtr< chemistry::MoleculeFragmentRecombination> recombinator_clone( recombinator.Clone());
      BCL_ExampleCheck( recombinator_clone.IsDefined(), true);

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // read in example molecules
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "recombinator_b.sdf"));
      chemistry::FragmentEnsemble mols( input);
      chemistry::FragmentComplete mol_a( mols.GetMolecules().FirstElement());
      io::File::CloseClearFStream( input);

      // Run
//      chemistry::FragmentEnsemble recombined_ensemble( recombinator( mol_a));

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryMoleculeFragmentRecombination

  const ExampleClass::EnumType ExampleChemistryMoleculeFragmentRecombination::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryMoleculeFragmentRecombination())
  );

} // namespace bcl
