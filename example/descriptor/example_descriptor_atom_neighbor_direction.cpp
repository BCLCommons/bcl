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
#include "descriptor/bcl_descriptor_atom_neighbor_direction.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_atom_neighbor_direction.cpp
  //!
  //! @author brownbp1, mendenjl
  //! @date Sep 26, 2017
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorAtomNeighborDirection :
    public ExampleInterface
  {
  public:

    ExampleDescriptorAtomNeighborDirection *Clone() const
    {
      return new ExampleDescriptorAtomNeighborDirection( *this);
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
      descriptor::AtomNeighborDirection atom_neighbor_direction;

      // test Clone using an atom property
      BCL_ExampleCheck( descriptor::CheminfoProperty( atom_neighbor_direction).IsDefined(), true);

    /////////////////
    // data access //
    /////////////////

      // get the name of the property without any parameters
      BCL_ExampleCheck( atom_neighbor_direction.GetAlias(), "Atom_NeighborDirection");

      // get the name of the property with any parameters
      BCL_ExampleCheck( atom_neighbor_direction.GetString(), "Atom_NeighborDirection");

    ///////////////
    // operators //
    ///////////////

      // create input stream for reading a smallmolecule ensemble
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "taxol.sdf"));
      // read in ensemble
      chemistry::FragmentEnsemble taxol_ensemble( input);

      chemistry::FragmentComplete taxol( *taxol_ensemble.Begin());

      // close stream
      io::File::CloseClearFStream( input);

      // check that the property we can get from the atoms of the small molecule is the same as we get from the property
      linal::Vector< float> taxol_atom_neighbor_direction( atom_neighbor_direction.CollectValuesOnEachElementOfObject( taxol));

      linal::Vector< float> atom_neighbor_direction_via_atom_property
      (
        descriptor::GetCheminfoProperties().calc_NeighborDirection->CollectValuesOnEachElementOfObject( taxol)
      );

      linal::Vector< float> atom_neighbor_direction_via_small_molecule( taxol.GetNumberAtoms() * 3, 0.0);

      // make sure that the vectors have the same size, otherwise the loop below will fail
      BCL_ExampleAssert( taxol_atom_neighbor_direction.GetSize(), atom_neighbor_direction_via_atom_property.GetSize());
      BCL_ExampleAssert( taxol_atom_neighbor_direction.GetSize(), atom_neighbor_direction_via_small_molecule.GetSize());
      BCL_ExampleCheckWithinAbsTolerance(taxol_atom_neighbor_direction(0),0.303338,0.00005);
      BCL_ExampleCheckWithinAbsTolerance(taxol_atom_neighbor_direction(6),-0.341041,0.00005);
      BCL_ExampleCheckWithinAbsTolerance(taxol_atom_neighbor_direction(22),0.642075,0.00005);

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( descriptor::CheminfoProperty( atom_neighbor_direction), descriptor::CheminfoProperty()),
        true,
        "AtomNeighborDirection I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorAtomNeighborDirection

  const ExampleClass::EnumType ExampleDescriptorAtomNeighborDirection::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorAtomNeighborDirection())
  );

} // namespace bcl
