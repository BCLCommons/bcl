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
#include "descriptor/bcl_descriptor_atom_type_number.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_atom_type_number.cpp
  //!
  //! @author mendenjl
  //! @date Jun 25, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorAtomTypeNumber :
    public ExampleInterface
  {
  public:

    ExampleDescriptorAtomTypeNumber *Clone() const
    {
      return new ExampleDescriptorAtomTypeNumber( *this);
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
      descriptor::AtomTypeNumber atomic_number;

    /////////////////
    // data access //
    /////////////////

      // get the name of the property without any parameters
      BCL_ExampleCheck( atomic_number.GetAlias(), "AtomTypeNumber");

      // get the name of the property with any parameters
      BCL_ExampleCheck( atomic_number.GetString(), "AtomTypeNumber");

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
      linal::Vector< float> taxol_atomic_number( atomic_number.CollectValuesOnEachElementOfObject( taxol));

      linal::Vector< float> taxol_atomic_number_via_atom_property
      (
        descriptor::GetCheminfoProperties().calc_AtomTypeNumber->CollectValuesOnEachElementOfObject( taxol)
      );

      linal::Vector< float> taxol_atomic_number_via_small_molecule( taxol.GetNumberAtoms(), 0.0);

      // make sure that the vectors have the same size, otherwise the loop below will fail
      BCL_ExampleAssert( taxol_atomic_number.GetSize(), taxol_atomic_number_via_atom_property.GetSize());
      BCL_ExampleAssert( taxol_atomic_number.GetSize(), taxol_atomic_number_via_small_molecule.GetSize());

      linal::Vector< float>::iterator taxol_itr( taxol_atomic_number_via_small_molecule.Begin());
      linal::Vector< float>::iterator taxol_itr_end( taxol_atomic_number_via_small_molecule.End());

      // check that we get the same values when accessing the atoms directly
      for
      (
        iterate::Generic< const chemistry::AtomConformationalInterface>
          itr_atoms( taxol.GetAtomsIterator());
        taxol_itr != taxol_itr_end;
        ++taxol_itr, ++itr_atoms
      )
      {
        *taxol_itr = itr_atoms->GetAtomType().GetIndex();
      }

      // check that the vectors are equal
      BCL_ExampleIndirectCheck
      (
        std::equal( taxol_atomic_number.Begin(), taxol_atomic_number.End(), taxol_atomic_number_via_atom_property.Begin()),
        true,
        "descriptor::AtomProperty( AtomTypeNumber) equivalence to AtomTypeNumber"
      );

      // check that the vectors are equal
      BCL_ExampleIndirectCheck
      (
        std::equal( taxol_atomic_number.Begin(), taxol_atomic_number.End(), taxol_atomic_number_via_small_molecule.Begin()),
        true,
        "descriptor::AtomProperty( AtomTypeNumber) equivalence to GetElementType()->GetAtomTypeNumber()"
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry
        (
          util::Implementation< descriptor::Base< chemistry::AtomConformationalInterface, float> >( atomic_number),
          util::Implementation< descriptor::Base< chemistry::AtomConformationalInterface, float> >()
        ),
        true,
        "AtomTypeNumber I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorAtomTypeNumber

  const ExampleClass::EnumType ExampleDescriptorAtomTypeNumber::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorAtomTypeNumber())
  );

} // namespace bcl
