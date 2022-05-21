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
#include "descriptor/bcl_descriptor_element_type_property_retriever.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_element_type_property_retriever.cpp
  //!
  //! @author kothiwsk
  //! @date Feb 11, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorElementTypePropertyRetriever :
    public ExampleInterface
  {
  public:

    ExampleDescriptorElementTypePropertyRetriever *Clone() const
    {
      return new ExampleDescriptorElementTypePropertyRetriever( *this);
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
      descriptor::ElementTypePropertyRetriever unspecified_property;

      BCL_ExampleCheck( unspecified_property.GetAlias(), "UnknownElementProperty");

      // constructor from property
      descriptor::ElementTypePropertyRetriever vdw_getter( chemistry::ElementTypeData::e_VDWaalsRadius);

      // copy constructor
      descriptor::ElementTypePropertyRetriever vdw_getter_copy( vdw_getter);

      // make an atom property that also gets the VDW radius
      descriptor::CheminfoProperty vdw_getter_from_atom_properties( descriptor::GetCheminfoProperties().calc_VDWaalsRadius);

    /////////////////
    // data access //
    /////////////////

      // get the name of the property without any parameters
      BCL_ExampleCheck( vdw_getter.GetAlias(), "Atom_VDWaalsRadius");

      // get the name of the property with any parameters (there shouldn't be any parameters in this case)
      BCL_ExampleCheck( vdw_getter.GetString(), "Atom_VDWaalsRadius");

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
      linal::Vector< float> taxol_vdw_radii( vdw_getter.CollectValuesOnEachElementOfObject( taxol));
      linal::Vector< float> taxol_vdw_radii_via_atom_property( vdw_getter_from_atom_properties->CollectValuesOnEachElementOfObject( taxol));
      linal::Vector< float> taxol_vdw_radii_via_small_molecule( taxol.GetNumberAtoms(), 0.0);

      // make sure that the vectors have the same size, otherwise the loop below will fail
      BCL_ExampleAssert( taxol_vdw_radii.GetSize(), taxol_vdw_radii_via_atom_property.GetSize());
      BCL_ExampleAssert( taxol_vdw_radii.GetSize(), taxol_vdw_radii_via_small_molecule.GetSize());

      linal::Vector< float>::iterator taxol_itr( taxol_vdw_radii_via_small_molecule.Begin());
      linal::Vector< float>::iterator taxol_itr_end( taxol_vdw_radii_via_small_molecule.End());

      // check that we get the same values when accessing the atoms directly
      for
      (
        iterate::Generic< const chemistry::AtomConformationalInterface> itr_atoms( taxol.GetAtomsIterator());
        taxol_itr != taxol_itr_end;
        ++taxol_itr, ++itr_atoms
      )
      {
        *taxol_itr = itr_atoms->GetElementType()->GetProperty( chemistry::ElementTypeData::e_VDWaalsRadius);
      }

      // check that the vectors are equal
      BCL_ExampleIndirectCheck
      (
        std::equal( taxol_vdw_radii.Begin(), taxol_vdw_radii.End(), taxol_vdw_radii_via_atom_property.Begin()),
        true,
        "descriptor::AtomProperty( ElementTypePropertyRetriever) equivalence to ElementTypePropertyRetriever"
      );

      // check that the vectors are equal
      BCL_ExampleIndirectCheck
      (
        std::equal( taxol_vdw_radii.Begin(), taxol_vdw_radii.End(), taxol_vdw_radii_via_small_molecule.Begin()),
        true,
        "descriptor::AtomProperty( ElementTypePropertyRetriever) equivalence to\n"
        + std::string( "GetElementType()->GetProperty( chemistry::ElementTypeData::e_VDWaalsRadius)")
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( descriptor::CheminfoProperty( vdw_getter), descriptor::CheminfoProperty()),
        true,
        "ElementTypePropertyRetriever I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
      } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorElementTypePropertyRetriever

  const ExampleClass::EnumType ExampleDescriptorElementTypePropertyRetriever::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorElementTypePropertyRetriever())
  );

} // namespace bcl
