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
#include "descriptor/bcl_descriptor_atom_type_property_retriever.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_atom_type_property_retriever.cpp
  //!
  //! @author kothiwsk
  //! @date Feb 11, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorAtomTypePropertyRetriever :
    public ExampleInterface
  {
  public:

    ExampleDescriptorAtomTypePropertyRetriever *Clone() const
    {
      return new ExampleDescriptorAtomTypePropertyRetriever( *this);
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
      descriptor::AtomTypePropertyRetriever unspecified_property;

      BCL_ExampleCheck( unspecified_property.GetAlias(), "UnknownAtomProperty");

      // constructor from property
      descriptor::AtomTypePropertyRetriever pien_getter( chemistry::AtomTypeData::e_PiOrbitalElectronegativityMulliken);

      // copy constructor
      descriptor::AtomTypePropertyRetriever pien_getter_copy( pien_getter);

      // make an atom property that also gets the raw pi-electronegativity
      descriptor::CheminfoProperty pien_getter_from_atom_properties
      (
        descriptor::GetCheminfoProperties().calc_PiOrbitalElectronegativityMulliken
      );

    /////////////////
    // data access //
    /////////////////

      // get the name of the property without any parameters
      BCL_ExampleCheck( pien_getter.GetAlias(), "Atom_PiOrbitalElectronegativityMulliken");

      // get the name of the property with any parameters (there shouldn't be any parameters in this case)
      BCL_ExampleCheck( pien_getter.GetString(), "Atom_PiOrbitalElectronegativityMulliken");

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
      linal::Vector< float> taxol_pien( pien_getter.CollectValuesOnEachElementOfObject( taxol));
      linal::Vector< float> taxol_pien_via_atom_property
      (
        pien_getter_from_atom_properties->CollectValuesOnEachElementOfObject( taxol)
      );
      linal::Vector< float> taxol_pien_via_small_molecule( taxol.GetNumberAtoms(), 0.0);

      // make sure that the vectors have the same size, otherwise the loop below will fail
      BCL_ExampleAssert( taxol_pien.GetSize(), taxol_pien_via_atom_property.GetSize());
      BCL_ExampleAssert( taxol_pien.GetSize(), taxol_pien_via_small_molecule.GetSize());

      linal::Vector< float>::iterator taxol_itr( taxol_pien_via_small_molecule.Begin());
      linal::Vector< float>::iterator taxol_itr_end( taxol_pien_via_small_molecule.End());

      // check that we get the same values when accessing the atoms directly
      for
      (
        iterate::Generic< const chemistry::AtomConformationalInterface>
          itr_atoms( taxol.GetAtomsIterator());
        taxol_itr != taxol_itr_end;
        ++taxol_itr, ++itr_atoms
      )
      {
        *taxol_itr
          = itr_atoms->GetAtomType()->GetAtomTypeProperty
            (
              chemistry::AtomTypeData::e_PiOrbitalElectronegativityMulliken
            );
      }

      // replace nans with 0's, because not all the atom types have defined pi-orbital electronegativities
      for
      (
        linal::Vector< float>::iterator itr_via_small_molecule( taxol_pien_via_small_molecule.Begin()),
          itr_via_small_molecule_end( taxol_pien_via_small_molecule.End()),
          itr_via_atom_property( taxol_pien_via_atom_property.Begin()),
          itr_via_property_retriever( taxol_pien.Begin());
        itr_via_small_molecule != itr_via_small_molecule_end;
        ++itr_via_small_molecule, ++itr_via_atom_property, ++itr_via_property_retriever
      )
      {
        if( !util::IsDefined( *itr_via_small_molecule))
        {
          *itr_via_small_molecule = 0.0;
        }

        if( !util::IsDefined( *itr_via_atom_property))
        {
          *itr_via_atom_property = 0.0;
        }

        if( !util::IsDefined( *itr_via_property_retriever))
        {
          *itr_via_property_retriever = 0.0;
        }
      }

      // check that the vectors are equal
      BCL_ExampleIndirectCheck
      (
        std::equal( taxol_pien.Begin(), taxol_pien.End(), taxol_pien_via_atom_property.Begin()),
        true,
        "descriptor::AtomProperty( AtomTypePropertyRetriever) equivalence to AtomTypePropertyRetriever"
      );

      // check that the vectors are equal
      BCL_ExampleIndirectCheck
      (
        std::equal( taxol_pien.Begin(), taxol_pien.End(), taxol_pien_via_small_molecule.Begin()),
        true,
        "descriptor::AtomProperty( AtomTypePropertyRetriever) equivalence to\n"
        + std::string( "GetAtomType()->GetAtomTypeProperty( chemistry::AtomTypeData::e_VDWaalsRadius)")
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( descriptor::CheminfoProperty( pien_getter), descriptor::CheminfoProperty()),
        true,
        "AtomTypePropertyRetriever I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
      } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorAtomTypePropertyRetriever

  const ExampleClass::EnumType ExampleDescriptorAtomTypePropertyRetriever::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorAtomTypePropertyRetriever())
  );

} // namespace bcl
