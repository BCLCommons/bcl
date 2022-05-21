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
#include "descriptor/bcl_descriptor_atom_polarizability.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "descriptor/bcl_descriptor_atom_misc_property.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "descriptor/bcl_descriptor_molecule_misc_property.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_linear_least_squares.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_atom_polarizability.cpp
  //!
  //! @author kothiwsk
  //! @date Feb 07, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorAtomPolarizability :
    public ExampleInterface
  {
  public:

    ExampleDescriptorAtomPolarizability *Clone() const
    {
      return new ExampleDescriptorAtomPolarizability( *this);
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
      descriptor::AtomPolarizability polarizability;

      // copy constructor
      descriptor::AtomPolarizability polarizability_copy( polarizability);

      // make an atom property that also gets polarizability
      descriptor::CheminfoProperty polarizability_from_atom_properties
      (
        descriptor::GetCheminfoProperties().calc_Polarizability
      );

      // test Clone using an atom property
      BCL_ExampleCheck( descriptor::CheminfoProperty( polarizability).IsDefined(), true);

    /////////////////
    // data access //
    /////////////////

      // get the name of the property without any parameters
      BCL_ExampleCheck( polarizability.GetAlias(), "Atom_Polarizability");

      // get the name of the property with any parameters
      BCL_ExampleCheck( polarizability.GetString(), "Atom_Polarizability");

    ///////////////
    // operators //
    ///////////////

      std::string filename( "small_molecule_descriptors_example_out_high_resolution.sdf");

      // create input stream for reading a smallmolecule ensemble
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, filename));
      // read in ensemble
      chemistry::FragmentEnsemble ensemble( input);
      // close stream
      io::File::CloseClearFStream( input);

      descriptor::AtomMiscProperty chem_property_getter( "bcl_Polarizability", 1);
      descriptor::MoleculeMiscProperty chem_mean_property_getter( "bcl_MeanMolecularPolarizability", 1);

      // keep track of the errors and counts of each atom type
      linal::Vector< float> sum_errors_for_atom_type( chemistry::GetAtomTypes().GetEnumCount(), 0.0);
      linal::Vector< size_t> atom_type_count( chemistry::GetAtomTypes().GetEnumCount(), size_t( 0));

      // keep track of every polarizability calculated with chem and chemistry so we can run
      // a correlation on them in the end
      storage::Vector< float> all_chem_values;
      storage::Vector< float> all_chemistry_values;

      for
      (
        storage::List< chemistry::FragmentComplete>::iterator itr( ensemble.Begin()), itr_end( ensemble.End());
        itr != itr_end;
        ++itr
      )
      {
        // get the (already calculated) polarizability in chem
        linal::Vector< float> chem_values( chem_property_getter.CollectValuesOnEachElementOfObject( *itr));

        // get the vcharge in chemistry
        linal::Vector< float> chemistry_values( polarizability.CollectValuesOnEachElementOfObject( *itr));

        // for each atom in the molecule
        size_t index( 0);
        for
        (
          iterate::Generic< const chemistry::AtomConformationalInterface>
          itr_atoms( itr->GetAtomsIterator());
          itr_atoms.NotAtEnd();
          ++itr_atoms, ++index
        )
        {
          if( util::IsDefined( chem_values( index)) && itr_atoms->GetAtomType()->IsGasteigerAtomType())
          {
            // record info about how far the chem and chemistry values were off
            atom_type_count( itr_atoms->GetAtomType().GetIndex())++;
            sum_errors_for_atom_type( itr_atoms->GetAtomType().GetIndex()) += math::Sqr( chemistry_values( index) - chem_values( index));
            all_chem_values.PushBack( chem_values( index));
            all_chemistry_values.PushBack( chemistry_values( index));
          }
        }
      }

      // write out the rmsd for each atom type
      size_t index( 0);
      for
      (
        chemistry::AtomTypes::const_iterator
          itr( chemistry::GetAtomTypes().Begin()),
          itr_end( chemistry::GetAtomTypes().End());
        itr != itr_end;
        ++itr, ++index
      )
      {
        BCL_MessageDbg
        (
          "Polariz: " + itr->GetName()
          + " Count:" + util::Format()( atom_type_count( index)) + " RMSD:"
          + util::Format()
          (
            atom_type_count( index)
            ? math::Sqrt( sum_errors_for_atom_type( index) / float( atom_type_count( index)))
            : 0.0
          )
        );
      }

      // record the chi-squared between the values from chem and the values from chemistry
      const float chem_to_chemistry_chi_sq
      (
        1.0 - math::LinearLeastSquares::SolutionAndChiSquared
        (
          linal::Matrix< float>( 1, all_chem_values.GetSize(), all_chem_values),
          linal::Vector< float>( all_chemistry_values)
        ).Second()
      );

      // check that chem and chemistry gave similar results
      BCL_ExampleIndirectCheck
      (
        chem_to_chemistry_chi_sq > 0.95,
        true,
        "correlation of chem and chemistry uncorrected polarizability (R^2="
        + util::Format()( chem_to_chemistry_chi_sq) + ")"
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( descriptor::CheminfoProperty( polarizability), descriptor::CheminfoProperty()),
        true,
        "AtomPolarizability I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorAtomPolarizability

  const ExampleClass::EnumType ExampleDescriptorAtomPolarizability::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorAtomPolarizability())
  );

} // namespace bcl
