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
#include "descriptor/bcl_descriptor_atom_sigma_charge.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "descriptor/bcl_descriptor_atom_misc_property.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_linear_least_squares.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_atom_sigma_charge.cpp
  //!
  //! @author kothiwsk, mendenjl, geanesar
  //! @date Dec 14, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorAtomSigmaCharge :
    public ExampleInterface
  {
  public:

    ExampleDescriptorAtomSigmaCharge *Clone() const
    {
      return new ExampleDescriptorAtomSigmaCharge( *this);
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
      descriptor::AtomSigmaCharge sigma_charge;

      // copy constructor
      descriptor::AtomSigmaCharge sigma_charge_copy( sigma_charge);

      // default constructor
      descriptor::AtomSigmaCharge sigma_en( false);

      // copy constructor
      descriptor::AtomSigmaCharge sigma_en_copy( sigma_en);

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( sigma_charge.GetAlias(), "Atom_SigmaCharge");

      BCL_ExampleCheck( sigma_charge.GetString(), "Atom_SigmaCharge");

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

      descriptor::AtomMiscProperty chem_value_getter( "bcl_SigmaCharge", 1);

      // keep track of the errors and counts of each atom type
      linal::Vector< float> sum_errors_for_atom_type( chemistry::GetAtomTypes().GetEnumCount(), 0.0);
      linal::Vector< size_t> atom_type_count( chemistry::GetAtomTypes().GetEnumCount(), size_t( 0));

      // keep track of every sigma charge calculated with chem and chemistry so we can run
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
        // get the sigma charge in chemistry
        linal::Vector< float> chemistry_values
        (
          descriptor::GetCheminfoProperties().calc_SigmaCharge->CollectValuesOnEachElementOfObject( *itr)
        );

        // get the (already calculated) sigma charge in chem
        linal::Vector< float> chem_values( chem_value_getter.CollectValuesOnEachElementOfObject( *itr));

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
          "SigChg: " + itr->GetName()
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
        chem_to_chemistry_chi_sq > 0.7,
        true,
        "Sigma charge similarity between chem and chemistry atom types was " + util::Format()( chem_to_chemistry_chi_sq)
      );

      // sigma electronegativity

      // create input stream for reading a smallmolecule ensemble
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "taxol_descriptors.out"));
      // read in ensemble
      chemistry::FragmentEnsemble ensemble_a( input);

      // close stream
      io::File::CloseClearFStream( input);

      // make a miscellaneous property retriever for sigma_en
      descriptor::AtomMiscProperty sigma_en_from_chem( "bcl_SigmaEN", 1);

      // check that the number of sigma_en from the old chem namespace is similar to the values we now get
      for
      (
        storage::List< chemistry::FragmentComplete>::iterator itr( ensemble_a.Begin()), itr_end( ensemble_a.End());
        itr != itr_end;
        ++itr
      )
      {
        BCL_MessageDbg
        (
          "atom types: " + itr->GetAtomTypesString()
        );

        linal::Vector< float> sig_en( sigma_en.CollectValuesOnEachElementOfObject( *itr));
        if( sigma_en_from_chem.CollectValuesOnEachElementOfObject( *itr).GetSize() > 0)
        {
          BCL_MessageDbg
          (
            "chem bcl_SigmaEN: "
            + util::Format()( sigma_en_from_chem.CollectValuesOnEachElementOfObject( *itr))
          );
        }

        BCL_MessageDbg
        (
          "chemistry SigmaEN: "
          + util::Format()( sig_en)
        );
      }

    //////////////////////
    // input and output //
    //////////////////////

      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( sigma_charge, sigma_charge_copy),
        true,
        "SigmaCharge I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorAtomSigmaCharge

  const ExampleClass::EnumType ExampleDescriptorAtomSigmaCharge::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorAtomSigmaCharge())
  );

} // namespace bcl
