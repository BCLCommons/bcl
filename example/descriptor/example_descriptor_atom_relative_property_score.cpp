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
#include "descriptor/bcl_descriptor_atom_relative_property_score.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "descriptor/bcl_descriptor_atom_misc_property.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_linear_least_squares.h"
#include "sdf/bcl_sdf_molecule_reading_pref.h"
// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_atom_relative_property_score.cpp
  //!
  //! @author brownbp1
  //! @date Dec 24, 2020
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorAtomRelativePropertyScore :
    public ExampleInterface
  {
  public:

    ExampleDescriptorAtomRelativePropertyScore *Clone() const
    {
      return new ExampleDescriptorAtomRelativePropertyScore( *this);
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

      // cheminfo property of interest
      descriptor::CheminfoProperty xlogp( "XLogP");

      // default constructor
      descriptor::AtomRelativePropertyScore score( xlogp);

      // copy constructor
      descriptor::AtomRelativePropertyScore score_copy( score);

    /////////////////
    // data access //
    /////////////////

      // get the alias
      BCL_ExampleCheck( score.GetAlias(), "Atom_RelativePropertyScore");

    ///////////////
    // operators //
    ///////////////

      // make IFstream
      io::IFStream input, input_base;

      // load input
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "diazepam.sdf"));
      BCL_ExampleMustOpenInputFile( input_base, AddExampleInputPathToFilename( e_Chemistry, "diazepam_base.sdf"));

      // get our molecules
      chemistry::FragmentEnsemble ensemble( input, sdf::e_Saturate);
      chemistry::FragmentEnsemble ensemble_base( input_base, sdf::e_Saturate);
      chemistry::FragmentComplete diazepam( ensemble.GetMolecules().FirstElement());

      // close IF stream
      io::File::CloseClearFStream( input);
      io::File::CloseClearFStream( input_base);

      // Set the reference ensemble
      score.SetReferenceMols( ensemble_base);

      // score
      linal::Vector< float> scores( score.SumOverObject( diazepam));

      // Sum
      BCL_ExampleIndirectCheckWithinAbsTolerance
      (
        scores( 0), 13.3664, 0.1, "Sum of per atom contributions to XLogP relative to the reference"
      );

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorAtomTopologicalPolarSurface

  const ExampleClass::EnumType ExampleDescriptorAtomRelativePropertyScore::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorAtomRelativePropertyScore())
  );

} // namespace bcl
