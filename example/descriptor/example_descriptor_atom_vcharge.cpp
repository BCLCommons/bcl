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
#include "descriptor/bcl_descriptor_atom_vcharge.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "descriptor/bcl_descriptor_atom_misc_property.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_linear_least_squares.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_atom_vcharge.cpp
  //!
  //! @author kothiwsk
  //! @date Feb 07, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorAtomVcharge :
    public ExampleInterface
  {
  public:

    ExampleDescriptorAtomVcharge *Clone() const
    {
      return new ExampleDescriptorAtomVcharge( *this);
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
      descriptor::AtomVcharge vcharge;

      // copy constructor
      descriptor::AtomVcharge vcharge_copy( vcharge);

    /////////////////
    // data access //
    /////////////////

      // get the name of the property without any parameters
      BCL_ExampleCheck( vcharge.GetAlias(), "Atom_Vcharge");

      // get the name of the property with any parameters
      BCL_ExampleCheck( vcharge.GetString(), "Atom_Vcharge");

    ///////////////
    // operators //
    ///////////////

      std::string filename( "histidine.sdf");

      // create input stream for reading a smallmolecule ensemble
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, filename));
      // read in ensemble
      chemistry::FragmentEnsemble ensemble( input, sdf::e_Saturate);

      // close stream
      io::File::CloseClearFStream( input);

      // keep track of every vcharge calculated with paper and chemistry so we can run
      // a correlation on them in the end
      storage::Vector< float> all_vcharge_values;
      storage::Vector< float> all_chemistry_values;

      // check that the property we can get from the atoms of the small molecule is the same as we get from the property
      linal::Vector< float> histidine_vcharge( vcharge.CollectValuesOnEachElementOfObject( ensemble.GetMolecules().FirstElement()));

      descriptor::AtomMiscProperty vcharge_property_getter( "Vcharge2003_FromPaper", 1);
      for
      (
        storage::List< chemistry::FragmentComplete>::iterator
          itr( ensemble.Begin()), itr_end( ensemble.End());
        itr != itr_end;
        ++itr
      )
      {
        // get the vcharge in chemistry
        linal::Vector< float> chemistry_values( vcharge.CollectValuesOnEachElementOfObject( ensemble.GetMolecules().FirstElement()));

        // get the (already calculated) vcharge from the paper
        linal::Vector< float> vcharge_values( vcharge_property_getter.CollectValuesOnEachElementOfObject( *itr));

        all_vcharge_values.Append( storage::Vector< float>( vcharge_values.GetSize(), vcharge_values.Begin()));
        all_chemistry_values.Append( storage::Vector< float>( chemistry_values.GetSize(), chemistry_values.Begin()));
      }

      // record the chi-squared between the values from paper and the values from chemistry
      const float vcharge_to_chemistry_chi_sq
      (
        1.0 - math::LinearLeastSquares::SolutionAndChiSquared
        (
          linal::Matrix< float>( size_t( 1), all_vcharge_values.GetSize(), all_chemistry_values),
          linal::Vector< float>( all_chemistry_values)
        ).Second()
      );

      // check that paper and chemistry gave similar results
      if
      (
        !BCL_ExampleIndirectCheck
        (
          vcharge_to_chemistry_chi_sq > 0.95,
          true,
          "Vcharge between paper and chemistry atom types was " + util::Format()( vcharge_to_chemistry_chi_sq)
        )
      )
      {
        // print out the actual vectors too:
        BCL_MessageDbg
        (
          "Chemistry Vcharge: " + util::Format()( all_chemistry_values)
          + "paper vcharge: " + util::Format()( all_vcharge_values)
        );
      }

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( descriptor::CheminfoProperty( vcharge), descriptor::CheminfoProperty()),
        true,
        "AtomVcharge I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorAtomVcharge

  const ExampleClass::EnumType ExampleDescriptorAtomVcharge::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorAtomVcharge())
  );

} // namespace bcl
