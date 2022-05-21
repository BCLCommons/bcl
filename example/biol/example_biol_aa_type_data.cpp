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
#include "biol/bcl_biol_aa_type_data.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_types.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_gnuplot_heatmap.h"
#include "restraint/bcl_restraint_sas_data_parameters.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_aa_type_data.cpp
  //!
  //! @author karakam
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by putnamdk on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolAATypeData :
    public ExampleInterface
  {
  public:

    ExampleBiolAATypeData *Clone() const
    {
      return new ExampleBiolAATypeData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
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
      biol::AATypeData aatypedata_default;

      // construct from just name information
      biol::AATypeData aatypedata_simple( "GAP", '-', false, "GAP");

      // construct from name information and atomtypes, first side chain atom type, and properties
      biol::AATypeData aatypedata_constr
      (
        "ASN",
        'N',
        true,
        "ASN",
        storage::Set< biol::AtomType>::Create
        (
          biol::GetAtomTypes().N, biol::GetAtomTypes().CA, biol::GetAtomTypes().C, biol::GetAtomTypes().O,
          biol::GetAtomTypes().CB, biol::GetAtomTypes().CG, biol::GetAtomTypes().OD1, biol::GetAtomTypes().ND2,
          biol::GetAtomTypes().H, biol::GetAtomTypes().HA, biol::GetAtomTypes().HB2, biol::GetAtomTypes().HB3,
          biol::GetAtomTypes().HD21, biol::GetAtomTypes().HD22
        ), biol::GetAtomTypes().CB,
        0.043, 1.600, 0.130, 2.950, -0.600, 6.520,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,
        0.00, 0.00, 0.00, 7.22, 3.64, 7.22, 3.64, 0.210, 0.220, 0.157, 0.081, -0.182, 0.850, 4.800, -3.500, -0.780,
        -0.200, 0.480, 0.500, 0.180, 0.220, 0.397, -0.143, -0.114, 0.759, 0.055, -0.048, 0.048, -0.220, 0.241, 0.113,
        -0.243, -0.175, -2.750, -2.167, 1.454, 0.00, 0.00, 259.845, 4.35, 12.577,  63.64, 150.80,    1,   1, false, 2.0
      );

      // copy constructor
      biol::AATypeData aatypedata_copy( aatypedata_constr);

      // clone
      util::ShPtr< util::ObjectInterface> ptr( aatypedata_constr.Clone());

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

      // Create glycine object
      const biol::AATypeData &glycine( *biol::GetAATypes().GLY);

      {

        // The format for the vector of ND4< double> is Qvalue, Sasa, ExcludedVolume and Hydration shell
        const double calculated_value( glycine.GetStructureFactor()->operator()
          (
            restraint::SasDataParameters( false, 0.0, 0.0, 1.0, 0.0, 0.0)
          )
        );
        const double expected_value( 9.96896);

        BCL_ExampleIndirectCheck
        (
          math::EqualWithinTolerance( calculated_value, expected_value),
          true,
          "incorrect result for glycine with scattering angle of 0 " + util::Format().FFP( 5)( calculated_value)
        );
      }

      {
        // The format for the vector of ND4< double> is Qvalue, Sasa, ExcludedVolume and Hydration shell
        const double calculated_value
        (
          glycine.GetStructureFactor()->operator()
          (
            restraint::SasDataParameters( false, 1.0, 0.0, 1.0, 0.0, 0.0)
          )
        );

        const double expected_value( 15.79023);

        BCL_ExampleIndirectCheck
        (
          math::EqualWithinTolerance( calculated_value, expected_value),
          true,
          " incorrect result for glycine with scattering angle of 1.0 " + util::Format().FFP( 5)( calculated_value)
        );
      }

      {
        // The format for the vector of ND4< double> is Qvalue, Sasa, ExcludedVolume and Hydration shell
        const double calculated_value
        (
          glycine.GetStructureFactor()->operator()
          (
            restraint::SasDataParameters( false, 3.0, 0.0, 1.0, 0.0, 0.0)
          )
        );

        const double expected_value( 16.14269);

        BCL_ExampleIndirectCheck
        (
          math::EqualWithinTolerance( calculated_value, expected_value),
          true,
          " incorrect result for glycine with scattering angle of 3.0 " + util::Format().FFP( 5)( calculated_value)
        );
      }

      {
        // The format for the vector of ND4< double> is Qvalue, Sasa, ExcludedVolume and Hydration shell
        const double calculated_value
        (
          glycine.GetStructureFactor()->operator()
          (
            restraint::SasDataParameters( false, 4.0, 0.0, 1.0, 0.0, 0.0)
          )
        );

        const double expected_value( 12.47393);

        BCL_ExampleIndirectCheck
        (
          math::EqualWithinTolerance( calculated_value, expected_value),
          true,
          " incorrect result for glycine with scattering angle of 4.0 " + util::Format().FFP( 5)( calculated_value)
        );
      }

    //////////////////////
    // input and output //
    //////////////////////

      // write object
      WriteBCLObject( aatypedata_constr);

      // read object
      biol::AATypeData aatypedata_read;
      ReadBCLObject( aatypedata_read);

      // compared
      BCL_ExampleIndirectCheck
      (
        aatypedata_constr.GetAllowedAtomTypes() == aatypedata_read.GetAllowedAtomTypes(),
        true,
        "write and read aa type data"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBiolAATypeData

  const ExampleClass::EnumType ExampleBiolAATypeData::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolAATypeData())
  );

} // namespace bcl
