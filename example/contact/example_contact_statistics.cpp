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
#include "contact/bcl_contact_statistics.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container_generator_protein_model.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_contact_statistics.cpp
  //!
  //! @author karakam
  //! @date Sep 23, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleContactStatistics :
    public ExampleInterface
  {
  public:

    ExampleContactStatistics *Clone() const
    {
      return new ExampleContactStatistics( *this);
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
      // read in 1UBI model
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      // read model
      assemble::ProteinModel model( Proteins::GetModel( pdb_filename));

      // initialize ranges
      math::Range< size_t> range_default( contact::GetDefaultSequenceSeparationRange());
      math::Range< size_t> range_short( contact::GetDefaultSequenceSeparationShortRange());
      math::Range< size_t> range_mid( contact::GetDefaultSequenceSeparationMidRange());
      math::Range< size_t> range_long( contact::GetDefaultSequenceSeparationLongRange());

      // neighbor generator

      // initialize a static instance
      const util::ShPtr< math::FunctionInterfaceSerializable< assemble::ProteinModel, assemble::AANeighborListContainer> >
      neighbor_generator
      (
        assemble::AANeighborListContainerGeneratorProteinModel::AANeighborListGenerator
        (
          contact::g_ContactCbDistanceCutoff,
          contact::g_ContactMinSequenceSeparation,
          true,
          false
        )
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      contact::Statistics stats_default;

      // constructor from all enum types
      contact::Statistics stats_enum_counts( contact::Statistics::e_NumberContacts, false);
      BCL_ExampleCheck( stats_enum_counts.GetUseRatio(), false);
      BCL_ExampleCheck( stats_enum_counts.GetSequenceSeparationRange(), range_default);

      contact::Statistics stats_enum_counts_short( contact::Statistics::e_NumberContactsShort, false);
      BCL_ExampleCheck( stats_enum_counts_short.GetUseRatio(), false);
      BCL_ExampleCheck( stats_enum_counts_short.GetSequenceSeparationRange(), range_short);

      contact::Statistics stats_enum_counts_mid( contact::Statistics::e_NumberContactsMid, false);
      BCL_ExampleCheck( stats_enum_counts_mid.GetUseRatio(), false);
      BCL_ExampleCheck( stats_enum_counts_mid.GetSequenceSeparationRange(), range_mid);

      contact::Statistics stats_enum_counts_long( contact::Statistics::e_NumberContactsLong, false);
      BCL_ExampleCheck( stats_enum_counts_long.GetUseRatio(), false);
      BCL_ExampleCheck( stats_enum_counts_long.GetSequenceSeparationRange(), range_long);

      contact::Statistics stats_enum_ratio_short( contact::Statistics::e_RatioContactsShort, false);
      BCL_ExampleCheck( stats_enum_ratio_short.GetUseRatio(), true);
      BCL_ExampleCheck( stats_enum_ratio_short.GetSequenceSeparationRange(), range_short);

      contact::Statistics stats_enum_ratio_mid( contact::Statistics::e_RatioContactsMid, false);
      BCL_ExampleCheck( stats_enum_ratio_mid.GetUseRatio(), true);
      BCL_ExampleCheck( stats_enum_ratio_mid.GetSequenceSeparationRange(), range_mid);

      contact::Statistics stats_enum_ratio_long( contact::Statistics::e_RatioContactsLong, false);
      BCL_ExampleCheck( stats_enum_ratio_long.GetUseRatio(), true);
      BCL_ExampleCheck( stats_enum_ratio_long.GetSequenceSeparationRange(), range_long);

      // constructor from range and ratio and neighbor generator
      contact::Statistics stats_counts( range_default, false, false);
      BCL_ExampleCheck( stats_counts.GetUseRatio(), false);
      BCL_ExampleCheck( stats_counts.GetSequenceSeparationRange(), range_default);

      contact::Statistics stats_counts_long( range_long, false, false);
      BCL_ExampleCheck( stats_counts_long.GetUseRatio(), false);
      BCL_ExampleCheck( stats_counts_long.GetSequenceSeparationRange(), range_long);

      contact::Statistics stats_ratio_long( range_long, true, false);
      BCL_ExampleCheck( stats_ratio_long.GetUseRatio(), true);
      BCL_ExampleCheck( stats_ratio_long.GetSequenceSeparationRange(), range_long);

    /////////////////
    // data access //
    /////////////////

      // test GetStatisticTypeDescriptor
      BCL_ExampleCheck
      (
        contact::Statistics::GetStatisticTypeDescriptor( contact::Statistics::e_NumberContacts), "NumberContacts"
      );
      BCL_ExampleCheck
      (
        contact::Statistics::GetStatisticTypeDescriptor( contact::Statistics::e_NumberContactsShort), "NumberContactsShort"
      );
      BCL_ExampleCheck
      (
        contact::Statistics::GetStatisticTypeDescriptor( contact::Statistics::e_NumberContactsMid), "NumberContactsMid"
      );
      BCL_ExampleCheck
      (
        contact::Statistics::GetStatisticTypeDescriptor( contact::Statistics::e_NumberContactsLong), "NumberContactsLong"
      );
      BCL_ExampleCheck
      (
        contact::Statistics::GetStatisticTypeDescriptor( contact::Statistics::e_RatioContactsShort), "RatioContactsShort"
      );
      BCL_ExampleCheck
      (
        contact::Statistics::GetStatisticTypeDescriptor( contact::Statistics::e_RatioContactsMid), "RatioContactsMid"
      );
      BCL_ExampleCheck
      (
        contact::Statistics::GetStatisticTypeDescriptor( contact::Statistics::e_RatioContactsLong), "RatioContactsLong"
      );

    ///////////////
    // operators //
    ///////////////

      // test operator()
      BCL_MessageStd( "Testing operator()");

      // expected values
      const double expected_counts_total( 134);
      const double expected_counts_short( 25);
      const double expected_counts_mid( 29);
      const double expected_counts_long( 80);
      const double expected_ratio_short( 100.0 * expected_counts_short/ expected_counts_total);
      const double expected_ratio_mid(   100.0 * expected_counts_mid  / expected_counts_total);
      const double expected_ratio_long(  100.0 * expected_counts_long / expected_counts_total);

      // check expected values with the calculated ones
      const double counts_total( stats_counts( model));
      BCL_MessageStd( "counts_total: " + util::Format()( counts_total));
      BCL_ExampleCheck( counts_total, expected_counts_total);

      const double counts_total_b( stats_enum_counts( model));
      BCL_MessageStd( "counts_total_b: " + util::Format()( counts_total_b));
      BCL_ExampleCheck( counts_total_b, expected_counts_total);

      const double counts_short( stats_enum_counts_short( model));
      BCL_MessageStd( "counts_short: " + util::Format()( counts_short));
      BCL_ExampleCheck( counts_short, expected_counts_short);

      const double counts_mid( stats_enum_counts_mid( model));
      BCL_MessageStd( "counts_mid: " + util::Format()( counts_mid));
      BCL_ExampleCheck( counts_mid, expected_counts_mid);

      const double counts_long( stats_enum_counts_long( model));
      BCL_MessageStd( "counts_long: " + util::Format()( counts_long));
      BCL_ExampleCheck( counts_long, expected_counts_long);

      const double ratio_short( stats_enum_ratio_short( model));
      BCL_MessageStd( "ratio_short: " + util::Format()( ratio_short));
      BCL_ExampleCheck( ratio_short, expected_ratio_short);

      const double ratio_mid( stats_enum_ratio_mid( model));
      BCL_MessageStd( "ratio_mid: " + util::Format()( ratio_mid));
      BCL_ExampleCheck( ratio_mid, expected_ratio_mid);

      const double ratio_long( stats_enum_ratio_long( model));
      BCL_MessageStd( "ratio_long: " + util::Format()( ratio_long));
      BCL_ExampleCheck( ratio_long, expected_ratio_long);

    //////////////////////
    // input and output //
    //////////////////////

      // write object
      WriteBCLObject( stats_counts);
      // read object
      contact::Statistics stats_read;
      ReadBCLObject( stats_read);

      // compare values
      BCL_ExampleCheck( stats_counts.GetUseRatio(), stats_read.GetUseRatio());
      BCL_ExampleCheck( stats_counts.GetSequenceSeparationRange(), stats_read.GetSequenceSeparationRange());

      // calculate the counts
      const size_t counts_total_read( stats_read( model));
      BCL_MessageStd( "counts read: " + util::Format()( counts_total_read));
      BCL_ExampleCheck( counts_total_read, expected_counts_total);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleContactStatistics

  const ExampleClass::EnumType ExampleContactStatistics::s_Instance
  (
    GetExamples().AddEnum( ExampleContactStatistics())
  );

} // namespace bcl
