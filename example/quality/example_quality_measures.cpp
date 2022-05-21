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
#include "quality/bcl_quality_measures.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_quality_measures.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleQualityMeasures :
    public ExampleInterface
  {
  public:

    ExampleQualityMeasures *Clone() const
    {
      return new ExampleQualityMeasures( *this);
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

      // construct qualities from each class type
      BCL_MessageStd
      (
        "Constructing Quality objects from all known Quality measures"
      );

      // construct all measures
      quality::Measure rmsd( quality::GetMeasures().e_RMSD);
      quality::Measure rmsd_no_superimpose( quality::GetMeasures().e_RMSD_NoSuperimposition);
      quality::Measure rmsd_xy_superimpose( quality::GetMeasures().e_RMSD_XYSuperimposition);
      quality::Measure dme( quality::GetMeasures().e_DME);
      quality::Measure dmf_ts( quality::GetMeasures().e_DMF_TS);
      quality::Measure dmf_ha( quality::GetMeasures().e_DMF_HA);
      quality::Measure lcs( quality::GetMeasures().e_LCS);
      quality::Measure gdt_ha( quality::GetMeasures().e_GDT_HA);
      quality::Measure gdt_ts( quality::GetMeasures().e_GDT_TS);
      quality::Measure gdt_1a( quality::GetMeasures().e_GDT_1A);
      quality::Measure gdt_2a( quality::GetMeasures().e_GDT_2A);
      quality::Measure gdt_4a( quality::GetMeasures().e_GDT_4A);
      quality::Measure gdt_8a( quality::GetMeasures().e_GDT_8A);
      quality::Measure maxsub( quality::GetMeasures().e_MaxSub);
      quality::Measure zero( quality::GetMeasures().e_Zero);

      // check that the constructors were correct
      BCL_ExampleCheck( rmsd, quality::GetMeasures().e_RMSD);
      BCL_ExampleCheck( rmsd_no_superimpose, quality::GetMeasures().e_RMSD_NoSuperimposition);
      BCL_ExampleCheck( rmsd_xy_superimpose, quality::GetMeasures().e_RMSD_XYSuperimposition);
      BCL_ExampleCheck( dme, quality::GetMeasures().e_DME);
      BCL_ExampleCheck( dmf_ts, quality::GetMeasures().e_DMF_TS);
      BCL_ExampleCheck( dmf_ha, quality::GetMeasures().e_DMF_HA);
      BCL_ExampleCheck( lcs, quality::GetMeasures().e_LCS);
      BCL_ExampleCheck( gdt_ha, quality::GetMeasures().e_GDT_HA);
      BCL_ExampleCheck( gdt_ts, quality::GetMeasures().e_GDT_TS);
      BCL_ExampleCheck( gdt_1a, quality::GetMeasures().e_GDT_1A);
      BCL_ExampleCheck( gdt_2a, quality::GetMeasures().e_GDT_2A);
      BCL_ExampleCheck( gdt_4a, quality::GetMeasures().e_GDT_4A);
      BCL_ExampleCheck( gdt_8a, quality::GetMeasures().e_GDT_8A);
      BCL_ExampleCheck( maxsub, quality::GetMeasures().e_MaxSub);
      BCL_ExampleCheck( zero, quality::GetMeasures().e_Zero);

      // construct undefined AAClass
      BCL_MessageStd( "Constructing an undefined measure");
      quality::Measure measure_undefined( util::GetUndefined< quality::Measure>());
      BCL_ExampleCheck( measure_undefined.IsDefined(), false);

      // use copy constructor
      BCL_MessageStd( "Calling copy constructor");
      quality::Measure measure_dme_copy( dme);
      BCL_ExampleCheck( measure_dme_copy, dme);

    /////////////////
    // data access //
    /////////////////

      // example check
      BCL_ExampleCheck( quality::GetMeasures().GetClassIdentifier(), "bcl::quality::Measures");

      // display total number of Measures
      BCL_MessageStd
      (
        "The total number of measures is " + util::Format()( quality::GetMeasures().GetEnumCount())
      );

      // check number of cutoffs for TS and HA
      BCL_ExampleCheck( quality::Measures::GetDistanceCutoffsHA().GetSize(), 4);
      BCL_ExampleCheck( quality::Measures::GetDistanceCutoffsTS().GetSize(), 4);

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleQualityMeasures

  const ExampleClass::EnumType ExampleQualityMeasures::s_Instance
  (
    GetExamples().AddEnum( ExampleQualityMeasures())
  );

} // namespace bcl
