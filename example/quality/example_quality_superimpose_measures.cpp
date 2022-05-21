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
#include "quality/bcl_quality_superimpose_measures.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_quality_superimpose_measures.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleQualitySuperimposeMeasures :
    public ExampleInterface
  {
  public:

    ExampleQualitySuperimposeMeasures *Clone() const
    {
      return new ExampleQualitySuperimposeMeasures( *this);
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

      // construct AAClass from each class type
      BCL_MessageStd
      (
        "Constructing objects from all known Quality superimpose measures"
      );

      // construct all measures
      quality::SuperimposeMeasure rmsd( quality::GetSuperimposeMeasures().e_RMSD);
      quality::SuperimposeMeasure rmsd_xy( quality::GetSuperimposeMeasures().e_RMSD_XYSuperimposition);
      quality::SuperimposeMeasure lcs( quality::GetSuperimposeMeasures().e_LCS);
      quality::SuperimposeMeasure gdt_1a( quality::GetSuperimposeMeasures().e_GDT_1A);
      quality::SuperimposeMeasure gdt_2a( quality::GetSuperimposeMeasures().e_GDT_2A);
      quality::SuperimposeMeasure gdt_4a( quality::GetSuperimposeMeasures().e_GDT_4A);
      quality::SuperimposeMeasure gdt_8a( quality::GetSuperimposeMeasures().e_GDT_8A);
      quality::SuperimposeMeasure maxsub( quality::GetSuperimposeMeasures().e_MaxSub);
      quality::SuperimposeMeasure no_superimpose( quality::GetSuperimposeMeasures().e_NoSuperimpose);

      // check that the constructors were correct
      BCL_ExampleCheck( rmsd, quality::GetSuperimposeMeasures().e_RMSD);
      BCL_ExampleCheck( rmsd_xy, quality::GetSuperimposeMeasures().e_RMSD_XYSuperimposition);
      BCL_ExampleCheck( lcs, quality::GetSuperimposeMeasures().e_LCS);
      BCL_ExampleCheck( gdt_1a, quality::GetSuperimposeMeasures().e_GDT_1A);
      BCL_ExampleCheck( gdt_2a, quality::GetSuperimposeMeasures().e_GDT_2A);
      BCL_ExampleCheck( gdt_4a, quality::GetSuperimposeMeasures().e_GDT_4A);
      BCL_ExampleCheck( gdt_8a, quality::GetSuperimposeMeasures().e_GDT_8A);
      BCL_ExampleCheck( maxsub, quality::GetSuperimposeMeasures().e_MaxSub);
      BCL_ExampleCheck( no_superimpose, quality::GetSuperimposeMeasures().e_NoSuperimpose);

      // construct undefined AAClass
      BCL_MessageStd( "Constructing an undefined measure");
      quality::SuperimposeMeasure measure_undefined( util::GetUndefined< quality::SuperimposeMeasure>());
      BCL_ExampleCheck( measure_undefined.IsDefined(), false);

      // use copy constructor
      BCL_MessageStd( "Calling copy constructor");
      quality::SuperimposeMeasure rmsd_copy( rmsd);
      BCL_ExampleCheck( rmsd_copy, rmsd);

    /////////////////
    // data access //
    /////////////////

      // example check
      BCL_ExampleCheck( quality::GetSuperimposeMeasures().GetClassIdentifier(), "bcl::quality::SuperimposeMeasures");

      // display total number of Measures
      BCL_MessageStd
      (
        "The total number of superimpose measures is " + util::Format()( quality::GetSuperimposeMeasures().GetEnumCount())
      );

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

  }; //end ExampleQualitySuperimposeMeasures

  const ExampleClass::EnumType ExampleQualitySuperimposeMeasures::s_Instance
  (
    GetExamples().AddEnum( ExampleQualitySuperimposeMeasures())
  );

} // namespace bcl
