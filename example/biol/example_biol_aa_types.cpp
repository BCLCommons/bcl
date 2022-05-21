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
#include "biol/bcl_biol_aa_types.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_aa_types.cpp
  //!
  //! @author karakam
  //! @date
  //! @remarks status empty
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolAATypes :
    public ExampleInterface
  {
  public:

    ExampleBiolAATypes *Clone() const
    { return new ExampleBiolAATypes( *this);}

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
      BCL_MessageStd( "name of the first amino acid = "              + biol::AAType(  1).GetName());
      BCL_MessageStd( "one letter code of the second amino acid = "  + util::Format()( biol::AAType(  2)->GetOneLetterCode()));
      BCL_MessageStd( "three letter code of the third amino acid = " + biol::AAType(  3)->GetThreeLetterCode());
      BCL_MessageStd( "AAType as enum of the 4th amino acid = "      + util::Format()( biol::AAType(  4)));

      BCL_MessageStd( "helix probability of th 6th amino acid = "
                + util::Format()( biol::AAType(  5)->GetAAProperty( biol::AATypeData::e_HelixProbability)));
      BCL_MessageStd( "hydrophobicity of the 7th amino acid = "
                + util::Format()( biol::AAType(  6)->GetAAProperty( biol::AATypeData::e_Hydrophobicity)));
      BCL_MessageStd( "Isoelectric point of the 8th amino acid = "
                + util::Format()( biol::AAType(  7)->GetAAProperty( biol::AATypeData::e_IsoelectricPoint)));
      BCL_MessageStd( "polarizebility of the 9th amino acid = "
                + util::Format()( biol::AAType(  8)->GetAAProperty( biol::AATypeData::e_Polarizability)));
      BCL_MessageStd( "Sterical parameter of the 10th amino acid = "
                + util::Format()( biol::AAType(  9)->GetAAProperty( biol::AATypeData::e_StericalParameter)));
      BCL_MessageStd( "strand probability of the 11th amino acid = "
                + util::Format()( biol::AAType( 10)->GetAAProperty( biol::AATypeData::e_StrandProbability)));
      BCL_MessageStd( "volume of the 12th amino acid = "
                + util::Format()( biol::AAType( 11)->GetAAProperty( biol::AATypeData::e_Volume)));
      BCL_MessageStd( "mass of the 13th amino acid = "
                + util::Format()( biol::AAType( 12)->GetAAProperty( biol::AATypeData::e_Mass)));

      // write enums to file
      WriteBCLObject( biol::GetAATypes());

      // a functional test to ensure that fragments of all AAs can be created as needed
      for
      (
        biol::AATypes::const_iterator itr( biol::GetAATypes().Begin()), itr_end( biol::GetAATypes().End());
        itr != itr_end;
        ++itr
      )
      {
        ( *itr)->GetFragment( true, true);
      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBiolAATypes

  const ExampleClass::EnumType ExampleBiolAATypes::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolAATypes())
  );

} // namespace bcl

