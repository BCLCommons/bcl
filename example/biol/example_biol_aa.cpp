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
#include "biol/bcl_biol_aa.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_classes.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_aa.cpp
  //!
  //! @author karakam
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolAA :
    public ExampleInterface
  {
  public:

    ExampleBiolAA *Clone() const
    { return new ExampleBiolAA( *this);}

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
      // construct AA from aatype, generated from string three letter code "LYS", with seqid 1 and pdbid 8
      biol::AA simpleAA1( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().LYS, 1, 8)));

      BCL_MessageStd( "AA created as LYS, with seqid 1 and pdbid 8 - pdbinsertion code should be \' \'");
      BCL_Example_Check
      (
        simpleAA1.GetSeqID()    == 1,
        "seq id is not 1 but: " + util::Format()( simpleAA1.GetSeqID())
      );
      BCL_Example_Check
      (
        simpleAA1.GetPdbID()    == 8,
        "pdb id is not 8 but: " + util::Format()( simpleAA1.GetPdbID())
      );
      BCL_Example_Check
      (
        simpleAA1.GetPdbICode() == ' ',
        "pdb icode is not \' \' but: \'" + util::Format()( simpleAA1.GetPdbICode()) + "\'"
      );

      // the number of atoms for an aa should be 0
      BCL_MessageStd( "number of atoms in aa" + util::Format()( simpleAA1.GetNumberOfAtoms()));
      BCL_Example_Check
      (
        simpleAA1.GetNumberOfAtoms() == 0,
        "number of atoms should be 0, but is: " + util::Format()( simpleAA1.GetNumberOfAtoms())
      );

      // type of amino acid
      BCL_MessageStd( "name of the aminoacid : " + simpleAA1.GetType().GetName());
      BCL_Example_Check
      (
        simpleAA1.GetType().GetName() == biol::GetAATypes().LYS.GetName(),
        "name is not \"LYS\" but: " + simpleAA1.GetType().GetName()
      );

      // aaclass of that amino acid
      BCL_MessageStd( "aaclass type of aa: " + simpleAA1.GetAAClass().GetName());
      BCL_Example_Check
      (
        simpleAA1.GetAAClass() == biol::GetAAClasses().e_AA,
        "aaclass type of this aa sholuld be AA but is: " + simpleAA1.GetAAClass().GetName()
      );

      // types of atoms
      BCL_MessageStd( "there should be no atoms types in that class: " + util::Format()( simpleAA1.GetTypesOfAtoms().GetSize()));
      BCL_Example_Check
      (
        simpleAA1.GetTypesOfAtoms().IsEmpty(),
        "there should be no atom types in that class"
      );

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBiolAA

  const ExampleClass::EnumType ExampleBiolAA::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolAA())
  );

} // namespace bcl
