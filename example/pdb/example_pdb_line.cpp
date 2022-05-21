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
#include "pdb/bcl_pdb_line.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_types.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_pdb_line.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExamplePdbLine :
    public ExampleInterface
  {
  public:

    ExamplePdbLine *Clone() const
    { return new ExamplePdbLine( *this);}

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
      //array with line like you would find them in pdb-files
      std::string lines[ 4] =
      {
        "SEQRES   1 A  160  MET PRO PRO MET LEU SER GLY LEU LEU ALA ARG LEU VAL",
        "ATOM      1  CA  MET A   1      -0.565  44.678  29.773  1.00  0.00           C",
        "HETATM14063  O   HOH     1      -3.456   2.956  19.529  1.00 47.47           O",
        "HETATM10183  N   MET D  91      45.662  -4.287  27.185  1.00 20.40           N"
      };

      util::ShPtrVector< pdb::Line> pdblines( 5);
      //construct Line from string
      BCL_MessageStd( "construct Lines from strings");
      for( size_t i( 0); i < 4; ++i)
      {
        pdblines( i) = util::ShPtr< pdb::Line>( new pdb::Line( lines[ i]));
      }

      //construct Line from LineType
      BCL_MessageStd( "construct new Line from LineType");
      pdblines( 4) = util::ShPtr< pdb::Line>( new pdb::Line( pdb::GetLineTypes().ATOM));

      //output of all Linetypes
      BCL_MessageStd( "These are the linetypes of the constructed pdb_Lines");
      for
      (
        util::ShPtrVector< pdb::Line>::const_iterator line_itr( pdblines.Begin()), line_itr_end( pdblines.End());
        line_itr != line_itr_end; ++line_itr
      )
      {
        BCL_MessageStd( util::Format()( ( *line_itr)->GetType()));
      }

      //whole pdb line
      BCL_MessageStd( "this is a complete ATOM line:\n" + pdblines( 1)->GetString());
      BCL_MessageStd( "this is the Atom Name: " + pdblines( 1)->GetString( pdb::GetEntryTypes().ATOMName));
      BCL_MessageStd( "this is the ResidueName containing this Atom: " + pdblines( 1)->GetString( pdb::GetEntryTypes().ATOMResidueName));
      BCL_MessageStd( "this is the Chain ID of this Atom: " + util::Format()( pdblines( 1)->GetChar( pdb::GetEntryTypes().ATOMChainID)));
      BCL_MessageStd( "this is the ResidueID of this Atom: " + util::Format()( pdblines( 1)->GetNumericalValue< int>( pdb::GetEntryTypes().ATOMResidueSequenceID)));
      //atom position
      BCL_MessageStd( "these are the x, y and z coordinates of the Atom: " +
                   util::Format()( pdblines( 1)->GetNumericalValue< double>( pdb::GetEntryTypes().ATOMX)) + " " +
                   util::Format()( pdblines( 1)->GetNumericalValue< double>( pdb::GetEntryTypes().ATOMY)) + " " +
                   util::Format()( pdblines( 1)->GetNumericalValue< double>( pdb::GetEntryTypes().ATOMZ)));
      BCL_MessageStd( "to get the Position of the ATOM you can use instead the function Position()\n" +
                   util::Format()( pdblines( 1)->RetrieveCoordinates()));

      //write things to pdbLine
      BCL_MessageStd( "this is an empty ATOM line: " + pdblines( 3)->GetString());
      BCL_MessageStd( "you can Put every Value in the line at the Entries position");
      BCL_MessageStd( "write X coordinate to ATOM line");
      pdblines( 4)->Put( pdb::GetEntryTypes().ATOMX, double( -1.234));
      BCL_MessageStd( pdblines( 4)->GetString());
      BCL_MessageStd( "write Y coordinate to ATOM line");
      pdblines( 4)->Put( pdb::GetEntryTypes().ATOMY, double( 9.876));
      BCL_MessageStd( pdblines( 4)->GetString());
      BCL_MessageStd( "write Atom name to ATOM line");
      pdblines( 4)->Put( pdb::GetEntryTypes().ATOMName, std::string( "CA"));
      BCL_MessageStd( pdblines( 4)->GetString());
      BCL_MessageStd( "write Residue name to ATOM line");
      pdblines( 4)->Put( pdb::GetEntryTypes().ATOMResidueName, biol::GetAATypes().LYS->GetThreeLetterCode());
      BCL_MessageStd( pdblines( 4)->GetString());
      BCL_MessageStd( "write vector3D to ATOM line");
      pdblines( 4)->PutCoordinates( linal::Vector3D( 1.111, 2.222, 3.333));
      BCL_MessageStd( pdblines( 4)->GetString());

      //ask Line if it matches certain criteria
      BCL_MessageStd( "does line contain CA in GetEntryTypes().AtomName entry? " +
                   util::Format()( pdblines( 4)->MatchesCriteria( storage::Pair< pdb::EntryType, std::string>( pdb::GetEntryTypes().ATOMName, "CA"))));

      BCL_MessageStd( "does line contain GLY in GetEntryTypes().AtomResidueName entry? " +
                   util::Format()( pdblines( 4)->MatchesCriteria( storage::Pair< pdb::EntryType, std::string>( pdb::GetEntryTypes().ATOMResidueName, "GLY"))));

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( pdb::GetLineTypes().HET->GetLastEntryType(), pdb::GetEntryTypes().HETDescription);
      BCL_ExampleCheck( pdb::GetLineTypes().HET->GetNumberOfEntryTypes(), 6);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExamplePdbLine

  const ExampleClass::EnumType ExamplePdbLine::s_Instance
  (
    GetExamples().AddEnum( ExamplePdbLine())
  );

} // namespace bcl
