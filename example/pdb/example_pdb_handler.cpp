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
#include "pdb/bcl_pdb_handler.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_pdb_handler.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExamplePdbHandler :
    public ExampleInterface
  {
  public:

    ExamplePdbHandler *Clone() const
    { return new ExamplePdbHandler( *this);}

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
      const std::string pdb_str( AddExampleInputPathToFilename( e_Biology, "1C1D.pdb"));

      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, pdb_str);
      pdb::Handler pdb_reader( read);
      io::File::CloseClearFStream( read);

      const std::string pdb_filename( AddExampleOutputPathToFilename( pdb::GetNamespaceIdentifier(), "test.pdb"));
      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, pdb_filename);
      pdb_reader.WriteLines( write);
      io::File::CloseClearFStream( write);
      const double x( 12.34);
      BCL_MessageStd( " Conversion of any data type into string of a size non smaller than the data:");
      BCL_MessageStd( ">" + util::Format().R().W( 7)( x) + "<");
      const std::string str( "abcd");
      BCL_MessageStd( ">" + util::Format().R().W( 7)( str) + "<");

      BCL_MessageStd
      (
        " Getting the line type of entry H_ChainID_Initial:  " + pdb::GetEntryTypes().HELIXChainID_Initial->GetLineType().GetName()
      );
      BCL_MessageStd
      (
        " Getting the line type of entry GetEntryTypes().AtomResidueSequenceID:  " + pdb::GetEntryTypes().ATOMResidueSequenceID->GetLineType().GetName()
      );
      BCL_MessageStd
      (
        " Getting the line type of entry GetEntryTypes().AtomCharge:  " + pdb::GetEntryTypes().ATOMCharge->GetLineType().GetName()
      );

      BCL_MessageStd( "Searching lines from a set of criteria. =================\n");
      storage::Vector< storage::Pair< std::string, pdb::EntryType> > crit( 2), criteria( 1);
      storage::Vector< pdb::Line> lines;

      crit( 0).First() = "CA";
      crit( 0).Second() = pdb::GetEntryTypes().ATOMName;
      crit( 1).First() = "ASN";
      crit( 1).Second() = pdb::GetEntryTypes().ATOMResidueName;
      BCL_MessageStd( "The search criteria: ");
      BCL_MessageStd( util::Format()( crit));
      BCL_MessageStd( "Lines in PDB fulfilling the criteria:  ");
      BCL_MessageStd( "Storing the positions in a Vector3D (well, only the first two are displayed):  ");
      util::ShPtrVector< linal::Vector3D> MVEC;
//      MVEC( 0)->Write( util::GetLogger());
//      MVEC( 1)->Write( util::GetLogger());

      criteria( 0).First() = "A";
      criteria( 0).Second() = pdb::GetEntryTypes().HELIXChainID_Initial;
      BCL_MessageStd( "Criteria: " + util::Format()( criteria) + "\nLines:");

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExamplePdbHandler

  const ExampleClass::EnumType ExamplePdbHandler::s_Instance
  (
    GetExamples().AddEnum( ExamplePdbHandler())
  );

} // namespace bcl
