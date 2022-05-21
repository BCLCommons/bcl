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
#include "pdb/bcl_pdb_site.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_pdb_site.cpp
  //!
  //! @author woetzen
  //! @date Feb 1, 2012
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExamplePdbSite :
    public ExampleInterface
  {
  public:

    ExamplePdbSite *Clone() const
    { return new ExamplePdbSite( *this);}

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
      const std::string pdb_str( AddExampleInputPathToFilename( e_Biology, "1lgh.pdb"));

      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, pdb_str);
      pdb::Handler pdb_reader( read);
      io::File::CloseClearFStream( read);

    ////////////////
    // operations //
    ////////////////

      // get all sites from the pdb
      const util::ShPtrList< pdb::Site> sites( pdb_reader.GetSites());

      BCL_ExampleAssert( sites.GetSize(), 64);
      util::ShPtrList< pdb::Site>::const_iterator bind_itr( sites.Begin());
      storage::AdvanceIterator( bind_itr, sites.End(), 12);
      const pdb::Site &binder( **bind_itr);
      BCL_MessageStd( "first binding site:\n" + util::Format()( binder));

    //////////////////////
    // input and output //
    //////////////////////

      WriteBCLObject( binder);
      pdb::Site site_read;
      ReadBCLObject( site_read);
      BCL_ExampleCheck( binder.GetName()              , site_read.GetName());
      BCL_ExampleCheck( binder.GetEvidenceCode()      , site_read.GetEvidenceCode());
      BCL_ExampleCheck( binder.GetDescription()       , site_read.GetDescription());
      BCL_ExampleCheck( binder.GetChainResidues()     , site_read.GetChainResidues());
      BCL_ExampleCheck( binder.GetLigand().IsDefined(), site_read.GetLigand().IsDefined());

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExamplePdbSite

  const ExampleClass::EnumType ExamplePdbSite::s_Instance
  (
    GetExamples().AddEnum( ExamplePdbSite())
  );

} // namespace bcl
