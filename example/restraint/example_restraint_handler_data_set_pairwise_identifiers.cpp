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
#include "restraint/bcl_restraint_handler_data_set_pairwise_identifiers.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "restraint/bcl_restraint_locator_coordinates_first_side_chain_atom.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_handler_data_set_pairwise_identifiers.cpp
  //! @brief TODO: add an detailed description for this example
  //!
  //! @author alexanns
  //! @date Oct 31, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintHandlerDataSetPairwiseIdentifiers :
    public ExampleInterface
  {
  public:

    ExampleRestraintHandlerDataSetPairwiseIdentifiers *Clone() const
    {
      return new ExampleRestraintHandlerDataSetPairwiseIdentifiers( *this);
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
      // make data set
      restraint::DataSetPairwise data_set;

      data_set.Insert
      (
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 32)),
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 48))
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      restraint::HandlerDataSetPairwiseIdentifiers def_constr;
      BCL_ExampleCheck( util::IsDefined( def_constr.GetScore()), false);

      // constructor taking score
      restraint::HandlerDataSetPairwiseIdentifiers param_constr( 3.0);
      BCL_ExampleCheck( param_constr.GetScore(), 3.0);

      // clone constructor
      util::ShPtr< restraint::HandlerDataSetPairwiseIdentifiers> clone_constr( param_constr.Clone());
      BCL_ExampleCheck( clone_constr->GetScore(), 3.0);

    /////////////////
    // data access //
    /////////////////

      // GetScore
      BCL_ExampleCheck( clone_constr->GetScore(), 3.0);

    ////////////////
    // operations //
    ////////////////

      // WriteDataSetPairwise
      {
        io::OFStream write;
        std::string data_set_filename( AddExampleOutputPathToFilename( def_constr, "data_set_a"));
        BCL_ExampleMustOpenOutputFile( write, data_set_filename);
        param_constr.WriteDataSetPairwise( write, data_set);
        io::File::CloseClearFStream( write);
        BCL_ExampleCheck( io::File::FilesMatch( data_set_filename, data_set_filename + ".correct"), true);
      }

      // ReadDataSetPairwise
      {
        io::IFStream read;
        std::string data_set_filename( AddExampleOutputPathToFilename( def_constr, "data_set_a"));
        BCL_ExampleMustOpenInputFile( read, data_set_filename);
        def_constr.ReadDataSetPairwise( read);
        BCL_ExampleCheck( def_constr.GetScore(), 3.0);
        BCL_ExampleCheck( def_constr.GetDataSetPairwise().GetSize(), 1);
      }

      // GetDataSetPairwise
      {
        BCL_ExampleCheck( def_constr.GetScore(), 3.0);
        BCL_ExampleCheck( def_constr.GetDataSetPairwise().GetSize(), 1);
      }

      BCL_MessageDbg( "try writing and reading dataset with more than one data pair");
      // try writing and reading dataset with more than one data pair
      BCL_ExampleAssert
      (
        data_set.Insert
        (
          util::ShPtr< assemble::LocatorAtomCoordinatesInterface>( new restraint::LocatorCoordinatesFirstSideChainAtom( 'B', 2)),
          util::ShPtr< assemble::LocatorAtomCoordinatesInterface>( new restraint::LocatorCoordinatesFirstSideChainAtom( 'C', 16))
        ).second, true
      );

      // WriteDataSetPairwise
      BCL_MessageDbg( "WriteDataSetPairwise");
      {
        io::OFStream write;
        std::string data_set_filename( AddExampleOutputPathToFilename( def_constr, "data_set_b"));
        BCL_ExampleMustOpenOutputFile( write, data_set_filename);
        param_constr.WriteDataSetPairwise( write, data_set);
        io::File::CloseClearFStream( write);
        BCL_ExampleCheck( io::File::FilesMatch( data_set_filename, data_set_filename + ".correct"), true);
      }
      // ReadDataSetPairwise
      BCL_MessageDbg( "ReadDataSetPairwise");
      {
        io::IFStream read;
        std::string data_set_filename( AddExampleOutputPathToFilename( def_constr, "data_set_b"));
        BCL_ExampleMustOpenInputFile( read, data_set_filename);
        restraint::HandlerDataSetPairwiseIdentifiers handler;
        handler.ReadDataSetPairwise( read);
        BCL_ExampleCheck( handler.GetScore(), 3.0);
        BCL_ExampleCheck( handler.GetDataSetPairwise().GetSize(), 2);
        def_constr = handler;
      }

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      // test Write and read
      {
        WriteBCLObject( def_constr);
        restraint::HandlerDataSetPairwiseIdentifiers read_handler;
        ReadBCLObject( read_handler);
        BCL_ExampleCheck( read_handler.GetScore(), 3.0);
        BCL_ExampleCheck( read_handler.GetDataSetPairwise().GetSize(), 2);
      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintDataSetPairwiseIdentifiers

  const ExampleClass::EnumType ExampleRestraintHandlerDataSetPairwiseIdentifiers::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintHandlerDataSetPairwiseIdentifiers())
  );

} // namespace bcl
