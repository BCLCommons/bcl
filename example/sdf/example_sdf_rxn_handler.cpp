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
#include "sdf/bcl_sdf_rxn_handler.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_sdf_rxn_handler.cpp
  //!
  //! @author geanesar, combss, mendenjl
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleSdfRXNHandler :
    public ExampleInterface
  {
  public:

    ExampleSdfRXNHandler *Clone() const
    {
      return new ExampleSdfRXNHandler( *this);
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

    /////////////////////////////
    // File read format checks //
    /////////////////////////////

      // Check terminal line is read as a terminal line windows and linux
      BCL_ExampleCheck( sdf::RXNHandler::IsRXNMolDelimiter( "$MOL\r"), true);
      BCL_ExampleCheck( sdf::RXNHandler::IsRXNMolDelimiter( "$MOL"),   true);

      // SDF delimiter is not a RXN delimiter 
      BCL_ExampleCheck( sdf::RXNHandler::IsRXNMolDelimiter( "$$$$"),    false);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // check default constructor
      sdf::RXNHandler handler_default;

      // check constructor with parameters
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Chemistry, "multi_reaction_1.rxn"));
      sdf::RXNHandler handler( read);
      io::File::CloseClearFStream( read);

      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Chemistry, "multi_reaction_2.rxn"));
      storage::List< std::string> rxn_buffer;
      std::string buf;
      while( read)
      {
        std::getline( read, buf);
        rxn_buffer.PushBack( buf);
      }
      io::File::CloseClearFStream( read);

      sdf::RXNHandler handler_second;
      handler_second.ReadFromRXN( rxn_buffer.Begin(), rxn_buffer.End());

      // check clone
      util::ShPtr< sdf::RXNHandler> sp_handler( handler.Clone());
      BCL_ExampleCheck( sp_handler.IsDefined(), 1);

      // close input file stream
      io::File::CloseClearFStream( read);

    /////////////////
    // data access //
    /////////////////

      //
      // Test reaction metadata and read success
      //

      // test GetMdlDescriptionLine
      BCL_ExampleCheck( handler.GetDescription(), "\n  Marvin       060501151315\n");
      BCL_ExampleCheck
      ( 
        handler_second.GetDescription(), 
        "\n  Marvin       010901150951\n  A retro Diels-Alder reaction"
      );

      // test the first reaction was read 
      BCL_ExampleCheck( handler.GetNumberReactants(), size_t( 2));
      BCL_ExampleCheck( handler.GetNumberProducts(), size_t( 2));

      // test that a second reaction was read
      BCL_ExampleCheck( handler_second.GetNumberReactants(), size_t( 1));
      BCL_ExampleCheck( handler_second.GetNumberProducts(), size_t( 2));

      //
      // Test internal MDL handlers
      //

      // Reactant descriptions 
      BCL_ExampleCheck( handler.GetReactantHandlers()( 0).GetDescription(), "\n  Mrv1551806051513152D\n");
      BCL_ExampleCheck( handler.GetReactantHandlers()( 1).GetDescription(), "\n  Mrv1551806051513152D\n");

      // Reactant/product counts and molecule sizes for first reaction
      BCL_ExampleAssert( handler.GetReactantHandlers().GetSize(), size_t( 2));
      BCL_ExampleAssert( handler.GetProductHandlers().GetSize(), size_t( 2));
      BCL_ExampleCheck( handler.GetReactantHandlers()( 0).GetAtomInfo().GetSize(), size_t( 3));
      BCL_ExampleCheck( handler.GetReactantHandlers()( 1).GetAtomInfo().GetSize(), size_t( 2));
      BCL_ExampleCheck( handler.GetProductHandlers()( 0).GetAtomInfo().GetSize(), size_t( 3));
      BCL_ExampleCheck( handler.GetProductHandlers()( 1).GetAtomInfo().GetSize(), size_t( 2));

      // Reactant/product counts/molecule sizes for the second reaction
      BCL_ExampleAssert( handler_second.GetReactantHandlers().GetSize(), size_t( 1));
      BCL_ExampleAssert( handler_second.GetProductHandlers().GetSize(), size_t( 2));
      BCL_ExampleCheck( handler_second.GetReactantHandlers()( 0).GetAtomInfo().GetSize(), size_t( 6));
      BCL_ExampleCheck( handler_second.GetProductHandlers()( 0).GetAtomInfo().GetSize(), size_t( 4));
      BCL_ExampleCheck( handler_second.GetProductHandlers()( 1).GetAtomInfo().GetSize(), size_t( 2));

      //
      // Atom mappings are set correctly
      //

      /*
      // first reaction, first reactant
      BCL_ExampleCheck( handler.GetReactantHandlers()( 0).GetAtomInfo()(0).GetAtomMapping(), 2);
      BCL_ExampleCheck( handler.GetReactantHandlers()( 0).GetAtomInfo()(1).GetAtomMapping(), 1);
      BCL_ExampleCheck( handler.GetReactantHandlers()( 0).GetAtomInfo()(2).GetAtomMapping(), 3);

      // first reaction second reactant
      BCL_ExampleAssert( handler.GetReactantHandlers()( 1).GetAtomInfo()( 0).GetAtomMapping(), 0);
      BCL_ExampleAssert( handler.GetReactantHandlers()( 1).GetAtomInfo()( 1).GetAtomMapping(), 4);

      // first reaction first product
      BCL_ExampleAssert( handler.GetProductHandlers()( 0).GetAtomInfo()( 0).GetAtomMapping(), 1);
      BCL_ExampleAssert( handler.GetProductHandlers()( 0).GetAtomInfo()( 1).GetAtomMapping(), 2);
      BCL_ExampleAssert( handler.GetProductHandlers()( 0).GetAtomInfo()( 2).GetAtomMapping(), 4);

      // first reaction second product
      BCL_ExampleAssert( handler.GetProductHandlers()( 1).GetAtomInfo()( 0).GetAtomMapping(), 0);
      BCL_ExampleAssert( handler.GetProductHandlers()( 1).GetAtomInfo()( 1).GetAtomMapping(), 3);

      // second reaction, only reactant
      BCL_ExampleCheck( handler_second.GetReactantHandlers()( 0).GetAtomInfo()( 0).GetAtomMapping(), 1);
      BCL_ExampleCheck( handler_second.GetReactantHandlers()( 0).GetAtomInfo()( 1).GetAtomMapping(), 2);
      BCL_ExampleCheck( handler_second.GetReactantHandlers()( 0).GetAtomInfo()( 2).GetAtomMapping(), 3);
      BCL_ExampleCheck( handler_second.GetReactantHandlers()( 0).GetAtomInfo()( 3).GetAtomMapping(), 4);
      BCL_ExampleCheck( handler_second.GetReactantHandlers()( 0).GetAtomInfo()( 4).GetAtomMapping(), 6);
      BCL_ExampleCheck( handler_second.GetReactantHandlers()( 0).GetAtomInfo()( 5).GetAtomMapping(), 5);
      
      // second reaction, first product
      BCL_ExampleCheck( handler_second.GetProductHandlers()( 0).GetAtomInfo()( 0).GetAtomMapping(), 1);
      BCL_ExampleCheck( handler_second.GetProductHandlers()( 0).GetAtomInfo()( 1).GetAtomMapping(), 2);
      BCL_ExampleCheck( handler_second.GetProductHandlers()( 0).GetAtomInfo()( 2).GetAtomMapping(), 3);
      BCL_ExampleCheck( handler_second.GetProductHandlers()( 0).GetAtomInfo()( 3).GetAtomMapping(), 4);

      // second reaction, second product
      BCL_ExampleCheck( handler_second.GetProductHandlers()( 1).GetAtomInfo()( 0).GetAtomMapping(), 5);
      BCL_ExampleCheck( handler_second.GetProductHandlers()( 1).GetAtomInfo()( 1).GetAtomMapping(), 6);
      */

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSdfRXNHandler

  const ExampleClass::EnumType ExampleSdfRXNHandler::s_Instance
  (
    GetExamples().AddEnum( ExampleSdfRXNHandler())
  );

} // namespace bcl
