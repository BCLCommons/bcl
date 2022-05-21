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
#include "assemble/bcl_assemble_protein_model_data.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_protein_model_data.cpp
  //!
  //! @author karakam
  //! @date Nov 27, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleProteinModelData :
    public ExampleInterface
  {
  public:

    ExampleAssembleProteinModelData *Clone() const
    {
      return new ExampleAssembleProteinModelData( *this);
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

      // test default constructor
      assemble::ProteinModelData data;

      // test clone constructor
      util::ShPtr< assemble::ProteinModelData> sp_data( data.Clone());

    /////////////////
    // data access //
    /////////////////

      // test GetStaticClassName
      BCL_ExampleCheck( GetStaticClassName< assemble::ProteinModelData>(), "bcl::assemble::ProteinModelData");

    ////////////////
    // operations //
    ////////////////

      // initialize data vector
      linal::Vector3D vector_a( 1.0, 2.0, 3.0);
      linal::Vector3D vector_b( 3.0, 4.0, 5.0);
      linal::Vector3D vector_c( 5.0, 6.0, 7.0);
      util::ShPtr< linal::Vector3D> sp_vector_a( vector_a.Clone());
      util::ShPtr< linal::Vector3D> sp_vector_b( vector_b.Clone());
      util::ShPtr< linal::Vector3D> sp_vector_c( vector_c.Clone());
      const assemble::ProteinModelData::Type vector_key( assemble::ProteinModelData::e_Membrane);

      // insert into the data
      BCL_ExampleCheck( data.Insert( vector_key, sp_vector_a), true);

      // try inserting again which should fail
      BCL_ExampleCheck( data.Insert( vector_key, sp_vector_a), false);

      // now get the data and cast it
      util::ShPtr< linal::Vector3D> sp_vector_a1( data.GetData( vector_key));

      // make sure it is defined and has correct values
      BCL_ExampleCheck( sp_vector_a1.IsDefined(), true);
      BCL_ExampleCheck( sp_vector_a1, sp_vector_a);

      // now replace with sp_vector_b
      BCL_ExampleCheck( data.Replace( vector_key, sp_vector_b), true);

      // get the data back
      util::ShPtr< linal::Vector3D> sp_vector_b1( data.GetData( vector_key));

      // make sure it is defined and has correct values
      BCL_ExampleCheck( sp_vector_b1.IsDefined(), true);
      BCL_ExampleCheck( sp_vector_b1, sp_vector_b);

      // now create a new key that does not exist yet
      const assemble::ProteinModelData::Type vector_key_c( assemble::ProteinModelData::e_LoopDomainLocators);

      // try to get the data with this key
      util::ShPtr< linal::Vector3D> sp_vector_c1( data.GetData( vector_key_c));

      // make sure it is not defined
      BCL_ExampleCheck( sp_vector_c1.IsDefined(), false);

      // also make sure replace also fails
      BCL_ExampleCheck( data.Replace( vector_key_c, sp_vector_c), false);

      // now insert vector_c with vector_key_c
      BCL_ExampleCheck( data.Insert( vector_key_c, sp_vector_c), true);
      util::ShPtr< linal::Vector3D> sp_vector_c2( data.GetData( vector_key_c));

      // try to get the data back and make sure it is correct
      BCL_ExampleCheck( sp_vector_c2.IsDefined(), true);
      BCL_ExampleCheck( sp_vector_c2, sp_vector_c);

    //////////////////////
    // input and output //
    //////////////////////

      // write object
      WriteBCLObject( data);
      // initialize a new data and read it in
      assemble::ProteinModelData data_read;
      ReadBCLObject( data_read);
      // get back the data stored for vector_key and vector_key_c
      util::ShPtr< linal::Vector3D> sp_vector_d1( data.GetData( vector_key));
      util::ShPtr< linal::Vector3D> sp_vector_d2( data.GetData( vector_key_c));
      // make sure both are defined
      BCL_ExampleCheck( sp_vector_d1.IsDefined(), true);
      BCL_ExampleCheck( sp_vector_d2.IsDefined(), true);
      // check the contents
      BCL_ExampleCheck( sp_vector_d1, sp_vector_b);
      BCL_ExampleCheck( sp_vector_d2, sp_vector_c);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleProteinModelData

  const ExampleClass::EnumType ExampleAssembleProteinModelData::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleProteinModelData())
  );

} // namespace bcl
