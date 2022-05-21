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
#include "storage/bcl_storage_symmetric_matrix.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_storage_symmetric_matrix.cpp
  //!
  //! @author teixeipl
  //! @date Aug 18, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleStorageSymmetricMatrix :
    public ExampleInterface
  {
  public:

    ExampleStorageSymmetricMatrix *Clone() const
    {
      return new ExampleStorageSymmetricMatrix( *this);
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
      size_t size( 6);
      storage::SymmetricMatrix< size_t> default_matrix( size);

      // Clone tested further down

      // Fill up the symmetric matrix
      size_t counter( 1);
      for( size_t i_count( 0); i_count < default_matrix.GetSize(); ++i_count)
      {
        for( size_t j_count( 0); j_count <= i_count; ++j_count, ++counter)
        {
          default_matrix( i_count, j_count) = counter;
        }
      }

    /////////////////
    // data access //
    /////////////////

      // Test GetSize()
      size_t size_chk( 6);
      BCL_ExampleIndirectCheck( default_matrix.GetSize(), size_chk, "Check Size");

    ///////////////
    // operators //
    ///////////////

      // Try to access first element
      size_t first_elem_chk( 1);
      BCL_ExampleIndirectCheck( default_matrix( 0, 0), first_elem_chk, "Check First Element");

      // Try to access last element
      size_t last_elem_chk( 21);
      BCL_ExampleIndirectCheck( default_matrix( 5, 5), last_elem_chk, "Check Last Element");

      // Try to access random position
      size_t pos_a_chk( 8);
      BCL_ExampleIndirectCheck( default_matrix( 3, 1), pos_a_chk, "Check Position Element");

      // Try to access position where j > i in i,j
      size_t pos_b_chk( 8);
      BCL_ExampleIndirectCheck( default_matrix( 1, 3), pos_b_chk, "Check Position Element");

      // Compare result of accessing i,j and j,i (should be identical) and Clone()
      util::ShPtr< storage::SymmetricMatrix< size_t> > clone_matrix( default_matrix.Clone());
      BCL_ExampleIndirectCheck( default_matrix( 1, 3), ( *clone_matrix)( 3, 1), "Compare Symmetry");

      // Try to obtain const reference
      const size_t temp_ref_chk( 3);
      const size_t temp_ref( default_matrix( 1, 1));
      BCL_ExampleIndirectCheck( temp_ref, temp_ref_chk, "Check const reference");

      // Check two equal matrices
      BCL_ExampleIndirectCheck( *clone_matrix == default_matrix, true, "Check == operator");

      // Change element at position i,j
      size_t i_index( 3);
      size_t j_index( 4);
      size_t test_value( 201);
      default_matrix( i_index, j_index) = test_value;
      BCL_ExampleIndirectCheck( default_matrix( i_index, j_index), test_value, "Check Change Element");

      // Check two unequal matrices
      BCL_ExampleIndirectCheck( *clone_matrix == default_matrix, false, "Check == operator");

    ////////////////
    // operations //
    ////////////////

      // Test Reset with no argument
      storage::SymmetricMatrix< size_t> copy_matrix( default_matrix);
      copy_matrix.Reset();
      size_t empty_size( 0);
      BCL_ExampleIndirectCheck( copy_matrix.GetSize(), empty_size, "Reset");

      // Test Reset with new size provided as param
      const size_t reset_size( 6);
      copy_matrix.Resize( reset_size);
      BCL_ExampleIndirectCheck( copy_matrix.GetSize(), reset_size, "Reset With Size");

      // Test Fill TODO:

    //////////////////////
    // input and output //
    //////////////////////

      // Test write and read
      BCL_ExampleCheck( TestBCLObjectIOForSymmetry( default_matrix, storage::SymmetricMatrix< size_t>()), true);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleStorageSymmetricMatrix

  const ExampleClass::EnumType ExampleStorageSymmetricMatrix::s_Instance
  (
    GetExamples().AddEnum( ExampleStorageSymmetricMatrix())
  );

} // namespace bcl
